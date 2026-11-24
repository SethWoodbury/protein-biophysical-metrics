"""Construct-relative motif counts and optional site annotations."""

from collections import Counter
from collections.abc import Iterable, Sequence
from functools import partial
from math import isfinite
from operator import index
import re

SITE_FLAGS_GROUP = "site-flags"

_DEAMIDATION_SOURCE = "https://pubmed.ncbi.nlm.nih.gov/1649998/"
_CHEMICAL_STRESS_SOURCE = "https://pubmed.ncbi.nlm.nih.gov/33475951/"
_CYCLIZATION_SOURCE = "https://pubmed.ncbi.nlm.nih.gov/21282104/"
_POLYPROLINE_SOURCE = "https://pubmed.ncbi.nlm.nih.gov/23239624/"
_TEV_SOURCE = "https://www.neb.com/en-us/products/p8112-tev-protease"
_WORKFLOW_PROTEASE_SOURCE = (
    "https://b2b.sigmaaldrich.com/SE/en/technical-documents/"
    "technical-article/protein-biology/protein-purification/cleavage-enzymes-fusion-protein-purification"
)
_CLPX_SOURCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC517522/"
_MUA_SOURCE = "https://doi.org/10.1016/j.jmb.2015.03.008"
_N_END_SOURCE = "https://doi.org/10.1038/s41467-026-79349-0"

_CHEMICAL_COUNTS = ("NG", "DG", "NS")
_CHEMICAL_FLAGS = {
    "NG": ("chemical_liability", _DEAMIDATION_SOURCE, "Asn deamidation liability"),
    "DG": ("chemical_liability", _CHEMICAL_STRESS_SOURCE, "Asp isomerization liability"),
    "NS": ("chemical_liability", _DEAMIDATION_SOURCE, "Asn deamidation liability"),
    "DP": ("chemical_annotation", _CHEMICAL_STRESS_SOURCE, "Condition-dependent cleavage liability"),
}
_CHEMICAL_PATTERN = re.compile(r"NG|DG|NS|DP")
_PROLINE_PATTERN = re.compile(r"(?=(PPP|PPG))")
_PROLINE_RUN = re.compile(r"P+")

_COMPLETE_TAGS = (
    ("AANDENYALAA", "ssrA-native", "https://doi.org/10.7554/eLife.61496"),
    ("AANDENYALVA", "ssrA-LVA", "https://doi.org/10.1186/s13036-024-00420-9"),
    ("AANDENYAAAV", "ssrA-AAV", "https://doi.org/10.1186/s13036-024-00420-9"),
    ("AANDENYAASV", "ssrA-ASV", "https://doi.org/10.1186/s13036-024-00420-9"),
)
_SHORT_TAGS = (
    ("RRKKAI", _MUA_SOURCE),
    ("LAA", _CLPX_SOURCE),
    ("VAA", _CLPX_SOURCE),
    ("AA", _CLPX_SOURCE),
)

# A lookahead preserves matches that share recognition residues.
# Offsets count matched residues up to the P1 residue, before the cut.
_PROTEASE_RULES = {
    "tev": (re.compile(r"(?=(ENLYFQ[SGAMCH]))"), 6, _TEV_SOURCE),
    "hrv-3c": (re.compile(r"(?=(LEVLFQGP))"), 6, _WORKFLOW_PROTEASE_SOURCE),
    "thrombin": (re.compile(r"(?=(LVPRGS))"), 4, _WORKFLOW_PROTEASE_SOURCE),
    "enterokinase": (re.compile(r"(?=(DDDDK[^P]))"), 5, _WORKFLOW_PROTEASE_SOURCE),
    "factor-xa": (re.compile(r"(?=(I[ED]GR[^PR]))"), 4, _WORKFLOW_PROTEASE_SOURCE),
}


def normalize_proteases(proteases_used: Iterable[str] | str) -> tuple[str, ...]:
    """Validate the selected workflow enzymes and retain their input order."""
    if isinstance(proteases_used, str):
        proteases_used = proteases_used.split(",") if proteases_used.strip() else ()
    selected = {}
    for name in proteases_used:
        canonical = name.strip().lower().replace(" ", "-") if isinstance(name, str) else None
        if canonical not in _PROTEASE_RULES:
            raise ValueError(f"unknown workflow protease {name!r}. Choose from {', '.join(_PROTEASE_RULES)}")
        selected[canonical] = None
    return tuple(selected)


def validate_positions(positions: Iterable[int], maximum: int, name: str) -> set[int]:
    """Validate 1-based positions for the supplied coordinate bound."""
    result = set()
    for position in positions:
        try:
            resolved = index(position)
        except TypeError as error:
            raise ValueError(f"{name} must contain integer positions") from error
        if isinstance(position, bool) or not 1 <= resolved <= maximum:
            raise ValueError(f"invalid {name} position {position!r}. Expected 1 through {maximum}.")
        result.add(resolved)
    return result


def _site(
    motif_id: str, category: str, matched: str, start: int, provenance: str,
    *, protected: set[int], exposure: tuple[float | None, ...] | None,
    system: str | None = None, protease: str | None = None, cleavage_after: int | None = None,
    intended: bool | None = None, interpretation: str,
) -> dict[str, object]:
    end = start + len(matched) - 1
    return {
        "motif_id": motif_id,
        "bin": category,
        "matched_sequence": matched,
        "start": start,
        "end": end,
        "system": system,
        "protease": protease,
        "cleavage_after": cleavage_after,
        "intended": intended,
        "protected": any(position in protected for position in range(start, end + 1)),
        "sasa_a2": list(exposure[start - 1:end]) if exposure is not None else None,
        "interpretation": interpretation,
        "provenance": provenance,
    }


def describe_sites(
    sequence: str,
    *,
    proteases_used: Iterable[str] | str = (),
    intended_cuts: Iterable[int] = (),
    protected_positions: Iterable[int] = (),
    mature_n_terminus: bool = False,
    free_c_terminus: bool = True,
    residue_sasa_a2: Sequence[float | None] | None = None,
    details: dict | None = None,
) -> dict[str, int | None]:
    """Count motifs on a validated canonical construct and record optional site details.

    Intended cuts are indexed by their P1 residue. Protection and site intervals
    use 1-based sequence coordinates. Optional total residue SASA is supplied in
    sequence order. A true mature_n_terminus declares the first residue exposed.
    """
    length = len(sequence)
    if not length:
        raise ValueError("sequence is empty")
    if not isinstance(mature_n_terminus, bool) or not isinstance(free_c_terminus, bool):
        raise ValueError("terminal declarations must be Boolean")
    selected = normalize_proteases(proteases_used)
    intended_positions = validate_positions(intended_cuts, length - 1, "intended cut")
    protected = validate_positions(protected_positions, length, "protected residue")
    exposure = None
    if residue_sasa_a2 is not None:
        if len(residue_sasa_a2) != length:
            raise ValueError("residue SASA must contain one entry per sequence position")
        try:
            exposure = tuple(None if area is None else float(area) for area in residue_sasa_a2)
        except (TypeError, ValueError) as error:
            raise ValueError("residue SASA must contain finite nonnegative areas or None") from error
        if any(area is not None and (not isfinite(area) or area < 0) for area in exposure):
            raise ValueError("residue SASA must contain finite nonnegative areas or None")

    # These three dimers cannot overlap themselves. str.count therefore counts
    # all starting positions for each motif independently.
    metrics = {f"{motif.lower()}_count": sequence.count(motif) for motif in _CHEMICAL_COUNTS}
    sites = []
    site = partial(_site, protected=protected, exposure=exposure)
    if details is not None:
        for match in _CHEMICAL_PATTERN.finditer(sequence):
            motif = match[0]
            category, source, interpretation = _CHEMICAL_FLAGS[motif]
            sites.append(site(motif, category, motif, match.start() + 1, source,
                              interpretation=interpretation))
        first = sequence[0]
        if mature_n_terminus and first in "QE":
            sites.append(site(f"N-terminal-{first}", "chemical_annotation", first, 1, _CYCLIZATION_SOURCE,
                              interpretation="Possible N-terminal cyclization"))
        if mature_n_terminus and first in "FLWY":
            sites.append(site("ClpS-primary", "n_end_signal", first, 1, _N_END_SOURCE,
                              system="ClpS", interpretation="Classical ClpS-recognized N-end residue"))
        if mature_n_terminus and first in "RK":
            sites.append(site("N-end-secondary", "n_end_signal", first, 1, _N_END_SOURCE,
                              system="LFTR", interpretation="Secondary signal requiring modification"))

    clusters = []
    cluster_count = 0
    cluster_end = 0
    for match in _PROLINE_PATTERN.finditer(sequence):
        start = match.start(1) + 1
        end = match.end(1)
        if start > cluster_end:
            cluster_count += 1
            if details is not None:
                clusters.append({"start": start, "end": end})
        elif details is not None:
            clusters[-1]["end"] = end
        cluster_end = end
        if details is not None:
            motif = match[1]
            sites.append(site(motif, "polyproline_match", motif, start, _POLYPROLINE_SOURCE,
                              interpretation="Translation pausing motif"))
    metrics["polyproline_cluster_count"] = cluster_count

    complete = None
    minor = None
    if free_c_terminus:
        complete = next((tag for tag in _COMPLETE_TAGS if sequence.endswith(tag[0])), None)
        if complete is None:
            minor = next((tag for tag in _SHORT_TAGS if sequence.endswith(tag[0])), None)
    metrics.update(degron_tag_count=int(complete is not None), degron_motif_count=int(minor is not None))
    if details is not None and complete is not None:
        matched, motif_id, source = complete
        sites.append(site(motif_id, "major", matched, length - len(matched) + 1, source,
                          system="ClpXP", interpretation="Complete C-terminal ssrA tag"))
    if details is not None and minor is not None:
        matched, source = minor
        interpretation = (
            "MuA-type signal whose recognition depends on additional contacts and context"
            if matched == "RRKKAI" else "Context-dependent C-terminal ClpX recognition motif"
        )
        sites.append(site(matched, "minor", matched, length - len(matched) + 1, source,
                          system="ClpX", interpretation=interpretation))

    unexpected = set()
    for enzyme in selected:
        pattern, offset, source = _PROTEASE_RULES[enzyme]
        for match in pattern.finditer(sequence):
            cut = match.start(1) + offset
            is_intended = cut in intended_positions
            if not is_intended:
                unexpected.add(cut)
            if details is not None:
                sites.append(site(enzyme, "workflow_protease", match[1], match.start(1) + 1, source,
                                  protease=enzyme, cleavage_after=cut, intended=is_intended,
                                  interpretation="Workflow-protease recognition motif"))
    metrics["protease_cut_site_count"] = len(unexpected) if selected else None

    if details is not None:
        details["site_flags"] = {
            "sites": sites,
            "category_counts": dict(Counter(item["bin"] for item in sites)),
            "polyproline_clusters": clusters,
            "longest_proline_run": max((len(match[0]) for match in _PROLINE_RUN.finditer(sequence)), default=0),
            "proteases_used": list(selected),
            "intended_cuts": sorted(intended_positions),
            "protected_positions": sorted(protected),
            "mature_n_terminus": mature_n_terminus,
            "free_c_terminus": free_c_terminus,
        }
    return metrics
