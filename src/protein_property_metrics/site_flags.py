"""Construct-relative motif counts and optional site annotations."""

from collections import Counter
from collections.abc import Iterable, Sequence
from functools import partial
from math import isfinite
from operator import index
import re

from .organisms import (
    CONTEXT_DEPENDENT_PROTEASE,
    DEFAULT_ORGANISM,
    WORKFLOW_PROTEASE,
    OrganismProfile,
    get_organism,
)

SITE_FLAGS_GROUP = "site-flags"

_DEAMIDATION_SOURCE = "https://pubmed.ncbi.nlm.nih.gov/1649998/"
_CHEMICAL_STRESS_SOURCE = "https://pubmed.ncbi.nlm.nih.gov/33475951/"
_CYCLIZATION_SOURCE = "https://pubmed.ncbi.nlm.nih.gov/21282104/"
_POLYPROLINE_SOURCE = "https://pubmed.ncbi.nlm.nih.gov/23239624/"

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

def normalize_proteases(
    proteases_used: Iterable[str] | str | None = None,
    *,
    organism: str | OrganismProfile = DEFAULT_ORGANISM,
) -> tuple[str, ...]:
    """Resolve host defaults or an explicit ordered enzyme selection."""
    profile = get_organism(organism)
    if proteases_used is None:
        proteases_used = profile.default_proteases
    if isinstance(proteases_used, str):
        proteases_used = proteases_used.split(",") if proteases_used.strip() else ()
    selected = {}
    for name in proteases_used:
        canonical = name.strip().lower().replace(" ", "-") if isinstance(name, str) else None
        if canonical != "none" and canonical not in profile.proteases:
            choices = ", ".join((*profile.proteases, "none"))
            raise ValueError(f"unknown protease {name!r} for {profile.name}. Choose from {choices}.")
        selected[canonical] = None
    if "none" in selected:
        if len(selected) > 1:
            raise ValueError("use 'none' alone to disable protease scanning")
        return ()
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
    proteases_used: Iterable[str] | str | None = None,
    organism: str | OrganismProfile = DEFAULT_ORGANISM,
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
    Omitted proteases use the host defaults. An empty selection disables scans.
    Workflow recognition motifs and context-dependent candidates are counted separately.
    """
    length = len(sequence)
    if not length:
        raise ValueError("sequence is empty")
    if not isinstance(mature_n_terminus, bool) or not isinstance(free_c_terminus, bool):
        raise ValueError("terminal declarations must be Boolean")
    profile = get_organism(organism)
    selected = normalize_proteases(proteases_used, organism=profile)
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
        if mature_n_terminus:
            for rule in profile.n_end_degrons:
                if first in rule.motif:
                    sites.append(site(rule.motif_id, "n_end_signal", first, 1, rule.source,
                                      system=rule.system, interpretation=rule.interpretation))

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
        complete = next((rule for rule in profile.complete_degrons if sequence.endswith(rule.motif)), None)
        if complete is None:
            minor = next((rule for rule in profile.short_degrons if sequence.endswith(rule.motif)), None)
    metrics.update(degron_tag_count=int(complete is not None), degron_motif_count=int(minor is not None))
    terminal_rule = complete or minor
    if details is not None and terminal_rule is not None:
        sites.append(site(
            terminal_rule.motif_id, "major" if complete else "minor", terminal_rule.motif,
            length - len(terminal_rule.motif) + 1, terminal_rule.source,
            system=terminal_rule.system, interpretation=terminal_rule.interpretation,
        ))

    unexpected: dict[str, set[int]] = {}
    for enzyme in selected:
        rule = profile.proteases[enzyme]
        bonds = unexpected.setdefault(rule.category, set())
        for match in rule.pattern.finditer(sequence):
            cut = match.start(1) + rule.cut_offset
            is_intended = cut in intended_positions
            if not is_intended:
                bonds.add(cut)
            if details is not None:
                sites.append(site(enzyme, rule.category, match[1], match.start(1) + 1, rule.source,
                                  protease=enzyme, cleavage_after=cut, intended=is_intended,
                                  interpretation=rule.interpretation))
    metrics["protease_cut_site_count"] = len(unexpected[WORKFLOW_PROTEASE]) if WORKFLOW_PROTEASE in unexpected else None
    if CONTEXT_DEPENDENT_PROTEASE in unexpected:
        metrics["context_dependent_protease_site_count"] = len(unexpected[CONTEXT_DEPENDENT_PROTEASE])

    if details is not None:
        details["site_flags"] = {
            "organism": profile.name,
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
