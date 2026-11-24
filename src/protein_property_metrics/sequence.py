"""Calculate descriptive metrics for an unmodified polypeptide sequence."""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from math import fsum, log10
from pathlib import Path

from Bio.Data.IUPACData import protein_letters
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .aggrescan import AGGRESCAN_GROUP, describe_aggrescan
from .helix import HELIX_GROUP, describe_helix
from .nfive import N_DEGRON_STABILITY_GROUP, describe_n_terminus
from .organisms import DEFAULT_ORGANISM, OrganismProfile, get_organism
from .pi_typicality import PI_TYPICALITY_GROUP, describe_pi_typicality
from .site_flags import SITE_FLAGS_GROUP, describe_sites, validate_positions

DEFAULT_PH = 7.4
DEFAULT_SEQUENCE_GROUPS = ("composition", "charge", "instability")
SEQUENCE_GROUPS = (
    *DEFAULT_SEQUENCE_GROUPS, "optics", N_DEGRON_STABILITY_GROUP,
    "swi", SITE_FLAGS_GROUP, AGGRESCAN_GROUP, HELIX_GROUP, PI_TYPICALITY_GROUP,
)
_SEQUENCE_SYMBOLS = frozenset(protein_letters + protein_letters.lower())

# Released three-decimal SWI weights from Bhandari, Lim and Gardner (2020),
# doi:10.1093/bioinformatics/btaa578. Independently implement the mean.
# Source: Gardner-BinfLab/TISIGNER-ReactJS, functions.py, hydrop_flex_swi column 3.
# SHA256: bd1808337c372d2dcd2721fc28c2d25fd4c85b8e23ba52361d4135d7cfd9b48c.
_SWI_WEIGHTS = {
    "A": 0.836, "C": 0.521, "D": 0.908, "E": 0.988, "F": 0.585,
    "G": 0.800, "H": 0.895, "I": 0.678, "K": 0.927, "L": 0.655,
    "M": 0.630, "N": 0.860, "P": 0.824, "Q": 0.789, "R": 0.771,
    "S": 0.744, "T": 0.810, "V": 0.736, "W": 0.637, "Y": 0.611,
}

MetricValue = str | int | float | None


@dataclass(frozen=True)
class SiteAnnotations:
    """Declarations in the supplied sequence's 1-based coordinates."""

    intended_cuts: tuple[int, ...] = ()
    protected_positions: tuple[int, ...] = ()
    mature_n_terminus: bool = False
    free_c_terminus: bool = True


def normalize_sequence(sequence: str) -> str:
    """Remove formatting whitespace and uppercase the 20 canonical residue codes."""
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    residues = "".join(sequence.split())
    if not residues:
        raise ValueError("sequence is empty")
    unsupported = set(residues) - _SEQUENCE_SYMBOLS
    if unsupported:
        symbols = ", ".join(repr(symbol) for symbol in sorted(unsupported))
        raise ValueError(f"sequence contains noncanonical residue codes ({symbols})")
    return residues.upper()


def describe(
    sequence: str,
    *,
    groups: Iterable[str] = DEFAULT_SEQUENCE_GROUPS,
    ph: float = DEFAULT_PH,
    organism: str | OrganismProfile = DEFAULT_ORGANISM,
    nfive_model: str | Path | None = None,
    nfive_lookup: str | Path | None = None,
    proteases_used: Iterable[str] | str | None = None,
    intended_cuts: Iterable[int] = (),
    protected_positions: Iterable[int] = (),
    mature_n_terminus: bool = False,
    free_c_terminus: bool = True,
    secondary_structure: str | None = None,
    residue_sasa_a2: Sequence[float | None] | None = None,
    details: dict | None = None,
) -> dict[str, MetricValue]:
    """Calculate the selected descriptor groups in canonical group order.

    Composition includes average mass in daltons, Kyte–Doolittle GRAVY,
    cysteine and methionine counts, aromatic residue fraction (F, W, and Y), and the
    Ikai aliphatic index.
    Charge uses Biopython's Bjellqvist pKa model with free termini. The
    supplied pH must be finite and within [0, 14]. The Guruprasad instability
    index is a sequence-derived score. Its formula gives zero for a single
    residue because the sequence contains no adjacent residue pairs.
    Optics reports the sequence-estimated 280 nm extinction coefficients
    for the fully reduced and maximum disulfide-paired cases.
    N-degron-stability-index reports the five-residue window and predicted
    Protein Stability Index (PSI). Lower PSI indicates a more destabilizing
    N-terminus. MSG starts also get a window and score after Met excision.
    Missing windows return None. The optional nfive_model is the path to
    the trusted upstream model, which defaults to the model in the checkout.
    Supply nfive_lookup to use a matching, precomputed fixed-MSG lookup.

    The organism profile supplies the pI reference cohort, degron rules, and
    protease settings. The pI-typicality group gives a local-density percentile.
    SWI is the mean of the original residue solubility weights. Site flags
    include chemical, proline and protease motifs. A None proteases_used uses
    the profile defaults. Pass an empty iterable to disable protease scanning.
    Intended cuts identify the residue before each bond, and protected
    positions identify residues. Both use 1-based coordinates. Declare
    mature_n_terminus for an exposed mature N-end and set free_c_terminus
    to False for a blocked C-end. Optional residue_sasa_a2 is total SASA
    in sequence order. Helix compatibility uses supplied per-residue DSSP
    labels. Pass a dict as details to collect structured annotations.
    """
    annotations = SiteAnnotations(tuple(intended_cuts), tuple(protected_positions), mature_n_terminus, free_c_terminus)
    return _describe_normalized(
        normalize_sequence(sequence), groups=groups, ph=ph, organism=organism,
        nfive_model=nfive_model, nfive_lookup=nfive_lookup,
        proteases_used=proteases_used, annotations=annotations, secondary_structure=secondary_structure,
        residue_sasa_a2=residue_sasa_a2, details=details,
    )


def _describe_normalized(
    sequence: str,
    *,
    groups: Iterable[str],
    ph: float,
    organism: str | OrganismProfile = DEFAULT_ORGANISM,
    nfive_model: str | Path | None = None,
    nfive_lookup: str | Path | None = None,
    nfive_scores: Mapping[str, float] | None = None,
    proteases_used: Iterable[str] | str | None = None,
    annotations: SiteAnnotations = SiteAnnotations(),
    secondary_structure: str | None = None,
    residue_sasa_a2: Sequence[float | None] | None = None,
    details: dict | None = None,
) -> dict[str, MetricValue]:
    selected = {group.lower() for group in groups}
    unknown = selected.difference(SEQUENCE_GROUPS)
    if unknown:
        raise ValueError(f"unknown sequence descriptor groups ({', '.join(sorted(unknown))})")
    ph = float(ph)
    if not 0.0 <= ph <= 14.0:
        raise ValueError("pH must be finite and within [0, 14]")
    profile = get_organism(organism)

    analysis = ProteinAnalysis(sequence)
    counts = analysis.count_amino_acids() if selected & {"composition", "charge", "swi", PI_TYPICALITY_GROUP} else {}
    metrics: dict[str, MetricValue] = {}
    if "composition" in selected:
        length = len(sequence)
        cysteine_count = counts["C"]
        metrics.update(
            length=length,
            molecular_weight_da=analysis.molecular_weight(),
            gravy=analysis.gravy(),
            cysteine_count=cysteine_count,
            cysteine_fraction=cysteine_count / length,
            methionine_count=counts["M"],
            aromaticity=analysis.aromaticity(),
            aliphatic_index=100 * (counts["A"] + 2.9 * counts["V"] + 3.9 * (counts["I"] + counts["L"])) / length,
        )
    if selected & {"charge", PI_TYPICALITY_GROUP}:
        charge_model = IsoelectricPoint(sequence, counts)
        # Extend the pI search by one decade beyond the maximum residue count.
        # This brackets the zero of the charge model for highly biased sequences.
        margin = log10(len(sequence) + 1) + 1
        pi = charge_model.pi(pH=7.0, min_=-margin, max_=14.0 + margin)
        if "charge" in selected:
            net_charge = charge_model.charge_at_pH(ph)
            metrics.update(isoelectric_point=pi, net_charge=net_charge, absolute_net_charge=abs(net_charge))
        if PI_TYPICALITY_GROUP in selected:
            metrics["pI_typicality_percentile"] = describe_pi_typicality(pi, profile, details=details)
    if "instability" in selected:
        metrics["instability_index"] = analysis.instability_index()
    if "optics" in selected:
        reduced, maximum_paired = analysis.molar_extinction_coefficient()
        metrics["extinction_280_reduced_m_inv_cm_inv"] = reduced
        metrics["extinction_280_max_paired_m_inv_cm_inv"] = maximum_paired
    if N_DEGRON_STABILITY_GROUP in selected:
        if not profile.nfive_supported:
            raise ValueError(f"{profile.name} has no registered N-FIVE model")
        metrics.update(describe_n_terminus(
            sequence, model_path=nfive_model, lookup_path=nfive_lookup, prepared_scores=nfive_scores,
        ))
    if "swi" in selected:
        metrics["swi_solubility_index"] = fsum(counts[aa] * weight for aa, weight in _SWI_WEIGHTS.items()) / len(sequence)
    if SITE_FLAGS_GROUP in selected:
        metrics.update(describe_sites(
            sequence, organism=profile, proteases_used=proteases_used, intended_cuts=annotations.intended_cuts,
            protected_positions=annotations.protected_positions, mature_n_terminus=annotations.mature_n_terminus,
            free_c_terminus=annotations.free_c_terminus, residue_sasa_a2=residue_sasa_a2, details=details,
        ))
    if AGGRESCAN_GROUP in selected:
        metrics.update(describe_aggrescan(sequence, details=details))
    if HELIX_GROUP in selected:
        metrics.update(describe_helix(sequence, secondary_structure, details=details))
    if details is not None:
        protected = validate_positions(annotations.protected_positions, len(sequence), "protected residue")
        if protected:
            details["protected_positions"] = sorted(protected)
    return metrics
