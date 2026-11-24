"""Calculate descriptive metrics for an unmodified polypeptide sequence."""

from collections.abc import Iterable
from math import log10

from Bio.Data.IUPACData import protein_letters
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils.ProtParam import ProteinAnalysis

DEFAULT_PH = 7.4
SEQUENCE_GROUPS = ("composition", "charge", "instability")
_SEQUENCE_SYMBOLS = frozenset(protein_letters + protein_letters.lower())


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
    groups: Iterable[str] = SEQUENCE_GROUPS,
    ph: float = DEFAULT_PH,
) -> dict[str, int | float]:
    """Calculate the selected descriptor groups in canonical group order.

    Composition includes average mass in daltons, Kyte–Doolittle GRAVY,
    cysteine content, and the aromatic residue fraction (F, W, and Y).
    Charge uses Biopython's Bjellqvist pKa model with free termini. The
    supplied pH must be finite and within [0, 14]. The Guruprasad instability
    index is a sequence-derived score. Its formula gives zero for a single
    residue because the sequence contains no adjacent residue pairs.
    """
    sequence = normalize_sequence(sequence)
    selected = set(groups)
    unknown = selected.difference(SEQUENCE_GROUPS)
    if unknown:
        raise ValueError(f"unknown sequence descriptor groups ({', '.join(sorted(unknown))})")
    ph = float(ph)
    if not 0.0 <= ph <= 14.0:
        raise ValueError("pH must be finite and within [0, 14]")

    analysis = ProteinAnalysis(sequence)
    metrics: dict[str, int | float] = {}
    if "composition" in selected:
        length = len(sequence)
        cysteine_count = sequence.count("C")
        metrics.update(
            length=length,
            molecular_weight_da=analysis.molecular_weight(),
            gravy=analysis.gravy(),
            cysteine_count=cysteine_count,
            cysteine_fraction=cysteine_count / length,
            aromaticity=analysis.aromaticity(),
        )
    if "charge" in selected:
        charge_model = IsoelectricPoint(sequence, analysis.count_amino_acids())
        net_charge = charge_model.charge_at_pH(ph)
        # Extend the pI search by one decade beyond the maximum residue count.
        # This brackets the zero of the charge model for highly biased sequences.
        margin = log10(len(sequence) + 1) + 1
        metrics.update(
            isoelectric_point=charge_model.pi(pH=7.0, min_=-margin, max_=14.0 + margin),
            net_charge=net_charge,
            absolute_net_charge=abs(net_charge),
        )
    if "instability" in selected:
        metrics["instability_index"] = analysis.instability_index()
    return metrics
