"""Organism-specific reference data, motif rules, and scan defaults."""

import re
from collections.abc import Mapping
from dataclasses import dataclass, field

DEFAULT_ORGANISM = "BL21_DE3_e_coli"
WORKFLOW_PROTEASE = "workflow_protease"
CONTEXT_DEPENDENT_PROTEASE = "context_dependent_protease"


@dataclass(frozen=True)
class ProteaseRule:
    """An overlapping recognition pattern and its candidate cut offset."""

    motif: str
    cut_offset: int
    source: str
    category: str = WORKFLOW_PROTEASE
    interpretation: str = "Workflow-protease recognition motif"
    pattern: re.Pattern[str] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "pattern", re.compile(f"(?=({self.motif}))"))


@dataclass(frozen=True)
class DegronRule:
    """A terminal motif or set of candidate mature N-terminal residues."""

    motif: str
    motif_id: str
    source: str
    system: str
    interpretation: str


@dataclass(frozen=True)
class OrganismProfile:
    """Settings for one host and its scoring references."""

    name: str
    proteases: Mapping[str, ProteaseRule] = field(default_factory=dict)
    default_proteases: tuple[str, ...] = ()
    complete_degrons: tuple[DegronRule, ...] = ()
    short_degrons: tuple[DegronRule, ...] = ()
    n_end_degrons: tuple[DegronRule, ...] = ()
    pi_reference: str | None = None
    nfive_supported: bool = False


_TEV_SOURCE = "https://www.neb.com/en-us/products/p8112-tev-protease"
_WORKFLOW_PROTEASE_SOURCE = (
    "https://b2b.sigmaaldrich.com/SE/en/technical-documents/"
    "technical-article/protein-biology/protein-purification/cleavage-enzymes-fusion-protein-purification"
)
_CLPX_SOURCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC517522/"
_N_END_SOURCE = "https://doi.org/10.1038/s41467-026-79349-0"

# Offset counts residues through P1, immediately before the candidate bond.
WORKFLOW_PROTEASE_RULES = {
    "tev": ProteaseRule("ENLYFQ[SGAMCH]", 6, _TEV_SOURCE),
    "hrv-3c": ProteaseRule("LEVLFQGP", 6, _WORKFLOW_PROTEASE_SOURCE),
    "thrombin": ProteaseRule("LVPRGS", 4, _WORKFLOW_PROTEASE_SOURCE),
    "enterokinase": ProteaseRule("DDDDK[^P]", 5, _WORKFLOW_PROTEASE_SOURCE),
    "factor-xa": ProteaseRule("I[ED]GR[^PR]", 4, _WORKFLOW_PROTEASE_SOURCE),
}

_BL21_PROTEASE_RULES = {
    **WORKFLOW_PROTEASE_RULES,
    # Ordinary BL21(DE3) is OmpT-deficient. These are context candidates.
    "ompt": ProteaseRule(
        "[KR][KR]", 1, "https://pmc.ncbi.nlm.nih.gov/articles/PMC516829/",
        category=CONTEXT_DEPENDENT_PROTEASE,
        interpretation="Potential dibasic site when active OmpT is present",
    ),
}


_BL21_DE3 = OrganismProfile(
    name=DEFAULT_ORGANISM,
    proteases=_BL21_PROTEASE_RULES,
    default_proteases=tuple(_BL21_PROTEASE_RULES),
    complete_degrons=tuple(
        DegronRule(motif, name, source, "ClpXP", "Complete C-terminal ssrA tag")
        for motif, name, source in (
            ("AANDENYALAA", "ssrA-native", "https://doi.org/10.7554/eLife.61496"),
            ("AANDENYALVA", "ssrA-LVA", "https://doi.org/10.1186/s13036-024-00420-9"),
            ("AANDENYAAAV", "ssrA-AAV", "https://doi.org/10.1186/s13036-024-00420-9"),
            ("AANDENYAASV", "ssrA-ASV", "https://doi.org/10.1186/s13036-024-00420-9"),
        )
    ),
    short_degrons=(
        DegronRule(
            "RRKKAI", "RRKKAI", "https://doi.org/10.1016/j.jmb.2015.03.008", "ClpX",
            "MuA-type signal whose recognition depends on additional contacts and context",
        ),
        *(DegronRule(motif, motif, _CLPX_SOURCE, "ClpX", "Context-dependent C-terminal ClpX recognition motif")
          for motif in ("LAA", "VAA", "AA")),
    ),
    n_end_degrons=(
        DegronRule("FLWY", "ClpS-primary", _N_END_SOURCE, "ClpS", "Classical ClpS-recognized N-end residue"),
        DegronRule("RK", "N-end-secondary", _N_END_SOURCE, "LFTR", "Secondary signal requiring modification"),
    ),
    pi_reference="BL21_DE3_e_coli.json",
    nfive_supported=True,
)

ORGANISMS = {_BL21_DE3.name: _BL21_DE3}


def get_organism(organism: str | OrganismProfile = DEFAULT_ORGANISM) -> OrganismProfile:
    """Resolve a registered organism ID or accept a supplied profile."""
    if isinstance(organism, OrganismProfile):
        return organism
    if isinstance(organism, str):
        if organism in ORGANISMS:
            return ORGANISMS[organism]
        for name, profile in ORGANISMS.items():
            if name.casefold() == organism.casefold():
                return profile
    raise ValueError(f"unknown organism {organism!r}. Choose from {', '.join(ORGANISMS)}.")
