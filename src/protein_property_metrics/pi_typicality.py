"""Organism-specific typicality of a calculated protein isoelectric point."""

import json
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from functools import lru_cache
from importlib.resources import files
from math import isfinite

from .organisms import OrganismProfile

PI_TYPICALITY_GROUP = "pi-typicality"


@dataclass(frozen=True)
class PIReference:
    """Reference pIs and a fixed, symmetric local density window in pH units."""

    pis: tuple[float, ...]
    half_window: float
    _native_counts: tuple[int, ...] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not self.pis or not all(isfinite(pi) for pi in self.pis):
            raise ValueError("pI reference requires finite pI values")
        if not isfinite(self.half_window) or self.half_window <= 0:
            raise ValueError("pI reference half-window must be positive and finite")
        object.__setattr__(self, "pis", tuple(sorted(self.pis)))
        object.__setattr__(self, "_native_counts", tuple(sorted(self._neighbor_count(pi) for pi in self.pis)))

    def _neighbor_count(self, pi: float) -> int:
        return bisect_right(self.pis, pi + self.half_window) - bisect_left(self.pis, pi - self.half_window)

    def percentile(self, pi: float) -> float:
        """Percent of reference pIs whose local density is at most the query's."""
        if not isfinite(pi):
            raise ValueError("pI typicality requires a finite pI")
        rank = bisect_right(self._native_counts, self._neighbor_count(pi))
        return 100.0 * rank / len(self.pis)


@lru_cache(maxsize=4)
def _load_reference(name: str) -> tuple[PIReference, str, str]:
    resource = files(__package__).joinpath("references").joinpath(name)
    with resource.open(encoding="utf-8") as handle:
        data = json.load(handle)
    # Each row stores the host sequence ID, annotation accession, and pI.
    reference = PIReference(tuple(row[2] for row in data["records"]), data["half_window_ph"])
    return reference, data["proteome_id"], data["uniprot_release"]


def describe_pi_typicality(pi: float, organism: OrganismProfile, *, details: dict | None = None) -> float:
    """Score a pI against the registered cohort and optionally record provenance."""
    if organism.pi_reference is None:
        raise ValueError(f"{organism.name} has no pI reference cohort")
    reference, proteome_id, release = _load_reference(organism.pi_reference)
    if details is not None:
        details["pi_typicality"] = {
            "organism": organism.name,
            "proteome_id": proteome_id,
            "reference_count": len(reference.pis),
            "uniprot_release": release,
            "half_window_ph": reference.half_window,
            "isoelectric_point": pi,
        }
    return reference.percentile(pi)
