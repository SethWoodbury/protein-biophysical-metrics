"""Charge compatibility at the ends of supplied alpha-helical assignments."""

import re

HELIX_GROUP = "helix"
_END_WIDTH = 4
_HELIX_RUN = re.compile(r"H+")
_ASSIGNMENT_CODES = frozenset("HBEGIPTSC.- ")
_CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}


def describe_helix(
    sequence: str,
    secondary_structure: str | None = None,
    *,
    details: dict | None = None,
) -> dict[str, float | None]:
    """Compare the last and first four residues of supplied DSSP H runs.

    Sequence is a normalized canonical chain. Supply one DSSP code per residue
    with C, period, hyphen, or space also accepted for coil. Histidine and other
    residues outside D/E/K/R contribute zero to this dimensionless proxy.
    """
    score = None
    reason = "assignment_unavailable"
    helices = []
    if secondary_structure is not None:
        if not isinstance(secondary_structure, str) or len(secondary_structure) != len(sequence):
            raise ValueError("secondary structure must be a string of the same length as the sequence")
        assignment = secondary_structure.upper()
        unsupported = set(assignment) - _ASSIGNMENT_CODES
        if unsupported:
            raise ValueError(f"secondary structure contains unsupported codes {sorted(unsupported)!r}")

        n_total = c_total = eligible = 0
        for match in _HELIX_RUN.finditer(assignment):
            start, stop = match.span()
            if stop - start < 2 * _END_WIDTH:
                continue
            n_charges = [_CHARGE.get(aa, 0) for aa in sequence[start:start + _END_WIDTH]]
            c_charges = [_CHARGE.get(aa, 0) for aa in sequence[stop - _END_WIDTH:stop]]
            n_total += sum(n_charges)
            c_total += sum(c_charges)
            eligible += 1
            if details is not None:
                helices.append({
                    "start": start + 1,
                    "end": stop,
                    "n_end_charges": n_charges,
                    "c_end_charges": c_charges,
                })
        reason = None if eligible else "no_eligible_helices"
        if eligible:
            score = (c_total - n_total) / (eligible * _END_WIDTH)

    if details is not None:
        details["helix"] = {"reason": reason, "eligible_helices": helices}
    return {"helix_end_charge_compatibility": score}
