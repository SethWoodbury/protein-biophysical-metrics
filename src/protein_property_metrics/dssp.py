"""Optional mkdssp secondary-structure assignment for one validated chain."""

import subprocess
from io import StringIO
from math import isfinite
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import MMCIFIO
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

DEFAULT_DSSP_EXECUTABLE = "mkdssp"
DEFAULT_DSSP_TIMEOUT = 30.0
_SUMMARY = "_dssp_struct_summary."


def _assignments(output: str, sequence: str) -> str:
    try:
        summary = MMCIF2Dict(StringIO(output))
        chains = summary[_SUMMARY + "label_asym_id"]
        rows = sorted(
            (int(position), protein_letters_3to1[residue], code)
            for _, position, residue, code in zip(
                chains, summary[_SUMMARY + "label_seq_id"], summary[_SUMMARY + "label_comp_id"],
                summary[_SUMMARY + "secondary_structure"], strict=True,
            )
        )
    except (KeyError, ValueError) as exc:
        raise ValueError("DSSP output has no usable per-residue mmCIF assignments in _dssp_struct_summary") from exc
    positions = [position for position, _, _ in rows]
    if len(set(chains)) != 1 or positions != list(range(1, len(sequence) + 1)):
        raise ValueError("DSSP assignments do not cover every residue of the selected chain exactly once")
    if "".join(residue for _, residue, _ in rows) != sequence:
        raise ValueError("DSSP assignments do not match the selected chain sequence")
    return "".join(code for _, _, code in rows)


def assign_dssp(
    protein: Chain, sequence: str, *, executable: str | Path, timeout: float = DEFAULT_DSSP_TIMEOUT,
) -> str:
    """Run mkdssp on an isolated chain and verify every assignment against its sequence."""
    if not isfinite(timeout) or timeout <= 0:
        raise ValueError("DSSP timeout must be finite and positive")
    command = str(executable)
    with TemporaryDirectory(prefix="protein-metrics-dssp-") as temp:
        # The exported chain has sequential label_seq_id values. Using mmCIF
        # supports generated PDBs without crystallographic headers and avoids
        # the classic DSSP format's one-character chain ID limit.
        path = Path(temp) / "chain.cif"
        writer = MMCIFIO()
        writer.set_structure(protein)
        writer.save(str(path))
        try:
            result = subprocess.run(
                [command, "--output-format=mmcif", str(path)],
                capture_output=True, text=True, check=True, timeout=timeout,
            )
        except FileNotFoundError as exc:
            raise ValueError(f"DSSP executable {command!r} was not found") from exc
        except subprocess.TimeoutExpired as exc:
            raise ValueError(f"DSSP exceeded the {timeout:g}-second timeout") from exc
        except subprocess.CalledProcessError as exc:
            diagnostic = (exc.stderr or "").strip() or f"exit status {exc.returncode}"
            raise ValueError(f"DSSP failed for the selected chain. {diagnostic}") from exc
    return _assignments(result.stdout, sequence)
