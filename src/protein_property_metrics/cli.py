"""Command-line interface for descriptive protein properties."""

from __future__ import annotations

import argparse
import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
from dataclasses import dataclass, replace
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO

from .nfive import DEFAULT_MODEL_PATH, N_DEGRON_STABILITY_GROUP
from .sequence import DEFAULT_PH, DEFAULT_SEQUENCE_GROUPS, SEQUENCE_GROUPS, describe, normalize_sequence
from .surface import DEFAULT_RADIUS, DEFAULT_SASA_POINTS, describe_structure

BASE_CONTEXT = "base"


@dataclass(frozen=True)
class Record:
    identifier: str
    sequence: str
    structure: Path | None = None
    chain: str | None = None


def _groups(value: str) -> tuple[str, ...]:
    groups = tuple(group.strip().lower() for group in value.split(","))
    allowed = (*SEQUENCE_GROUPS, "sap")
    if any(group not in allowed for group in groups) or len(set(groups)) != len(groups):
        raise argparse.ArgumentTypeError(f"use each group once from {', '.join(allowed)}")
    return groups


def _context(value: str) -> tuple[str, str, str]:
    parts = value.split(":")
    if len(parts) != 3 or not parts[0].strip() or parts[0].strip() == BASE_CONTEXT:
        raise argparse.ArgumentTypeError("use NAME:PREFIX:SUFFIX with a name other than 'base'")
    name, prefix, suffix = (part.strip() for part in parts)
    try:
        return (
            name,
            normalize_sequence(prefix) if prefix else "",
            normalize_sequence(suffix) if suffix else "",
        )
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid residues for context {name!r}: {exc}") from exc


def _positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("expected a positive integer")
    return number


def _positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("expected a finite positive number")
    return number


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="protein-property-metrics", description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--sequence", help="One protein amino acid sequence")
    inputs.add_argument("--fasta", type=Path, help="Protein FASTA file")
    inputs.add_argument("--manifest", type=Path, help="CSV with id, sequence, and optional structure_path and chain")
    parser.add_argument("--structure", type=Path, help="PDB or mmCIF file for a single input record")
    parser.add_argument("--chain", help="Structure chain ID for a single input record")
    parser.add_argument(
        "--groups", type=_groups, default=DEFAULT_SEQUENCE_GROUPS,
        help=f"Comma-separated descriptor groups. The default is {','.join(DEFAULT_SEQUENCE_GROUPS)}. "
        f"Add sap for structures or {N_DEGRON_STABILITY_GROUP} for N-FIVE's predicted PSI. "
        "Lower PSI indicates a more destabilizing N-terminus.",
    )
    parser.add_argument(
        "--nfive-model", type=Path,
        help="Path to the trusted N-FIVE WT model.pkl. The default is the model in the submodule source checkout.",
    )
    parser.add_argument("--ph", type=float, default=DEFAULT_PH, help=f"pH for net charge. The default is {DEFAULT_PH}.")
    parser.add_argument(
        "--radius", type=_positive_float, default=DEFAULT_RADIUS,
        help=f"Static SAP neighborhood radius in angstroms. The default is {DEFAULT_RADIUS}.",
    )
    parser.add_argument(
        "--sasa-points", type=_positive_int, default=DEFAULT_SASA_POINTS,
        help=f"Solvent-accessible surface sample points per atom. The default is {DEFAULT_SASA_POINTS}.",
    )
    parser.add_argument(
        "--context", metavar="NAME:PREFIX:SUFFIX", action="append", type=_context, default=[],
        help="Add a named sequence context with an optional prefix and suffix. May be repeated.",
    )
    parser.add_argument("--workers", type=_positive_int, default=1, help="Worker processes. The default is 1.")
    parser.add_argument("--output", "--csv", required=True, type=Path, help="Output CSV path")
    parser.add_argument("--jsonl", type=Path, help="Optional output JSON Lines path")
    return parser


def _manifest_records(path: Path) -> list[Record]:
    records = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, strict=True)
        fields = reader.fieldnames or []
        if not {"id", "sequence"}.issubset(fields) or len(set(fields)) != len(fields):
            raise ValueError(f"{path} needs unique column names including 'id' and 'sequence'")
        for row in reader:
            identifier = (row["id"] or "").strip()
            if None in row or not identifier or not row["sequence"]:
                raise ValueError(f"{path} line {reader.line_num} must have an ID and sequence in the declared columns")
            structure = (row.get("structure_path") or "").strip()
            chain = (row.get("chain") or "").strip() or None
            structure_path = path.parent / Path(structure).expanduser() if structure else None
            records.append(Record(identifier, row["sequence"], structure_path, chain))
    return records


def _records(args: argparse.Namespace) -> list[Record]:
    if args.sequence is not None:
        records = [Record("sequence", args.sequence)]
    elif args.fasta:
        with args.fasta.open(encoding="utf-8-sig") as handle:
            records = [Record(record.id, str(record.seq)) for record in SeqIO.parse(handle, "fasta")]
        if any(not record.identifier for record in records):
            raise ValueError(f"{args.fasta} has a FASTA record with an empty ID")
    else:
        records = _manifest_records(args.manifest)
    if not records:
        raise ValueError("the input contains no protein records")
    if args.structure is not None or args.chain is not None:
        if len(records) != 1:
            raise ValueError("--structure and --chain require exactly one input record")
        record = records[0]
        records[0] = replace(
            record,
            structure=args.structure or record.structure,
            chain=args.chain if args.chain is not None else record.chain,
        )
        if records[0].structure is None:
            raise ValueError("--chain requires a structure path")
    return records


def _describe(
    record: Record,
    *,
    groups: tuple[str, ...],
    ph: float,
    radius: float,
    sasa_points: int,
    contexts: list[tuple[str, str, str]],
    nfive_model: Path | None,
) -> list[dict[str, str | int | float | None]]:
    sequence_groups = tuple(group for group in groups if group != "sap")
    try:
        sequence = normalize_sequence(record.sequence)
        properties = describe(sequence, groups=sequence_groups, ph=ph, nfive_model=nfive_model)
        base = {"id": record.identifier, "context": BASE_CONTEXT, **properties}
        if "sap" in groups and record.structure is not None:
            surface = describe_structure(
                record.structure,
                chain=record.chain,
                expected_sequence=sequence,
                radius=radius,
                sasa_points=sasa_points,
            )
            base.update(surface)
        rows = [base]
        for name, prefix, suffix in contexts:
            properties = describe(prefix + sequence + suffix, groups=sequence_groups, ph=ph, nfive_model=nfive_model)
            rows.append({"id": record.identifier, "context": name, **properties})
    except (OSError, ValueError) as exc:
        raise ValueError(f"record {record.identifier!r}: {exc}") from exc
    return rows


def _write_outputs(rows: list[dict[str, str | int | float | None]], csv_path: Path, jsonl_path: Path | None) -> None:
    columns = list(dict.fromkeys(column for row in rows for column in row))
    with ExitStack() as stack:
        staged_csv = Path(stack.enter_context(TemporaryDirectory(dir=csv_path.parent))) / "properties.csv"
        with staged_csv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
        if jsonl_path is not None:
            staged_jsonl = Path(stack.enter_context(TemporaryDirectory(dir=jsonl_path.parent))) / "properties.jsonl"
            with staged_jsonl.open("w", encoding="utf-8") as handle:
                for row in rows:
                    handle.write(json.dumps(row, allow_nan=False) + "\n")
            staged_jsonl.replace(jsonl_path)
        staged_csv.replace(csv_path)


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if len({name for name, _, _ in args.context}) != len(args.context):
            raise ValueError("context names must be unique")
        records = _records(args)
        structures = {record.structure for record in records if record.structure is not None}
        if "sap" in args.groups and not structures:
            raise ValueError("the sap group requires at least one structure path")
        csv_path = args.output.resolve()
        jsonl_path = args.jsonl.resolve() if args.jsonl is not None else None
        outputs = [path for path in (csv_path, jsonl_path) if path is not None]
        inputs = {path.resolve() for path in (*structures, args.fasta, args.manifest) if path is not None}
        if N_DEGRON_STABILITY_GROUP in args.groups:
            inputs.add((args.nfive_model or DEFAULT_MODEL_PATH).expanduser().resolve())
        if len(set(outputs)) != len(outputs) or inputs.intersection(outputs):
            raise ValueError("output paths must be distinct from each other and from inputs")
        if any(path.exists() and not path.is_file() for path in outputs):
            raise ValueError("outputs must be regular file paths")

        describe_record = partial(
            _describe, groups=args.groups, ph=args.ph, radius=args.radius, sasa_points=args.sasa_points,
            contexts=args.context, nfive_model=args.nfive_model,
        )
        if args.workers == 1:
            rows = [row for record in records for row in describe_record(record)]
        else:
            with ProcessPoolExecutor(max_workers=min(args.workers, len(records))) as pool:
                rows = [row for batch in pool.map(describe_record, records) for row in batch]
        _write_outputs(rows, csv_path, jsonl_path)
    except (OSError, ValueError, csv.Error) as exc:
        parser.exit(2, f"{parser.prog}: error: {exc}\n")
    return 0
