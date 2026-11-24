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
from itertools import repeat
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO

from .nfive import (
    DEFAULT_MODEL_PATH, MET_EXCISION_PREFIX, N_DEGRON_STABILITY_GROUP, nfive_windows, prepare_nfive_scores,
)
from .organisms import DEFAULT_ORGANISM, ORGANISMS, get_organism
from .pi_typicality import PI_TYPICALITY_GROUP
from .sequence import (
    DEFAULT_PH, DEFAULT_SEQUENCE_GROUPS, SEQUENCE_GROUPS, SiteAnnotations, _describe_normalized, normalize_sequence,
)
from .site_flags import SITE_FLAGS_GROUP, normalize_proteases, validate_positions
from .surface import DEFAULT_RADIUS, DEFAULT_SASA_POINTS, SURFACE_CHEMISTRY_FIELDS, SURFACE_GROUPS, describe_structure

BASE_CONTEXT = "base"


@dataclass(frozen=True)
class Record:
    identifier: str
    sequence: str
    structure: Path | None = None
    chain: str | None = None
    annotations: SiteAnnotations = SiteAnnotations()
    secondary_structure: str | None = None
    context: str = BASE_CONTEXT


@dataclass(frozen=True)
class Context:
    name: str
    prefix: str
    suffix: str
    after_msg: bool = False


def _groups(value: str) -> tuple[str, ...]:
    groups = tuple(group.strip().lower() for group in value.split(","))
    allowed = (*SEQUENCE_GROUPS, *SURFACE_GROUPS)
    if any(group not in allowed for group in groups) or len(set(groups)) != len(groups):
        raise argparse.ArgumentTypeError(f"use each group once from {', '.join(allowed)}")
    return groups


def _context(value: str) -> Context:
    parts = value.split(":")
    if len(parts) != 3 or not parts[0].strip() or parts[0].strip() == BASE_CONTEXT:
        raise argparse.ArgumentTypeError("use NAME:PREFIX:SUFFIX with a name other than 'base'")
    name, prefix, suffix = (part.strip() for part in parts)
    insertion_marker = f"{MET_EXCISION_PREFIX}+"
    after_msg = prefix.upper().startswith(insertion_marker)
    if after_msg:
        prefix = prefix[len(insertion_marker):]
    try:
        return Context(
            name,
            normalize_sequence(prefix) if prefix else "",
            normalize_sequence(suffix) if suffix else "",
            after_msg,
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


def _positions(value: str) -> tuple[int, ...]:
    try:
        positions = tuple(int(part) for part in value.split(",")) if value.strip() else ()
    except ValueError as exc:
        raise ValueError("expected comma-separated 1-based integer positions") from exc
    if any(position < 1 for position in positions):
        raise ValueError("positions must be positive 1-based integers")
    return positions


def _boolean(value: str | None, *, default: bool) -> bool:
    normalized = (value or "").strip().lower()
    if not normalized:
        return default
    if normalized in {"true", "yes", "1"}:
        return True
    if normalized in {"false", "no", "0"}:
        return False
    raise ValueError(f"expected a Boolean, got {value!r}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="protein-property-metrics", description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--sequence", help="One protein amino acid sequence")
    inputs.add_argument("--fasta", type=Path, help="Protein FASTA file")
    inputs.add_argument("--manifest", type=Path, help="CSV with id, sequence, and optional structure_path and chain")
    parser.add_argument("--structure", type=Path, help="PDB or mmCIF file for a single input record")
    parser.add_argument("--chain", help="Structure chain ID for a single input record")
    parser.add_argument(
        "--organism", default=DEFAULT_ORGANISM,
        help=f"Registered host profile. Available profiles are {', '.join(ORGANISMS)}. Default is {DEFAULT_ORGANISM}.",
    )
    parser.add_argument(
        "--groups", type=_groups, default=DEFAULT_SEQUENCE_GROUPS,
        help=f"Comma-separated groups from {', '.join((*SEQUENCE_GROUPS, *SURFACE_GROUPS))}. "
        f"The default is {','.join(DEFAULT_SEQUENCE_GROUPS)}. "
        "Lower N-FIVE PSI indicates a more destabilizing N-terminus.",
    )
    parser.add_argument(
        "--nfive-model", type=Path,
        help="Path to the trusted N-FIVE WT model.pkl. The default is the model in the submodule source checkout.",
    )
    parser.add_argument("--nfive-lookup", type=Path, help="Optional precomputed fixed-MSG N-FIVE lookup (.npz).")
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
        help="Add a named context. Use MSG+INSERT as the prefix to insert after an existing MSG. May be repeated.",
    )
    parser.add_argument(
        "--proteases-used",
        help="Override the host's protease list for site-flags with comma-separated names or 'none'.",
    )
    parser.add_argument("--intended-cuts", help="Comma-separated 1-based positions immediately before intended cuts.")
    parser.add_argument("--protected-positions", help="Comma-separated 1-based protected residue positions.")
    parser.add_argument(
        "--mature-n-terminus", action="store_true", default=None, help="Declare the supplied exposed mature N-end.",
    )
    parser.add_argument(
        "--blocked-c-terminus", action="store_false", dest="free_c_terminus", default=None,
        help="Declare a blocked C-end. The default is free.",
    )
    parser.add_argument("--secondary-structure", help="Supplied per-residue DSSP assignment string for the input.")
    parser.add_argument("--workers", type=_positive_int, default=1, help="Worker processes. The default is 1.")
    parser.add_argument("--output", "--csv", required=True, type=Path, help="Output CSV path")
    parser.add_argument("--jsonl", type=Path, help="Optional output JSON Lines path")
    parser.add_argument("--details", action="store_true", help="Include structured annotations in --jsonl.")
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
            try:
                annotations = SiteAnnotations(
                    intended_cuts=_positions(row.get("intended_cuts") or ""),
                    protected_positions=_positions(row.get("protected_positions") or ""),
                    mature_n_terminus=_boolean(row.get("mature_n_terminus"), default=False),
                    free_c_terminus=_boolean(row.get("free_c_terminus"), default=True),
                )
            except ValueError as exc:
                raise ValueError(f"{path} line {reader.line_num}: {exc}") from exc
            records.append(Record(
                identifier, row["sequence"], structure_path, chain,
                annotations=annotations, secondary_structure=row.get("secondary_structure") or None,
            ))
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
    overrides = {
        "intended_cuts": _positions(args.intended_cuts) if args.intended_cuts is not None else None,
        "protected_positions": _positions(args.protected_positions) if args.protected_positions is not None else None,
        "mature_n_terminus": args.mature_n_terminus,
        "free_c_terminus": args.free_c_terminus,
    }
    overrides = {key: value for key, value in overrides.items() if value is not None}
    if overrides or args.secondary_structure is not None:
        if len(records) != 1:
            raise ValueError("site declarations and --secondary-structure require one input record or per-row manifest columns")
        record = records[0]
        records[0] = replace(
            record, annotations=replace(record.annotations, **overrides),
            secondary_structure=(
                record.secondary_structure if args.secondary_structure is None else args.secondary_structure
            ),
        )
    return records


def _constructs(record: Record, contexts: list[Context]) -> tuple[Record, ...]:
    sequence = normalize_sequence(record.sequence)
    base = replace(record, sequence=sequence)
    rows = [base]
    annotations = base.annotations
    validate_positions(annotations.intended_cuts, len(sequence) - 1, "intended cut")
    validate_positions(annotations.protected_positions, len(sequence), "protected residue")
    for context in contexts:
        insert_at = len(MET_EXCISION_PREFIX) if context.after_msg else 0
        if context.after_msg and not sequence.startswith(MET_EXCISION_PREFIX):
            raise ValueError(f"context {context.name!r} requires an input starting with MSG for MSG+INSERT")
        if (
            not context.after_msg and context.prefix.startswith(MET_EXCISION_PREFIX)
            and sequence.startswith(MET_EXCISION_PREFIX)
        ):
            raise ValueError(
                f"input already starts with MSG in context {context.name!r}. "
                "Use MSG+INSERT to add residues after that MSG or submit the complete construct."
            )
        construct = sequence[:insert_at] + context.prefix + sequence[insert_at:] + context.suffix
        added = len(context.prefix)
        changed_c_end = bool(context.suffix or (added and insert_at == len(sequence)))
        # An insertion breaks the bond immediately before it. Shift only existing bonds.
        cuts = tuple(
            position + added if position > insert_at else position
            for position in annotations.intended_cuts if not (added and position == insert_at)
        )
        protected = tuple(
            position + added if position > insert_at else position for position in annotations.protected_positions
        )
        rows.append(Record(
            record.identifier, construct, context=context.name,
            annotations=replace(
                annotations, intended_cuts=cuts, protected_positions=protected,
                mature_n_terminus=annotations.mature_n_terminus and (not added or context.after_msg),
                free_c_terminus=changed_c_end or annotations.free_c_terminus,
            ),
            secondary_structure=base.secondary_structure if construct == sequence else None,
        ))
    return tuple(rows)


def _describe(
    records: tuple[Record, ...],
    nfive_scores: dict[str, float] | None,
    *,
    groups: tuple[str, ...],
    ph: float,
    radius: float,
    sasa_points: int,
    organism: str,
    proteases_used: tuple[str, ...],
    include_details: bool,
) -> list[dict[str, object]]:
    sequence_groups = tuple(group for group in groups if group in SEQUENCE_GROUPS)
    surface_groups = tuple(group for group in groups if group in SURFACE_GROUPS)
    rows = []
    for record in records:
        try:
            details = {"sequence": record.sequence, "organism": organism} if include_details else None
            surface = dict.fromkeys(SURFACE_CHEMISTRY_FIELDS) if "surface-chemistry" in surface_groups else {}
            exposure = None
            if surface_groups and record.structure is not None:
                exposure = [] if include_details and SITE_FLAGS_GROUP in sequence_groups else None
                surface = describe_structure(
                    record.structure, chain=record.chain, expected_sequence=record.sequence,
                    radius=radius, sasa_points=sasa_points, groups=surface_groups,
                    details=details, residue_sasa_a2=exposure,
                    protected_positions=record.annotations.protected_positions,
                )
            elif details is not None and "surface-chemistry" in surface_groups:
                details["surface"] = {"reason": "structure_unavailable"}
            properties = _describe_normalized(
                record.sequence, groups=sequence_groups, ph=ph, organism=organism, nfive_scores=nfive_scores,
                annotations=record.annotations, proteases_used=proteases_used, secondary_structure=record.secondary_structure,
                residue_sasa_a2=exposure, details=details,
            )
            row = {"id": record.identifier, "context": record.context, **properties, **surface}
            if details is not None:
                row["details"] = details
            rows.append(row)
        except (OSError, ValueError) as exc:
            label = f"record {record.identifier!r}"
            if record.context != BASE_CONTEXT:
                label += f", context {record.context!r}"
            raise ValueError(f"{label}: {exc}") from exc
    return rows


def _write_outputs(rows: list[dict[str, object]], csv_path: Path, jsonl_path: Path | None) -> None:
    columns = list(dict.fromkeys(column for row in rows for column in row if column != "details"))
    with ExitStack() as stack:
        staged_csv = Path(stack.enter_context(TemporaryDirectory(dir=csv_path.parent))) / "properties.csv"
        with staged_csv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
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
        if len({context.name for context in args.context}) != len(args.context):
            raise ValueError("context names must be unique")
        if args.details and args.jsonl is None:
            raise ValueError("--details requires --jsonl")
        profile = get_organism(args.organism)
        if args.proteases_used is not None and SITE_FLAGS_GROUP not in args.groups:
            raise ValueError(f"--proteases-used requires the {SITE_FLAGS_GROUP} group")
        if N_DEGRON_STABILITY_GROUP in args.groups and not profile.nfive_supported:
            raise ValueError(f"{profile.name} has no registered N-FIVE model")
        if PI_TYPICALITY_GROUP in args.groups and profile.pi_reference is None:
            raise ValueError(f"{profile.name} has no pI reference cohort")
        proteases_used = normalize_proteases(args.proteases_used, organism=profile) if SITE_FLAGS_GROUP in args.groups else ()
        records = _records(args)
        structures = {record.structure for record in records if record.structure is not None}
        if "sap" in args.groups and not structures:
            raise ValueError("the sap group requires at least one structure path")
        csv_path = args.output.resolve()
        jsonl_path = args.jsonl.resolve() if args.jsonl is not None else None
        outputs = [path for path in (csv_path, jsonl_path) if path is not None]
        inputs = {path.resolve() for path in (*structures, args.fasta, args.manifest) if path is not None}
        if N_DEGRON_STABILITY_GROUP in args.groups or args.nfive_model is not None:
            inputs.add((args.nfive_model or DEFAULT_MODEL_PATH).expanduser().resolve())
        if args.nfive_lookup is not None:
            inputs.add(args.nfive_lookup.expanduser().resolve())
        if len(set(outputs)) != len(outputs) or inputs.intersection(outputs):
            raise ValueError("output paths must be distinct from each other and from inputs")
        if any(path.exists() and not path.is_file() for path in outputs):
            raise ValueError("outputs must be regular file paths")

        constructs = []
        for record in records:
            try:
                constructs.append(_constructs(record, args.context))
            except ValueError as exc:
                raise ValueError(f"record {record.identifier!r}: {exc}") from exc
        if N_DEGRON_STABILITY_GROUP in args.groups:
            all_scores = prepare_nfive_scores(
                (record.sequence for batch in constructs for record in batch),
                model_path=args.nfive_model, lookup_path=args.nfive_lookup,
            )
            # Send workers only the windows for their record rather than the full batch map.
            record_scores = (
                {window: all_scores[window] for record in batch for window in nfive_windows(record.sequence).values()}
                for batch in constructs
            )
        else:
            record_scores = repeat(None, len(constructs))
        describe_record = partial(
            _describe, groups=args.groups, ph=args.ph, radius=args.radius, sasa_points=args.sasa_points,
            organism=profile.name, proteases_used=proteases_used, include_details=args.details,
        )
        if args.workers == 1:
            rows = [
                row for batch, scores in zip(constructs, record_scores, strict=True)
                for row in describe_record(batch, scores)
            ]
        else:
            with ProcessPoolExecutor(max_workers=min(args.workers, len(records))) as pool:
                rows = [row for batch in pool.map(describe_record, constructs, record_scores) for row in batch]
        _write_outputs(rows, csv_path, jsonl_path)
    except (OSError, ValueError, csv.Error) as exc:
        parser.exit(2, f"{parser.prog}: error: {exc}\n")
    return 0
