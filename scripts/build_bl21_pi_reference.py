"""Build the BL21(DE3) pI reference from downloaded, same-release UniProt data."""

import argparse
import hashlib
import json
from pathlib import Path

import Bio
from Bio import SeqIO

from protein_property_metrics.pi_typicality import PIReference
from protein_property_metrics.sequence import describe, normalize_sequence

_CYTOPLASM = "SL-0086"


def build(native_fasta: Path, reviewed_json: Path, release: str, half_window: float) -> dict:
    native = {str(record.seq): record.id for record in SeqIO.parse(native_fasta, "fasta")}
    reviewed = json.loads(reviewed_json.read_text(encoding="utf-8"))["results"]
    records = []
    skipped_noncanonical = 0
    for entry in sorted(reviewed, key=lambda entry: entry["primaryAccession"]):
        if entry.get("entryType") != "UniProtKB reviewed (Swiss-Prot)":
            raise ValueError(f"Expected reviewed UniProtKB entries. Check {entry['primaryAccession']}.")
        locations = {
            location["location"].get("id")
            for comment in entry.get("comments", ()) if comment["commentType"] == "SUBCELLULAR LOCATION"
            for location in comment.get("subcellularLocations", ())
            if "location" in location
        }
        if locations != {_CYTOPLASM}:
            continue
        sequence = entry["sequence"]["value"]
        native_id = native.pop(sequence, None)
        if native_id is None:
            continue
        try:
            sequence = normalize_sequence(sequence)
        except ValueError:
            skipped_noncanonical += 1
            continue
        pi = describe(sequence, groups=("charge",))["isoelectric_point"]
        records.append([native_id, entry["primaryAccession"], pi])

    reference = PIReference(tuple(record[2] for record in records), half_window)
    return {
        "proteome_id": "UP000002032",
        "annotation_proteome_id": "UP000000625",
        "uniprot_release": release,
        "half_window_ph": reference.half_window,
        "pi_model": f"Biopython {Bio.__version__} Bjellqvist, free termini, package charge group",
        "native_fasta_sha256": hashlib.sha256(native_fasta.read_bytes()).hexdigest(),
        "reviewed_json_sha256": hashlib.sha256(reviewed_json.read_bytes()).hexdigest(),
        "noncanonical_matches_skipped": skipped_noncanonical,
        "records": sorted(records),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-fasta", type=Path, required=True, help="BL21(DE3) UniParc proteome FASTA")
    parser.add_argument("--reviewed-json", type=Path, required=True, help="Reviewed K-12 UniProtKB cytoplasm JSON")
    parser.add_argument("--uniprot-release", required=True, help="Release shared by the two downloads")
    parser.add_argument("--half-window", type=float, default=0.25, help="Density half-window in pH units. Default is 0.25.")
    parser.add_argument("--output", type=Path, required=True, help="New JSON resource path")
    args = parser.parse_args()
    reference = build(args.native_fasta, args.reviewed_json, args.uniprot_release, args.half_window)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(reference, handle, separators=(",", ":"))
        handle.write("\n")
    print(f"Wrote {len(reference['records'])} reference pIs to {args.output}")
