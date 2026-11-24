"""Static solvent exposure and hydrophobic patch descriptors."""

from __future__ import annotations

from collections.abc import Iterable
from functools import lru_cache
from io import StringIO
from itertools import islice
from math import fsum, isfinite
from operator import index
from pathlib import Path
from typing import NamedTuple

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import MMCIFParser, PDBParser, is_aa
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.Residue import Residue
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.kdtrees import KDTree
from Bio.SeqUtils.ProtParamData import bm

from .site_flags import validate_positions

DEFAULT_RADIUS = 5.0
DEFAULT_SASA_POINTS = 960
SURFACE_GROUPS = ("sap", "surface-chemistry")
SURFACE_CHEMISTRY_FIELDS = ("cys_sulfur_sasa_a2", "met_sulfur_sasa_a2", "hydrophobic_burial_fraction")
_PROBE_RADIUS = 1.4
_STRUCTURE_SUFFIXES = {".cif", ".mmcif", ".pdb", ".ent"}

_BACKBONE = {"N", "CA", "C", "O"}
_TERMINAL = {"OXT"}
_SULFUR_ATOMS = {"C": "SG", "M": "SD"}
_PDB_LINK_METALS = frozenset("LI NA K RB CS MG CA SR BA AL V CR MN FE CO NI CU ZN MO CD W HG PB".split())
_SIDECHAIN_ATOMS = {
    "A": "CB",
    "C": "CB SG",
    "D": "CB CG OD1 OD2",
    "E": "CB CG CD OE1 OE2",
    "F": "CB CG CD1 CD2 CE1 CE2 CZ",
    "G": "",
    "H": "CB CG ND1 CD2 CE1 NE2",
    "I": "CB CG1 CG2 CD1",
    "K": "CB CG CD CE NZ",
    "L": "CB CG CD1 CD2",
    "M": "CB CG SD CE",
    "N": "CB CG OD1 ND2",
    "P": "CB CG CD",
    "Q": "CB CG CD OE1 NE2",
    "R": "CB CG CD NE CZ NH1 NH2",
    "S": "CB OG",
    "T": "CB OG1 CG2",
    "V": "CB CG1 CG2",
    "W": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "Y": "CB CG CD1 CD2 CE1 CE2 CZ OH",
}

# Side-chain area (Å²) of the central residue in Ala–X–Ala, Samanta et al. (2002),
# Table II, doi:10.1093/protein/15.8.659. Glycine has no heavy side-chain atoms.
_REFERENCE_SIDECHAIN_SASA = {
    "A": 55.28,
    "C": 81.91,
    "D": 97.66,
    "E": 132.42,
    "F": 164.14,
    "H": 141.17,
    "I": 130.58,
    "K": 147.83,
    "L": 141.47,
    "M": 150.32,
    "N": 109.87,
    "P": 98.23,
    "Q": 129.52,
    "R": 190.04,
    "S": 68.97,
    "T": 88.45,
    "V": 103.00,
    "W": 209.57,
    "Y": 180.01,
}


def _is_polymer_record(residue: Residue) -> bool:
    # Modified amino acids can be encoded as HETATM records.
    return residue.id[0] == " " or (
        residue.resname not in protein_letters_3to1 and is_aa(residue)
    )


class _MMCIFStructureBuilder(StructureBuilder):
    def init_model(self, model_id: int, serial_num: int | None = None) -> None:
        # Biopython can initialize the same implicit model on each atom row
        # when the optional _atom_site.pdbx_PDB_model_num column is absent.
        if serial_num is None and self.model is not None and self.model.id == model_id:
            return
        super().init_model(model_id, serial_num)


class _CoordinateMMCIFParser(MMCIFParser):
    def _build_structure(self, structure_id):
        # Generated coordinates may omit temperature factors. SASA uses
        # coordinates and elements, so supply the parser's required placeholder.
        data = self._mmcif_dict
        if "_atom_site.B_iso_or_equiv" not in data:
            data["_atom_site.B_iso_or_equiv"] = ["0"] * len(data["_atom_site.id"])
        super()._build_structure(structure_id)


class _AtomAddress(NamedTuple):
    chain: str
    residue_number: int
    insertion_code: str
    residue_name: str
    atom: str


class _Connection(NamedTuple):
    kind: str
    source: str
    first: _AtomAddress
    second: _AtomAddress


def _atom_address(residue: Residue, atom: str) -> _AtomAddress:
    return _AtomAddress(residue.parent.id, residue.id[1], residue.id[2].strip(), residue.resname, atom)


def _pdb_connections(content: str, model: Model) -> tuple[_Connection, ...]:
    # Fixed-width connectivity fields follow the wwPDB coordinate format 3.3.
    # https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html
    connections = []
    metal_atoms = None
    for line in StringIO(content):
        try:
            if line.startswith("SSBOND"):
                first = _AtomAddress(line[15:16], int(line[17:21]), line[21:22].strip(), line[11:14].strip(), "SG")
                second = _AtomAddress(line[29:30], int(line[31:35]), line[35:36].strip(), line[25:28].strip(), "SG")
                connections.append(_Connection("disulfide", "PDB SSBOND", first, second))
            elif line.startswith("LINK  "):
                first = _AtomAddress(line[21:22], int(line[22:26]), line[26:27].strip(), line[17:20].strip(), line[12:16].strip())
                second = _AtomAddress(line[51:52], int(line[52:56]), line[56:57].strip(), line[47:50].strip(), line[42:46].strip())
                if first.atom == second.atom == "SG" and first.residue_name == second.residue_name == "CYS":
                    connections.append(_Connection("disulfide", "PDB LINK", first, second))
                elif first.atom in _SULFUR_ATOMS.values() or second.atom in _SULFUR_ATOMS.values():
                    if metal_atoms is None:
                        metal_atoms = {
                            _atom_address(atom.parent, atom.id)
                            for atom in model.get_atoms()
                            if atom.element in _PDB_LINK_METALS
                        }
                    if first in metal_atoms or second in metal_atoms:
                        connections.append(_Connection("metal_coordination", "PDB LINK", first, second))
        except ValueError:
            # Incomplete annotations do not identify a coordinate site.
            continue
    return tuple(connections)


def _mmcif_connections(data: dict) -> tuple[_Connection, ...]:
    connections = []
    for row, connection_type in enumerate(data.get("_struct_conn.conn_type_id", ())):
        kind = {"disulf": "disulfide", "metalc": "metal_coordination"}.get(connection_type)
        if kind is None:
            continue
        partners = []
        for partner in (1, 2):
            prefix = f"_struct_conn.ptnr{partner}_"
            try:
                chain_id = data[prefix + "auth_asym_id"][row]
                number = int(data[prefix + "auth_seq_id"][row])
                residue = data[prefix + "label_comp_id"][row]
                atom = data[prefix + "label_atom_id"][row]
                insertions = data.get(f"_struct_conn.pdbx_ptnr{partner}_PDB_ins_code")
                insertion = insertions[row] if insertions is not None else ""
            except (KeyError, IndexError, ValueError):
                break
            partners.append(_AtomAddress(chain_id, number, "" if insertion in {".", "?"} else insertion, residue, atom))
        if len(partners) == 2:
            connections.append(_Connection(kind, "mmCIF _struct_conn", *partners))
    return tuple(connections)


def _load_protein_chain(path: Path, chain_id: str | None, content: str) -> tuple[Chain, str, tuple[_Connection, ...]]:
    if path.suffix.lower() in {".cif", ".mmcif"}:
        parser = _CoordinateMMCIFParser(structure_builder=_MMCIFStructureBuilder(), QUIET=True)
    else:
        parser = PDBParser(PERMISSIVE=False, QUIET=True)

    try:
        structure = parser.get_structure(path.stem, StringIO(content))
    except (PDBConstructionException, IndexError, KeyError, ValueError) as error:
        raise ValueError(f"Cannot parse structure {path}: {error}") from error
    models = list(structure.get_models())
    if len(models) != 1:
        raise ValueError(f"Expected one coordinate model in {path}, found {len(models)}")

    protein_chains = {
        chain.id: chain
        for chain in models[0]
        if any(_is_polymer_record(residue) and is_aa(residue) for residue in chain)
    }
    if chain_id is None:
        if len(protein_chains) != 1:
            choices = ", ".join(repr(ch) for ch in protein_chains) or "none"
            raise ValueError(f"Expected one protein chain in {path}. Protein chains: {choices}")
        chain = next(iter(protein_chains.values()))
    else:
        if chain_id not in protein_chains:
            raise ValueError(f"Protein chain {chain_id!r} is absent from {path}")
        chain = protein_chains[chain_id]

    if isinstance(parser, _CoordinateMMCIFParser):
        connections = _mmcif_connections(parser._mmcif_dict)
    else:
        connections = _pdb_connections(content, models[0])
    sequence = []
    for residue in list(chain):
        if not _is_polymer_record(residue):
            chain.detach_child(residue.id)
            continue
        aa = protein_letters_3to1.get(residue.resname)
        label = (
            f"{path}, chain {chain.id!r}, residue {residue.resname} "
            f"{residue.id[1]}{residue.id[2].strip()}"
        )
        if aa is None or residue.is_disordered() == 2:
            raise ValueError(f"Unsupported or ambiguous protein residue in {label}")

        for atom in list(residue):
            if atom.element in {"H", "D"}:
                residue.detach_child(atom.id)
            elif atom.is_disordered():
                raise ValueError(f"Alternate coordinates for {atom.id} in {label}")

        atoms = {atom.id: atom for atom in residue}
        required = _BACKBONE.union(_SIDECHAIN_ATOMS[aa].split())
        missing = required - atoms.keys()
        if missing:
            raise ValueError(f"Missing heavy atoms {', '.join(sorted(missing))} in {label}")
        unexpected = atoms.keys() - required - _TERMINAL
        if unexpected:
            raise ValueError(f"Unexpected heavy atoms {', '.join(sorted(unexpected))} in {label}")
        for atom_id, atom in atoms.items():
            if atom.element != atom_id[0]:
                raise ValueError(f"Unexpected element {atom.element!r} for {atom_id} in {label}")
        sequence.append(aa)

    if not sequence:
        raise ValueError(f"Protein chain {chain.id!r} has no standard coordinate residues in {path}")
    # The exposure cache should retain only the selected and validated chain.
    models[0].detach_child(chain.id)
    return chain, "".join(sequence), connections


def describe_structure(
    path: str | Path,
    *,
    chain: str | None = None,
    expected_sequence: str | None = None,
    radius: float = DEFAULT_RADIUS,
    sasa_points: int = DEFAULT_SASA_POINTS,
    groups: Iterable[str] = ("sap",),
    details: dict | None = None,
    residue_sasa_a2: list[float] | None = None,
    protected_positions: Iterable[int] = (),
) -> dict[str, int | float | None]:
    """Describe one isolated protein chain from a single PDB or mmCIF model.

    ``chain`` is the author chain ID. When omitted, the file must contain one
    protein chain. ``expected_sequence`` checks the entire coordinate-derived
    sequence exactly after uppercasing and stripping outer whitespace. Each
    scored residue must have its canonical heavy atoms with blank alternate
    location identifiers. Waters, ligands, other chains, and explicit
    hydrogens are excluded from exposure.

    This static approximation to Chennamsetty et al. (2009, doi:10.1073/
    pnas.0904191106) uses Shrake–Rupley atomic SASA in Å² with a 1.4 Å probe and
    ``sasa_points`` surface samples per atom. For each heavy-atom center, it
    sums over side-chain atoms at distances up to ``radius`` Å, including the
    center itself if it is a side-chain atom. Each contribution is the atom's
    SASA divided by its residue type's Ala–X–Ala reference side-chain SASA
    (Samanta et al., 2002), multiplied by Black–Mould hydrophobicity with the
    glycine value subtracted (doi:10.1016/0003-2697(91)90045-U). Glycine can
    center a neighborhood and contributes no side-chain area. The exposure
    ratio is uncapped. Scores describe the supplied conformation.

    Residue scores average their heavy-atom center scores. ``sap_static_total``
    sums the positive residue scores. ``sap_static_per_residue`` divides that
    sum by ``sap_static_covered_residues``, the number of scored residues.
    ``sap_static_max`` is the maximum signed residue score. All SAP values are
    dimensionless. ``surface_total_sasa_a2`` and ``surface_sidechain_sasa_a2``
    report the total and side-chain heavy-atom SASA in Å².
    ``surface_ca_radius_of_gyration_a`` uses the unweighted Cα centroid.

    Add ``surface-chemistry`` to ``groups`` for sulfur-atom SASA and a bounded
    hydrophobic burial fraction. Both groups share atomic SASA. An optional
    ``residue_sasa_a2`` list receives summed heavy-atom areas in sequence order.
    A ``details`` dict receives structure scope and reported sulfur sites with
    any usable PDB or mmCIF connection annotation, independent of isolation
    for SASA. Supply 1-based ``protected_positions`` to mark those sites.
    Absence of a connection leaves the sulfur's chemical state unknown.
    """
    selected = set(groups)
    if unknown := selected.difference(SURFACE_GROUPS):
        raise ValueError(f"unknown structure descriptor groups ({', '.join(sorted(unknown))})")
    if not isfinite(radius) or radius <= 0:
        raise ValueError("radius must be finite and positive")
    if isinstance(sasa_points, bool):
        raise ValueError("sasa_points must be a positive integer")
    try:
        sasa_points = index(sasa_points)
    except TypeError as error:
        raise ValueError("sasa_points must be a positive integer") from error
    if sasa_points <= 0:
        raise ValueError("sasa_points must be a positive integer")

    path = Path(path)
    if path.suffix.lower() not in _STRUCTURE_SUFFIXES:
        raise ValueError(f"Expected a PDB or mmCIF file, got {path}")
    expected = expected_sequence.strip().upper() if expected_sequence is not None else None
    protein, sequence, coords, connections = _measure_exposure(path, path.read_text(), chain, expected, sasa_points)
    protected = validate_positions(protected_positions, len(sequence), "protected residue")
    if residue_sasa_a2 is not None:
        residue_sasa_a2.extend(fsum(atom.sasa for atom in residue) for residue in protein)
    surface_details = None
    if details is not None:
        surface_details = {
            "atom_scope": "isolated_canonical_heavy_atom_chain",
            "chain": protein.id,
            "sasa_points": sasa_points,
            "probe_radius_a": _PROBE_RADIUS,
        }
        details["surface"] = surface_details

    metrics = {}
    if "sap" in selected:
        metrics.update(_sap_metrics(protein, sequence, coords, radius))
    if "surface-chemistry" in selected:
        metrics.update(_surface_chemistry(protein, sequence, connections, surface_details, protected))
    return metrics


@lru_cache(maxsize=4)
def _measure_exposure(
    path: Path, content: str, chain: str | None, expected: str | None, sasa_points: int,
) -> tuple[Chain, str, np.ndarray, tuple[_Connection, ...]]:
    """Cache only private exposure data, keyed by file text and selection."""
    protein, sequence, connections = _load_protein_chain(path, chain, content)
    if expected is not None and sequence != expected:
        first_difference = next(
            (i for i, (actual, wanted) in enumerate(zip(sequence, expected), 1) if actual != wanted),
            min(len(sequence), len(expected)) + 1,
        )
        raise ValueError(
            f"Sequence mismatch for chain {protein.id!r} in {path}. "
            f"First difference at position {first_difference}. "
            f"Expected {len(expected)} residues, found {len(sequence)}."
        )
    coords = np.asarray([atom.coord for atom in protein.get_atoms()], dtype=np.float64)
    if not np.isfinite(coords).all():
        raise ValueError(f"Non-finite coordinates in chain {protein.id!r} of {path}")
    ShrakeRupley(probe_radius=_PROBE_RADIUS, n_points=sasa_points).compute(protein, level="A")
    return protein, sequence, coords, connections


def _sidechain_sasa(atom) -> float:
    return 0.0 if atom.id in _BACKBONE or atom.id in _TERMINAL else atom.sasa


def _sap_metrics(protein: Chain, sequence: str, coords: np.ndarray, radius: float) -> dict[str, int | float]:
    tree = KDTree(coords)
    weighted_areas = []
    sidechain_areas = []
    for residue, aa in zip(protein, sequence, strict=True):
        scale = (bm[aa] - bm["G"]) / _REFERENCE_SIDECHAIN_SASA[aa] if aa != "G" else 0.0
        for atom in residue:
            sidechain_area = _sidechain_sasa(atom)
            sidechain_areas.append(sidechain_area)
            weighted_areas.append(scale * sidechain_area)

    atom_scores = (
        fsum(weighted_areas[neighbor.index] for neighbor in tree.search(center, radius))
        for center in coords
    )
    residue_scores = [fsum(islice(atom_scores, len(residue))) / len(residue) for residue in protein]
    total = fsum(max(score, 0.0) for score in residue_scores)
    covered = len(residue_scores)
    ca_coords = np.array([residue["CA"].coord for residue in protein], dtype=np.float64)
    centered = ca_coords - ca_coords.mean(axis=0)
    total_sasa = fsum(atom.sasa for atom in protein.get_atoms())
    return {
        "sap_static_total": total,
        "sap_static_per_residue": total / covered,
        "sap_static_max": max(residue_scores),
        "sap_static_covered_residues": covered,
        "surface_total_sasa_a2": total_sasa,
        "surface_sasa_per_residue_a2": total_sasa / covered,
        "surface_sidechain_sasa_a2": fsum(sidechain_areas),
        "surface_ca_radius_of_gyration_a": float(np.sqrt(np.square(centered).sum(axis=1).mean())),
    }


def _surface_chemistry(
    protein: Chain, sequence: str, connections: tuple[_Connection, ...], details: dict | None,
    protected: set[int],
) -> dict[str, float | None]:
    sulfur_areas: dict[str, list[float]] = {aa: [] for aa in _SULFUR_ATOMS}
    weights = []
    buried_weights = []
    site_details = []
    connections_by_atom: dict[_AtomAddress, list[dict]] = {}
    if details is not None:
        for kind, source, first, second in connections:
            for site, partner in ((first, second), (second, first)):
                connections_by_atom.setdefault(site, []).append({
                    "kind": kind, "source": source, "partner": partner._asdict(),
                })

    for position, (residue, aa) in enumerate(zip(protein, sequence, strict=True), 1):
        atom_id = _SULFUR_ATOMS.get(aa)
        if atom_id is not None:
            area = float(residue[atom_id].sasa)
            sulfur_areas[aa].append(area)
            if details is not None:
                address = _atom_address(residue, atom_id)
                site_details.append({
                    "position": position,
                    "residue": aa,
                    "atom": atom_id,
                    "residue_number": address.residue_number,
                    "insertion_code": address.insertion_code,
                    "sasa_a2": area,
                    "protected": position in protected,
                    "connections": connections_by_atom.get(address, []),
                })
        weight = bm[aa] - bm["G"]
        if weight > 0.0:
            sidechain_area = fsum(_sidechain_sasa(atom) for atom in residue)
            ratio = min(1.0, max(0.0, sidechain_area / _REFERENCE_SIDECHAIN_SASA[aa]))
            weights.append(weight)
            buried_weights.append(weight * (1.0 - ratio))

    if details is not None:
        details["sulfur_sites"] = site_details
    denominator = fsum(weights)
    return {
        "cys_sulfur_sasa_a2": fsum(sulfur_areas["C"]),
        "met_sulfur_sasa_a2": fsum(sulfur_areas["M"]),
        "hydrophobic_burial_fraction": fsum(buried_weights) / denominator if denominator else None,
    }
