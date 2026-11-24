"""Regression checks for shared exposure and helix-end descriptors."""

import json
from math import fsum
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

from Bio.SeqUtils.ProtParamData import bm

from protein_property_metrics.cli import main
from protein_property_metrics.helix import describe_helix
from protein_property_metrics.surface import ShrakeRupley, describe_structure


def _write_chain(path, residues, *, header="", ligand=""):
    atom_index = 0
    lines = [header]
    for position, (resname, sidechain) in enumerate(residues, 1):
        for atom_id in ("N", "CA", "C", "O", *sidechain):
            atom_index += 1
            lines.append(
                f"ATOM  {atom_index:5d} {atom_id:^4s} {resname} A{10 * position:4d}    "
                f"{20.0 * atom_index:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 10.00          {atom_id[0]:>2}\n"
            )
    path.write_text("".join(lines) + ligand + "END\n")


class SurfaceChemistryTests(TestCase):
    def test_shared_sasa_sulfur_selection_and_burial_clipping(self):
        areas = {("CYS", "SG"): 20.0, ("MET", "SD"): 80.0, ("MET", "CG"): 200.0}

        def assign_areas(protein, *, level):
            self.assertEqual(level, "A")
            for atom in protein.get_atoms():
                atom.sasa = areas.get((atom.parent.resname, atom.id), 0.0)

        with TemporaryDirectory() as temp:
            path = Path(temp) / "cys-met.pdb"
            _write_chain(path, [("CYS", ("CB", "SG")), ("MET", ("CB", "CG", "SD", "CE"))])
            exposure = []
            details = {}
            with patch.object(ShrakeRupley, "compute", side_effect=assign_areas) as compute:
                both = describe_structure(
                    path, groups=("sap", "surface-chemistry"), expected_sequence="CM",
                    details=details, residue_sasa_a2=exposure,
                )
                legacy = describe_structure(path, expected_sequence="CM")
                self.assertEqual(compute.call_count, 1)
                path.write_text(path.read_text().replace("  20.000", "  21.000"))
                describe_structure(path, expected_sequence="CM")
                self.assertEqual(compute.call_count, 2)

        self.assertEqual({key: both[key] for key in legacy}, legacy)
        cys_weight, met_weight = (bm[aa] - bm["G"] for aa in "CM")
        expected_sap = fsum((cys_weight * 20 / 81.91 / 6, met_weight * 280 / 150.32 / 8))
        self.assertAlmostEqual(both["sap_static_total"], expected_sap)
        self.assertEqual(both["cys_sulfur_sasa_a2"], 20.0)
        self.assertEqual(both["met_sulfur_sasa_a2"], 80.0)
        self.assertAlmostEqual(
            both["hydrophobic_burial_fraction"], cys_weight * (1 - 20 / 81.91) / (cys_weight + met_weight),
        )
        self.assertEqual(exposure, [20.0, 280.0])
        self.assertEqual([site["atom"] for site in details["surface"]["sulfur_sites"]], ["SG", "SD"])

    def test_covered_absence_and_incomplete_structure(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "gly.pdb"
            _write_chain(path, [("GLY", ())])
            metrics = describe_structure(path, groups=("surface-chemistry",), expected_sequence="G", sasa_points=128)
            self.assertEqual(metrics, {
                "cys_sulfur_sasa_a2": 0.0, "met_sulfur_sasa_a2": 0.0, "hydrophobic_burial_fraction": None,
            })
            with self.assertRaisesRegex(ValueError, "Sequence mismatch"):
                describe_structure(path, groups=("surface-chemistry",), expected_sequence="C", sasa_points=128)
            _write_chain(path, [("CYS", ("CB",))])
            with self.assertRaisesRegex(ValueError, "Missing heavy atoms SG"):
                describe_structure(path, groups=("surface-chemistry",), expected_sequence="C", sasa_points=128)
            _write_chain(path, [("GLY", ()), ("KCX", ())])
            path.write_text("".join(
                "HETATM" + line[6:] if line[17:20] == "KCX" else line
                for line in path.read_text().splitlines(keepends=True)
            ))
            with self.assertRaisesRegex(ValueError, "Unsupported or ambiguous protein residue.*KCX"):
                describe_structure(path, groups=("surface-chemistry",), sasa_points=128)

    def test_cli_protected_sulfur_without_site_flags(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            path = root / "cys.pdb"
            _write_chain(path, [("CYS", ("CB", "SG"))])
            jsonl = root / "metrics.jsonl"
            main([
                "--sequence", "C", "--structure", str(path), "--groups", "surface-chemistry",
                "--protected-positions", "1", "--details", "--jsonl", str(jsonl),
                "--output", str(root / "metrics.csv"),
            ])
            row = json.loads(jsonl.read_text())
        self.assertEqual(row["details"]["protected_positions"], [1])
        sulfur, = row["details"]["surface"]["sulfur_sites"]
        self.assertTrue(sulfur["protected"])
        self.assertGreater(sulfur["sasa_a2"], 0.0)
        self.assertEqual(sulfur["sasa_a2"], row["cys_sulfur_sasa_a2"])

    def test_annotated_disulfide_context(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "disulfide.pdb"
            _write_chain(
                path, [("CYS", ("CB", "SG")), ("CYS", ("CB", "SG"))],
                header="SSBOND   1 CYS A   10    CYS A   20\n",
            )
            details = {}
            describe_structure(path, groups=("surface-chemistry",), details=details, sasa_points=128)
        first, second = details["surface"]["sulfur_sites"]
        self.assertEqual(first["connections"][0]["kind"], "disulfide")
        self.assertEqual(first["connections"][0]["source"], "PDB SSBOND")
        self.assertEqual(first["connections"][0]["partner"]["residue_number"], 20)
        self.assertEqual(second["connections"][0]["partner"]["residue_number"], 10)

    def test_annotated_metal_context_keeps_isolated_exposure(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "zinc.pdb"
            residues = [("CYS", ("CB", "SG"))]
            _write_chain(path, residues)
            isolated = describe_structure(path, groups=("surface-chemistry",), sasa_points=128)
            link = f"LINK        {'SG':^4} CYS A{10:4d}{' ' * 16}{'ZN':^4}  ZN B{90:4d}\n"
            metal = f"HETATM{99:5d} {'ZN':^4}  ZN B{90:4d}    "
            metal += f"{121.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 10.00          ZN\n"
            _write_chain(path, residues, header=link, ligand=metal)
            details = {}
            annotated = describe_structure(path, groups=("surface-chemistry",), details=details, sasa_points=128)
        self.assertEqual(isolated, annotated)
        cysteine, = details["surface"]["sulfur_sites"]
        self.assertEqual(cysteine["connections"][0]["kind"], "metal_coordination")
        self.assertEqual(cysteine["connections"][0]["partner"]["atom"], "ZN")


class HelixTests(TestCase):
    def test_end_sign_disjoint_windows_and_histidine(self):
        details = {}
        self.assertEqual(describe_helix("DEHDKRRH", "HHHHHHHH", details=details), {
            "helix_end_charge_compatibility": 1.5,
        })
        helix, = details["helix"]["eligible_helices"]
        self.assertEqual((helix["start"], helix["end"]), (1, 8))
        self.assertEqual(helix["n_end_charges"], [-1, -1, 0, -1])
        self.assertEqual(helix["c_end_charges"], [1, 1, 1, 0])
        self.assertEqual(describe_helix("RRRRDDDDEEEE", "H" * 12)["helix_end_charge_compatibility"], -2.0)

    def test_missing_no_eligible_and_invalid_assignments(self):
        details = {}
        self.assertIsNone(describe_helix("GGGGGGGG", details=details)["helix_end_charge_compatibility"])
        self.assertEqual(details["helix"]["reason"], "assignment_unavailable")
        self.assertIsNone(describe_helix("GGGGGGGG", "HHHHHHHG", details=details)["helix_end_charge_compatibility"])
        self.assertEqual(details["helix"]["reason"], "no_eligible_helices")
        with self.assertRaisesRegex(ValueError, "length"):
            describe_helix("GGGGGGGG", "HHHHHHH")
        with self.assertRaisesRegex(ValueError, "unsupported"):
            describe_helix("GGGGGGGG", "HHHHHHHX")
