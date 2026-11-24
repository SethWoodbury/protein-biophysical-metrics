"""Regression checks for shared exposure and helix-end descriptors."""

import json
from math import fsum
from pathlib import Path
from subprocess import CompletedProcess, TimeoutExpired
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

from Bio.SeqUtils.ProtParamData import bm

from protein_property_metrics import surface
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

    @staticmethod
    def dssp_summary(codes: str) -> str:
        header = (
            "data_dssp\nloop_\n_dssp_struct_summary.label_asym_id\n_dssp_struct_summary.label_seq_id\n"
            "_dssp_struct_summary.label_comp_id\n_dssp_struct_summary.secondary_structure\n"
        )
        return header + "".join(
            f"A {position} {'ASP' if position <= 4 else 'LYS'} {code}\n"
            for position, code in enumerate(codes, 1)
        )

    @staticmethod
    def write_helix_input(path: Path) -> None:
        _write_chain(path, [("ASP", ("CB", "CG", "OD1", "OD2"))] * 4 + [("LYS", ("CB", "CG", "CD", "CE", "NZ"))] * 4)

    def test_cli_dssp_uses_coordinates_without_sasa_and_preserves_caller_assignments(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            path = root / "input.pdb"
            self.write_helix_input(path)
            jsonl = root / "metrics.jsonl"
            args = [
                "--sequence", "DDDDKKKK", "--structure", str(path), "--groups", "helix", "--dssp",
                "--context", "tagged:MSG:", "--jsonl", str(jsonl), "--details", "--output", str(root / "metrics.csv"),
            ]
            result = CompletedProcess([], 0, stdout=self.dssp_summary("H" * 8), stderr="")
            with (
                patch("protein_property_metrics.dssp.subprocess.run", return_value=result) as dssp,
                patch.object(ShrakeRupley, "compute") as compute,
            ):
                main(args)
                base, tagged = map(json.loads, jsonl.read_text().splitlines())
                self.assertEqual(base["helix_end_charge_compatibility"], 2.0)
                self.assertEqual(base["details"]["helix"]["assignment_source"], "dssp")
                self.assertIsNone(tagged["helix_end_charge_compatibility"])
                main([*args, "--secondary-structure", "........"])
                supplied, _ = map(json.loads, jsonl.read_text().splitlines())
                self.assertIsNone(supplied["helix_end_charge_compatibility"])
                dssp.assert_called_once()
                compute.assert_not_called()

    def test_dssp_requires_complete_assignments_and_keeps_structure_validation(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "input.pdb"
            self.write_helix_input(path)
            result = CompletedProcess([], 0, stdout=self.dssp_summary("H" * 7), stderr="")
            with patch("protein_property_metrics.dssp.subprocess.run", return_value=result) as dssp:
                with self.assertRaisesRegex(ValueError, "DSSP.*cover"):
                    describe_structure(path, groups=(), dssp_executable="mkdssp", expected_sequence="DDDDKKKK")
                with self.assertRaisesRegex(ValueError, "Sequence mismatch"):
                    describe_structure(path, groups=(), dssp_executable="mkdssp", expected_sequence="GGGGGGGG")
                dssp.assert_called_once()

    def test_combined_dssp_shares_structure_parsing_and_enforces_timeout(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "input.pdb"
            self.write_helix_input(path)
            result = CompletedProcess([], 0, stdout=self.dssp_summary("H" * 8), stderr="")
            with (
                patch("protein_property_metrics.dssp.subprocess.run", return_value=result) as dssp,
                patch.object(surface, "_load_protein_chain", wraps=surface._load_protein_chain) as load,
            ):
                metrics = describe_structure(path, dssp_executable="mkdssp", dssp_timeout=2.0, sasa_points=64)
                self.assertEqual(metrics["helix_end_charge_compatibility"], 2.0)
                self.assertIn("sap_static_total", metrics)
                self.assertEqual(dssp.call_args.kwargs["timeout"], 2.0)
                load.assert_called_once()
            error = TimeoutExpired(["mkdssp"], timeout=2.0)
            with patch("protein_property_metrics.dssp.subprocess.run", side_effect=error):
                with self.assertRaisesRegex(ValueError, "DSSP exceeded the 2-second timeout"):
                    describe_structure(path, groups=(), dssp_executable="mkdssp", dssp_timeout=2.0)
