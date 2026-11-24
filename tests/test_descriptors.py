"""Small numerical checks for the descriptor conventions."""

from math import dist, isclose, pi
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from Bio.PDB import MMCIFIO, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from protein_property_metrics.sequence import describe
from protein_property_metrics.surface import describe_structure


class SequenceTests(TestCase):
    def test_free_terminal_charge_and_pi(self):
        metrics = describe("GG", groups=("charge",), ph=7.5)
        self.assertAlmostEqual(metrics["isoelectric_point"], 5.525, places=4)
        self.assertAlmostEqual(metrics["net_charge"], 0.5 - 1 / (1 + 10 ** (3.55 - 7.5)))
        self.assertAlmostEqual(metrics["absolute_net_charge"], abs(metrics["net_charge"]))
        self.assertGreater(describe("GG", groups=("charge",), ph=4)["net_charge"], 0)
        self.assertLess(describe("GG", groups=("charge",), ph=7)["net_charge"], 0)

    def test_terminal_residue_and_wide_pi_bracket(self):
        self.assertAlmostEqual(describe("GA", groups=("charge",))["isoelectric_point"], 5.525, places=4)
        self.assertAlmostEqual(describe("AG", groups=("charge",))["isoelectric_point"], 5.570, places=4)
        self.assertLess(describe("DDD", groups=("charge",))["isoelectric_point"], 4.05)
        self.assertGreater(describe("RRR", groups=("charge",))["isoelectric_point"], 12)

    def test_optics_and_aliphatic_index(self):
        optics = describe("CCCYW", groups=("optics",))
        self.assertEqual(optics["extinction_280_reduced_m_inv_cm_inv"], 6990)
        self.assertEqual(optics["extinction_280_max_paired_m_inv_cm_inv"], 7115)
        self.assertAlmostEqual(describe("AVIL", groups=("composition",))["aliphatic_index"], 292.5)


class SurfaceTests(TestCase):
    atoms = (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"))
    radii_with_probe = (2.95, 3.10, 3.10, 2.92)

    def write_glycine(self, path, coordinates):
        lines = []
        for i, ((name, element), (x, y, z)) in enumerate(zip(self.atoms, coordinates), 1):
            lines.append(
                f"ATOM  {i:5d} {name:^4s} GLY A   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{10.0:6.2f}          {element:>2s}\n"
            )
        path.write_text("".join(lines) + "END\n")

    def test_disjoint_spheres(self):
        with TemporaryDirectory() as temp:
            path = Path(temp) / "separate.pdb"
            self.write_glycine(path, [(20 * i, 0, 0) for i in range(4)])
            metrics = describe_structure(path, expected_sequence="G")
        expected = sum(4 * pi * r * r for r in self.radii_with_probe)
        self.assertTrue(isclose(metrics["surface_total_sasa_a2"], expected, rel_tol=1e-6))
        self.assertEqual(metrics["surface_sidechain_sasa_a2"], 0)
        self.assertEqual(metrics["surface_ca_radius_of_gyration_a"], 0)

    def test_overlapping_spheres_and_mmcif(self):
        coords = [(0, 0, 0), (2.3, 1.4, 0.8), (20, 0, 0), (40, 0, 0)]
        with TemporaryDirectory() as temp:
            pdb = Path(temp) / "overlap.pdb"
            cif = Path(temp) / "overlap.cif"
            self.write_glycine(pdb, coords)
            source = PDBParser(QUIET=True).get_structure("test", pdb)
            exporter = MMCIFIO()
            exporter.set_structure(source)
            exporter.save(str(cif))
            fields = MMCIF2Dict(cif)
            for column in ("B_iso_or_equiv", "pdbx_PDB_model_num"):
                fields.pop(f"_atom_site.{column}")
            exporter.set_dict(fields)
            exporter.save(str(cif))
            metrics = describe_structure(pdb, expected_sequence="G", sasa_points=4096)
            self.assertEqual(metrics, describe_structure(cif, expected_sequence="G", sasa_points=4096))

        expected = sum(4 * pi * r * r for r in self.radii_with_probe)
        distance = dist(*coords[:2])
        r1, r2 = self.radii_with_probe[:2]
        for r, other in ((r1, r2), (r2, r1)):
            plane = (distance * distance + r * r - other * other) / (2 * distance)
            expected -= 2 * pi * r * (r - plane)
        self.assertTrue(isclose(metrics["surface_total_sasa_a2"], expected, rel_tol=0.002))
