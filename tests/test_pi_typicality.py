"""Check pI typicality, host selection, and public output."""

import csv
import json
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

from protein_property_metrics.cli import main
from protein_property_metrics.sequence import describe


class PITypicalityTests(TestCase):
    def test_both_local_density_modes_rank_above_a_gap(self):
        from protein_property_metrics.pi_typicality import PIReference

        reference = PIReference((4.9, 5.0, 5.1, 7.0, 9.9, 10.0, 10.1), 0.25)
        self.assertEqual(reference.percentile(5.0), 100.0)
        self.assertEqual(reference.percentile(10.0), 100.0)
        self.assertAlmostEqual(reference.percentile(7.0), 100 / 7)
        self.assertEqual(reference.percentile(2.0), 0.0)

    def test_pi_calculated_once_for_both_groups(self):
        from protein_property_metrics.sequence import IsoelectricPoint

        sequence = "ACDEKLRW"
        expected = describe(sequence, groups=("charge",))["isoelectric_point"]
        calls = []
        original = IsoelectricPoint.pi

        def pi(self, *args, **kwargs):
            calls.append(self)
            return original(self, *args, **kwargs)

        with patch.object(IsoelectricPoint, "pi", pi):
            both = describe(sequence, groups=("charge", "pi-typicality"))
        # Biopython recurses internally while locating the pI. Compare the
        # call count against calculating charge with the same search bounds.
        typicality_calls = len(calls)
        calls.clear()
        with patch.object(IsoelectricPoint, "pi", pi):
            describe(sequence, groups=("charge",))
        self.assertEqual(typicality_calls, len(calls))
        self.assertEqual(both["isoelectric_point"], expected)
        self.assertGreaterEqual(both["pI_typicality_percentile"], 0)
        self.assertLessEqual(both["pI_typicality_percentile"], 100)
        self.assertNotIn("pI_typicality_percentile", describe(sequence))

    def test_other_profiles_select_their_own_capabilities(self):
        from protein_property_metrics.organisms import get_organism

        organism = replace(get_organism("BL21_DE3_e_coli"), name="test_host", pi_reference=None, nfive_supported=False)
        with self.assertRaisesRegex(ValueError, "test_host.*pI reference"):
            describe("ACDEKLRW", groups=("pi-typicality",), organism=organism)
        with self.assertRaisesRegex(ValueError, "test_host.*N-FIVE"):
            describe("ACDEKLRW", groups=("n-degron-stability-index",), organism=organism)

    def test_cli_emits_reference_and_separates_protease_counts(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            output, detail_output = root / "metrics.csv", root / "metrics.jsonl"
            main([
                "--organism", "BL21_DE3_e_coli", "--sequence", "KKENLYFQG", "--proteases-used", "tev,ompt",
                "--groups", "charge,pi-typicality,site-flags", "--output", str(output),
                "--jsonl", str(detail_output), "--details",
            ])
            with output.open(newline="") as handle:
                row, = csv.DictReader(handle)
            metadata = json.loads(detail_output.read_text())["details"]["pi_typicality"]
        self.assertEqual(row["protease_cut_site_count"], "1")
        self.assertEqual(row["context_dependent_protease_site_count"], "1")
        self.assertTrue(row["pI_typicality_percentile"])
        self.assertEqual(metadata["organism"], "BL21_DE3_e_coli")
        self.assertEqual(metadata["proteome_id"], "UP000002032")
        self.assertGreater(metadata["reference_count"], 100)
