"""Check the additive public API, context handling, and CSV/JSONL outputs."""

import csv
import json
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from protein_property_metrics.cli import main
from protein_property_metrics.sequence import DEFAULT_SEQUENCE_GROUPS, describe


class AddedMetricsTests(TestCase):
    def test_swi_weights_and_composition(self):
        residues = "ACDEFGHIKLMNPQRSTVWY"
        weights = [
            .836, .521, .908, .988, .585, .800, .895, .678, .927, .655,
            .630, .860, .824, .789, .771, .744, .810, .736, .637, .611,
        ]
        for residue, weight in zip(residues, weights, strict=True):
            self.assertAlmostEqual(describe(residue, groups=("swi",))["swi_solubility_index"], weight)
        self.assertAlmostEqual(describe(" ac de ", groups=("swi",))["swi_solubility_index"], .81325)
        default = describe("MCM")
        self.assertEqual(DEFAULT_SEQUENCE_GROUPS, ("composition", "charge", "instability"))
        self.assertEqual(default["methionine_count"], 2)
        self.assertEqual(default["cysteine_count"], 1)
        self.assertNotIn("swi_solubility_index", default)

    def test_jsonl_details_and_explicit_constructs(self):
        with TemporaryDirectory() as temp:
            output, jsonl = Path(temp) / "properties.csv", Path(temp) / "properties.jsonl"
            main([
                "--sequence", "MSGNGQ", "--groups", "site-flags,surface-chemistry",
                "--protected-positions", "4", "--context", "inserted:MSG+AA:LAA",
                "--output", str(output), "--jsonl", str(jsonl), "--details",
            ])
            with output.open(newline="") as handle:
                base_csv, context_csv = csv.DictReader(handle)
            base, context = map(json.loads, jsonl.read_text().splitlines())

        self.assertEqual(base_csv["ng_count"], "1")
        self.assertEqual(context_csv["degron_motif_count"], "1")
        self.assertEqual(base_csv["cys_sulfur_sasa_a2"], "")
        self.assertNotIn("details", base_csv)
        self.assertEqual(base["details"]["sequence"], "MSGNGQ")
        self.assertEqual(context["details"]["sequence"], "MSGAANGQLAA")
        self.assertIsNone(context["cys_sulfur_sasa_a2"])
        self.assertEqual(context["protease_cut_site_count"], 0)
        site_details = context["details"]["site_flags"]
        self.assertEqual(context["details"]["protected_positions"], [6])
        self.assertEqual(site_details["protected_positions"], [6])
        ng = next(site for site in site_details["sites"] if site["matched_sequence"] == "NG")
        self.assertEqual((ng["start"], ng["end"], ng["protected"]), (6, 7, True))

    def test_duplicate_literal_msg_requires_explicit_context(self):
        with TemporaryDirectory() as temp, redirect_stderr(StringIO()) as stderr:
            output = Path(temp) / "properties.csv"
            with self.assertRaises(SystemExit) as error:
                main([
                    "--sequence", "MSGACD", "--context", "vector:MSG:",
                    "--output", str(output),
                ])
            self.assertEqual(error.exception.code, 2)
            self.assertFalse(output.exists())
        self.assertIn("already starts with MSG", stderr.getvalue())

    def test_new_suffix_defines_its_own_c_terminus(self):
        with TemporaryDirectory() as temp:
            output = Path(temp) / "properties.csv"
            main([
                "--sequence", "MPLAA", "--groups", "site-flags", "--blocked-c-terminus",
                "--context", "tagged::VAA", "--output", str(output),
            ])
            with output.open(newline="") as handle:
                base, context = csv.DictReader(handle)
        self.assertEqual(base["degron_motif_count"], "0")
        self.assertEqual(context["degron_motif_count"], "1")
