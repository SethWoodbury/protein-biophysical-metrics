"""Check the sequence windows sent to the N-FIVE model."""

import csv
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

import numpy as np

from protein_property_metrics.cli import main
from protein_property_metrics.sequence import describe


class NFiveTests(TestCase):
    @patch("protein_property_metrics.nfive._load_model")
    def test_msg_windows_in_base_and_tagged_context(self, load_model):
        # Distinct outputs verify the row and window associated with each prediction.
        predict = load_model.return_value.predict
        predict.side_effect = [np.array([2.0]), np.array([2.25, 2.75])]
        with TemporaryDirectory() as temp:
            output = Path(temp) / "properties.csv"
            main([
                "--sequence", "acdefgh", "--groups", "N-degron-stability-index",
                "--context", "tagged:msg:", "--output", str(output),
            ])
            with output.open(newline="") as handle:
                base, tagged = csv.DictReader(handle)

        self.assertEqual(base["n_degron_sequence"], "ACDEF")
        self.assertEqual(base["n_degron_stability_index"], "2.0")
        self.assertEqual(base["n_degron_sequence_met_excised"], "")
        self.assertEqual(base["n_degron_stability_index_met_excised"], "")
        self.assertEqual(tagged["n_degron_sequence"], "MSGAC")
        self.assertEqual(tagged["n_degron_stability_index"], "2.25")
        self.assertEqual(tagged["n_degron_sequence_met_excised"], "SGACD")
        self.assertEqual(tagged["n_degron_stability_index_met_excised"], "2.75")
        encoded = predict.call_args.args[0]
        self.assertEqual(encoded.shape, (2, 100))
        np.testing.assert_array_equal(np.flatnonzero(encoded[0]), [12, 35, 47, 60, 84])
        np.testing.assert_array_equal(np.flatnonzero(encoded[1]), [15, 27, 40, 64, 83])

    @patch("protein_property_metrics.nfive._load_model")
    def test_short_windows_and_literal_msg(self, load_model):
        group = ("n-degron-stability-index",)
        short = describe("MSGK", groups=group)
        self.assertTrue(all(value is None for value in short.values()))
        load_model.assert_not_called()

        load_model.return_value.predict.side_effect = [np.array([2.0]), np.array([2.0, 2.5])]
        five = describe("MSGKA", groups=group)
        self.assertEqual(five["n_degron_sequence"], "MSGKA")
        self.assertIsNone(five["n_degron_sequence_met_excised"])
        self.assertIsNone(five["n_degron_stability_index_met_excised"])
        six = describe("MSGKAL", groups=group)
        self.assertEqual(six["n_degron_sequence_met_excised"], "SGKAL")
        self.assertEqual(six["n_degron_stability_index_met_excised"], 2.5)
