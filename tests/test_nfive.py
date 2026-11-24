"""Check the sequence windows sent to the N-FIVE model."""

import csv
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

import numpy as np

from protein_property_metrics.cli import main
from protein_property_metrics.nfive import build_msg_lookup, describe_n_terminus, prepare_nfive_scores
from protein_property_metrics.sequence import describe


class NFiveTests(TestCase):
    def setUp(self):
        scratch = TemporaryDirectory()
        self.addCleanup(scratch.cleanup)
        self.directory = Path(scratch.name)
        self.model = self.directory / "model.pkl"
        # Distinct fake model content isolates the process-local score cache.
        self.model.write_bytes(self.id().encode())

    @patch("protein_property_metrics.nfive._load_model")
    def test_msg_windows_in_base_and_tagged_context(self, load_model):
        # Distinct outputs verify the row and window associated with each prediction.
        predict = load_model.return_value.predict
        predict.return_value = np.array([2.0, 2.25, 2.75])
        output = self.directory / "properties.csv"
        main([
            "--sequence", "acdefgh", "--groups", "N-degron-stability-index",
            "--context", "tagged:msg:", "--output", str(output),
            "--nfive-model", str(self.model),
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
        predict.assert_called_once()
        encoded = predict.call_args.args[0]
        self.assertEqual(encoded.shape, (3, 100))
        np.testing.assert_array_equal(np.flatnonzero(encoded[1]), [12, 35, 47, 60, 84])
        np.testing.assert_array_equal(np.flatnonzero(encoded[2]), [15, 27, 40, 64, 83])

    @patch("protein_property_metrics.nfive._load_model")
    def test_short_windows_and_literal_msg(self, load_model):
        group = ("n-degron-stability-index",)
        short = describe("MSGK", groups=group)
        self.assertTrue(all(value is None for value in short.values()))
        load_model.assert_not_called()

        predict = load_model.return_value.predict
        predict.side_effect = [np.array([2.0]), np.array([2.5])]
        five = describe("MSGKA", groups=group, nfive_model=self.model)
        self.assertEqual(five["n_degron_sequence"], "MSGKA")
        self.assertIsNone(five["n_degron_sequence_met_excised"])
        self.assertIsNone(five["n_degron_stability_index_met_excised"])
        six = describe("MSGKAL", groups=group, nfive_model=self.model)
        self.assertEqual(six["n_degron_stability_index"], 2.0)
        self.assertEqual(six["n_degron_sequence_met_excised"], "SGKAL")
        self.assertEqual(six["n_degron_stability_index_met_excised"], 2.5)
        self.assertEqual(predict.call_args.args[0].shape, (1, 100))

    @patch("protein_property_metrics.nfive._load_model")
    def test_unique_batch_and_model_change_invalidate_cache(self, load_model):
        predict = load_model.return_value.predict
        predict.side_effect = [np.array([1.5, 2.5, 3.5]), np.array([4.0])]
        scores = prepare_nfive_scores(["ACDEFG", "MSGACDE", "MSGACDM", "SGACDL"], model_path=self.model)
        self.assertEqual(scores, {"ACDEF": 1.5, "MSGAC": 2.5, "SGACD": 3.5})
        self.assertEqual(predict.call_args.args[0].shape, (3, 100))
        row = describe_n_terminus("MSGACDE", prepared_scores=scores)
        self.assertEqual(row["n_degron_stability_index_met_excised"], 3.5)
        self.assertEqual(prepare_nfive_scores(["ACDEFM"], model_path=self.model), {"ACDEF": 1.5})
        predict.assert_called_once()

        self.model.write_bytes(b"changed model")
        self.assertEqual(prepare_nfive_scores(["ACDEFG"], model_path=self.model), {"ACDEF": 4.0})
        self.assertEqual(predict.call_count, 2)

    @patch("protein_property_metrics.nfive._load_model")
    def test_msg_lookup_and_direct_score_parity(self, load_model):
        lookup = self.directory / "msg.npz"
        predict = load_model.return_value.predict
        predict.side_effect = lambda encoded: encoded @ np.arange(100, dtype=np.float32)
        build_msg_lookup(lookup, model_path=self.model)
        with np.load(lookup, allow_pickle=False) as data:
            self.assertEqual(data["alphabet"].item(), "ARNDCEQGHILKMFPSTWYV")
            self.assertEqual(data["msg"].shape, (20, 20))
            self.assertEqual(data["sg"].shape, (20, 20, 20))
            self.assertEqual(data["msg"].dtype, np.float32)
            self.assertEqual(data["sg"].dtype, np.float32)

        tagged = ["MSGAAAE", "MSGVVVF", "MSGACDE", "SGACDK"]
        predict.reset_mock()
        scores = prepare_nfive_scores(tagged, model_path=self.model, lookup_path=lookup)
        predict.assert_not_called()
        self.assertEqual(scores, prepare_nfive_scores(tagged, model_path=self.model))
        predict.assert_called_once()
        mixed = prepare_nfive_scores(["ACDEFG", "MSGACDE"], model_path=self.model, lookup_path=lookup)
        self.assertEqual(mixed["SGACD"], scores["SGACD"])
        self.assertEqual(predict.call_args.args[0].shape, (1, 100))

        self.model.write_bytes(b"other model")
        with self.assertRaisesRegex(ValueError, "model SHA256"):
            prepare_nfive_scores(tagged, model_path=self.model, lookup_path=lookup)
