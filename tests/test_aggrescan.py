"""Regression fixtures for the pinned, web-compatible psico convention."""

from unittest import TestCase

from protein_property_metrics.aggrescan import describe_aggrescan


class AggrescanTests(TestCase):
    def test_length_three_uses_one_padded_window_and_keeps_negative_score(self):
        detail = {}
        metrics = describe_aggrescan("AAA", details=detail)
        self.assertAlmostEqual(metrics["aggrescan_na4vss"], -56.35)
        self.assertEqual(metrics["aggrescan_thsar"], 0)
        self.assertEqual(detail["aggrescan"]["hotspots"], [])

    def test_short_sequences_have_no_full_padded_window(self):
        for sequence in ("I", "II"):
            with self.subTest(sequence=sequence):
                detail = {}
                self.assertEqual(
                    describe_aggrescan(sequence, details=detail),
                    {"aggrescan_na4vss": None, "aggrescan_thsar": None},
                )
                self.assertIsNone(detail["aggrescan"]["hotspots"])
                self.assertIn("window", detail["aggrescan"]["reason"])
        with self.assertRaisesRegex(ValueError, "noncanonical"):
            describe_aggrescan("X")

    def test_positive_profile_without_five_consecutive_positions_has_no_hotspot(self):
        metrics = describe_aggrescan("IIII")
        self.assertAlmostEqual(metrics["aggrescan_na4vss"], 118.665)
        self.assertEqual(metrics["aggrescan_thsar"], 0)

    def test_terminal_hotspots_and_proline_break(self):
        # Five positive non-Pro positions on each side of a positive-profile Pro.
        # The pinned reference adds half the initial excess to the N-terminal
        # hotspot, with no corresponding bonus at the C terminus.
        detail = {}
        metrics = describe_aggrescan("IIIIIPIIIII", details=detail)
        self.assertAlmostEqual(metrics["aggrescan_na4vss"], 15.3446 / 11 * 100)
        self.assertAlmostEqual(metrics["aggrescan_thsar"], 14.7302 / 11)
        n_terminal, c_terminal = detail["aggrescan"]["hotspots"]
        self.assertEqual((n_terminal["start"], n_terminal["end"]), (1, 5))
        self.assertAlmostEqual(n_terminal["area"], 7.5456)
        self.assertEqual((c_terminal["start"], c_terminal["end"]), (7, 11))
        self.assertAlmostEqual(c_terminal["area"], 7.1846)

    def test_default_window_size_boundaries(self):
        for length, window in ((75, 5), (76, 7), (175, 7), (176, 9), (300, 9), (301, 11)):
            with self.subTest(length=length):
                detail = {}
                metrics = describe_aggrescan("A" * length, details=detail)
                # For uniform A, only the boundary windows differ from -0.036.
                # Their combined excess is -2.6375 / window and each is reused
                # window // 2 times by the pinned off-center convention.
                expected = -3.6 - 263.75 * (window // 2) / (length * window)
                self.assertEqual(detail["aggrescan"]["window_size"], window)
                self.assertAlmostEqual(metrics["aggrescan_na4vss"], expected)
                self.assertEqual(metrics["aggrescan_thsar"], 0)
