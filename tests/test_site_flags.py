"""Checks for construct-relative site flags and compact occurrence counts."""

from dataclasses import replace
from unittest import TestCase

from protein_property_metrics.organisms import DEFAULT_ORGANISM, OrganismProfile, get_organism
from protein_property_metrics.site_flags import describe_sites


class SiteFlagTests(TestCase):
    def test_chemical_counts_protection_exposure_and_declared_n_terminus(self):
        sequence = "QNGNSDGDPQENG"
        sasa = [10.0] * len(sequence)
        sasa[1] = 0.0
        sasa[3] = None
        details = {}
        result = describe_sites(
            sequence,
            protected_positions=(2, 8),
            mature_n_terminus=True,
            residue_sasa_a2=sasa,
            details=details,
        )
        self.assertEqual((result["ng_count"], result["dg_count"], result["ns_count"]), (2, 1, 1))
        sites = details["site_flags"]["sites"]
        by_start = {site["start"]: site for site in sites}
        self.assertEqual(set(by_start), {1, 2, 4, 6, 8, 12})
        self.assertEqual([(by_start[i]["motif_id"], by_start[i]["end"]) for i in (2, 4, 6, 12)],
                         [("NG", 3), ("NS", 5), ("DG", 7), ("NG", 13)])
        self.assertTrue(by_start[2]["protected"])
        self.assertEqual(by_start[2]["sasa_a2"][0], 0.0)
        self.assertIsNone(by_start[4]["sasa_a2"][0])
        self.assertEqual(by_start[8]["motif_id"], "DP")
        self.assertTrue(by_start[8]["protected"])
        self.assertEqual(by_start[1]["motif_id"], "N-terminal-Q")
        details = {}
        describe_sites(sequence, details=details)
        self.assertNotIn(1, {s["start"] for s in details["site_flags"]["sites"]})

    def test_n_degron_signal_classes_and_msg(self):
        for residue, motif_id in (("F", "ClpS-primary"), ("R", "N-end-secondary")):
            details = {}
            describe_sites(residue + "SGA", mature_n_terminus=True, details=details)
            site, = details["site_flags"]["sites"]
            self.assertEqual(site["motif_id"], motif_id)
            self.assertEqual((site["start"], site["end"]), (1, 1))
        for sequence in ("MSGF", "SGFR", "QG"):
            details = {}
            describe_sites(sequence, mature_n_terminus=sequence != "QG", details=details)
            self.assertFalse(details["site_flags"]["sites"])

    def test_polyproline_overlaps_and_adjacent_clusters(self):
        details = {}
        result = describe_sites("PPPGAPPPAPPPPG", details=details)
        self.assertEqual(result["polyproline_cluster_count"], 3)
        flags = details["site_flags"]
        self.assertEqual(flags["longest_proline_run"], 4)
        self.assertEqual([(site["start"], site["end"]) for site in flags["polyproline_clusters"]],
                         [(1, 4), (6, 8), (10, 14)])
        self.assertEqual([(s["motif_id"], s["start"], s["end"]) for s in flags["sites"]], [
            ("PPP", 1, 3), ("PPG", 2, 4), ("PPP", 6, 8),
            ("PPP", 10, 12), ("PPP", 11, 13), ("PPG", 12, 14),
        ])
        self.assertEqual(describe_sites("PPGPPG")["polyproline_cluster_count"], 2)

    def test_terminal_degron_tags_longest_minor_and_blocked_or_internal_sites(self):
        for variant, suffix in (
            ("ssrA-native", "AANDENYALAA"),
            ("ssrA-LVA", "AANDENYALVA"),
            ("ssrA-AAV", "AANDENYAAAV"),
            ("ssrA-ASV", "AANDENYAASV"),
        ):
            with self.subTest(suffix=suffix):
                details = {}
                result = describe_sites("SG" + suffix, details=details)
                self.assertEqual((result["degron_tag_count"], result["degron_motif_count"]), (1, 0))
                site, = details["site_flags"]["sites"]
                self.assertEqual((site["motif_id"], site["matched_sequence"], site["start"], site["end"]),
                                 (variant, suffix, 3, 13))
                self.assertIsNone(site["cleavage_after"])
                for sequence, free_end in (("SG" + suffix + "GG", True), ("SG" + suffix, False)):
                    result = describe_sites(sequence, free_c_terminus=free_end)
                    self.assertEqual((result["degron_tag_count"], result["degron_motif_count"]), (0, 0))
        for motif in ("LAA", "VAA", "AA", "RRKKAI"):
            details = {}
            result = describe_sites("SG" + motif, details=details)
            self.assertEqual((result["degron_tag_count"], result["degron_motif_count"]), (0, 1))
            site, = details["site_flags"]["sites"]
            self.assertEqual(site["matched_sequence"], motif)
            self.assertEqual(describe_sites("SG" + motif + "G")["degron_motif_count"], 0)

    def test_workflow_offsets_intention_protection_and_metadata(self):
        sequence = "AENLYFQGLEVLFQGPLVPRGSDDDDKAIEGRV"
        details = {}
        result = describe_sites(
            sequence,
            proteases_used="TEV,HRV 3C,Thrombin,Enterokinase,Factor Xa",
            intended_cuts=(14, 27),
            protected_positions=(7, 20, 28),
            details=details,
        )
        flags = details["site_flags"]
        by_cut = {
            site["cleavage_after"]: (site["protease"], site["matched_sequence"], site["intended"], site["protected"])
            for site in flags["sites"]
        }
        self.assertEqual(flags["proteases_used"], ["tev", "hrv-3c", "thrombin", "enterokinase", "factor-xa"])
        self.assertEqual(by_cut, {
            7: ("tev", "ENLYFQG", False, True),
            14: ("hrv-3c", "LEVLFQGP", True, False),
            20: ("thrombin", "LVPRGS", False, True),
            27: ("enterokinase", "DDDDKA", True, True),
            32: ("factor-xa", "IEGRV", False, False),
        })
        self.assertEqual(result["protease_cut_site_count"], 3)
        self.assertIsNone(describe_sites(sequence, proteases_used=())["protease_cut_site_count"])
        self.assertEqual(describe_sites("AAA", proteases_used=("tev",))["protease_cut_site_count"], 0)

    def test_overlapping_protease_matches_and_unique_bonds(self):
        sequence = "DDDDKDDDDKA"
        details = {}
        result = describe_sites(sequence, proteases_used=("enterokinase",), details=details)
        self.assertEqual(result["protease_cut_site_count"], 2)
        self.assertEqual([site["cleavage_after"] for site in details["site_flags"]["sites"]], [5, 10])
        # Two rules share the same offset and match here to exercise deduplication.
        profile = get_organism(DEFAULT_ORGANISM)
        rules = dict(profile.proteases, tev=profile.proteases["enterokinase"])
        result = describe_sites(
            sequence, organism=replace(profile, proteases=rules),
            proteases_used=("tev", "enterokinase"), intended_cuts=(5,), details=details,
        )
        self.assertEqual(result["protease_cut_site_count"], 1)
        hits = details["site_flags"]["sites"]
        self.assertEqual({site["protease"] for site in hits if site["cleavage_after"] == 10}, {"tev", "enterokinase"})

    def test_bl21_default_tev_and_optional_ompt_candidates(self):
        sequence = "AENLYFQGAKKRRA"
        details = {}
        result = describe_sites(sequence, details=details)
        self.assertEqual(result["protease_cut_site_count"], 1)
        self.assertNotIn("context_dependent_protease_site_count", result)
        flags = details["site_flags"]
        self.assertEqual(flags["organism"], DEFAULT_ORGANISM)
        self.assertEqual(flags["proteases_used"], ["tev"])
        tev, = flags["sites"]
        self.assertEqual(tev["bin"], "workflow_protease")

        result = describe_sites(sequence, proteases_used="TEV,OmpT", intended_cuts=(11,), details=details)
        self.assertEqual(result["protease_cut_site_count"], 1)
        self.assertEqual(result["context_dependent_protease_site_count"], 2)
        candidates = [site for site in details["site_flags"]["sites"] if site["protease"] == "ompt"]
        self.assertEqual([site["matched_sequence"] for site in candidates], ["KK", "KR", "RR"])
        self.assertEqual([site["cleavage_after"] for site in candidates], [10, 11, 12])
        self.assertEqual([site["intended"] for site in candidates], [False, True, False])
        self.assertTrue(all(site["bin"] == "context_dependent_protease" for site in candidates))

        result = describe_sites(sequence, proteases_used="ompt")
        self.assertIsNone(result["protease_cut_site_count"])
        self.assertEqual(result["context_dependent_protease_site_count"], 3)
        self.assertEqual(describe_sites("AA", proteases_used="ompt")["context_dependent_protease_site_count"], 0)
        self.assertIsNone(describe_sites(sequence, proteases_used="NONE")["protease_cut_site_count"])

    def test_organism_specific_rules_can_be_replaced(self):
        profile = OrganismProfile(name="other_host")
        details = {}
        result = describe_sites("FENLYFQGAANDENYALAA", organism=profile, mature_n_terminus=True, details=details)
        self.assertEqual(result["degron_tag_count"], 0)
        self.assertIsNone(result["protease_cut_site_count"])
        self.assertEqual(details["site_flags"]["organism"], "other_host")
        self.assertEqual(details["site_flags"]["sites"], [])
        self.assertEqual(get_organism("bl21_de3_e_coli"), get_organism(DEFAULT_ORGANISM))
        with self.assertRaisesRegex(ValueError, "unknown organism"):
            get_organism("unsupported")

    def test_required_following_residues_and_specificity(self):
        for residue in "SGAMCH":
            self.assertEqual(describe_sites("ENLYFQ" + residue, proteases_used=("tev",))["protease_cut_site_count"], 1)
        self.assertEqual(describe_sites("IDGRA", proteases_used=("factor-xa",))["protease_cut_site_count"], 1)
        for protease, prefix, rejected in (
            ("tev", "ENLYFQ", "P"),
            ("hrv-3c", "LEVLFQ", "A"),
            ("thrombin", "LVPR", "AS"),
            ("enterokinase", "DDDDK", "P"),
            ("factor-xa", "IDGR", "P"),
            ("factor-xa", "IEGR", "R"),
        ):
            for sequence in (prefix, prefix + rejected):
                with self.subTest(protease=protease, sequence=sequence):
                    self.assertEqual(describe_sites(sequence, proteases_used=(protease,))["protease_cut_site_count"], 0)

    def test_bad_metadata_is_rejected(self):
        for kwargs in (
            {"intended_cuts": (3,)},
            {"protected_positions": (0,)},
            {"protected_positions": (True,)},
            {"proteases_used": "unsupported"},
            {"proteases_used": "none,tev"},
            {"residue_sasa_a2": [1.0, 2.0]},
            {"residue_sasa_a2": [1.0, float("nan"), 2.0]},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                describe_sites("ADC", **kwargs)
