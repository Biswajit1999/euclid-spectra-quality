import unittest
from group_dependence_20261110 import source_level_pull_summary, design_effect


class GroupDependenceTests(unittest.TestCase):
    def test_duplicate_source_not_double_counted(self):
        a = source_level_pull_summary(["A","A","B"],[1.1,1.1,1.0],[1,1,1],[.1,.1,.1])
        self.assertEqual(a["rows"],3)
        self.assertEqual(a["unique_sources"],2)
        self.assertAlmostEqual(a["source_mean_pull"],0.5)

    def test_zero_pull(self):
        self.assertEqual(source_level_pull_summary(["A"],[1],[1],[.1])["source_rms_pull"],0)

    def test_bad_sigma(self):
        with self.assertRaises(ValueError):
            source_level_pull_summary(["A"],[1],[1],[0])

    def test_independent_limit(self):
        self.assertEqual(design_effect(20,5,0)["effective_n"],20)

    def test_perfect_cluster_limit(self):
        self.assertEqual(design_effect(20,5,1)["effective_n"],4)

    def test_invalid_cluster(self):
        with self.assertRaises(ValueError):
            design_effect(19,5,0)


if __name__ == "__main__":
    unittest.main()
