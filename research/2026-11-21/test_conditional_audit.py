import unittest
import numpy as np
from pathlib import Path
import importlib.util

spec = importlib.util.spec_from_file_location("audit",Path(__file__).with_name("conditional_audit.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class TestAudit(unittest.TestCase):
    def test_normalization(self):
        np.testing.assert_allclose(m.normalized_redshift_error([1.2],[1.]), [.1])

    def test_stratified_exact(self):
        e = [0,.3,0,.3,0,.3,0,.3]
        f = [0,1,0,1,0,1,0,1]
        g = [0,0,0,0,1,1,1,1]
        val, ci = m.stratified_outlier_difference(e,f,g,n_boot=100,seed=3)
        self.assertEqual(val,1.)
        np.testing.assert_allclose(ci,[1,1])

    def test_missing_cell_rejected(self):
        with self.assertRaises(ValueError):
            m.stratified_outlier_difference([0,.2],[False,True],[0,1])

    def test_coverage(self):
        self.assertEqual(m.interval_coverage([1,2],[1,2.5],[.1,.1]),.5)

    def test_nonpositive_sigma_rejected(self):
        with self.assertRaises(ValueError):
            m.interval_coverage([1],[1],[0])


if __name__ == "__main__":
    unittest.main()
