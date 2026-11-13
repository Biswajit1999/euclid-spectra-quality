import unittest
from research.redshift_threshold_20261113 import residuals,wilson,threshold_sweep,leave_one_source_out
class TestRedshift(unittest.TestCase):
 def test_residual(self): self.assertAlmostEqual(residuals([1.1],[1.])[0],0.05)
 def test_wilson(self): self.assertGreater(wilson(0,10)[1],0)
 def test_threshold(self): self.assertEqual(threshold_sweep([0.01,0.1],[0.01])[0][1],1)
 def test_influence(self): self.assertEqual(leave_one_source_out([0,0.2,0],["a","b","c"],0.1)["b"],0)
 def test_duplicates(self):
  with self.assertRaises(ValueError): leave_one_source_out([0,1],["a","a"],0.1)
