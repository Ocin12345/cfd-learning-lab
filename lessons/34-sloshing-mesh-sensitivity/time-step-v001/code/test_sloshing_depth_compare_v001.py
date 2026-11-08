import unittest,numpy as np
from sloshing_depth_compare_v001 import depth_comparison
class Depth(unittest.TestCase):
 def test_physical_strip_average_and_statistics(self):
  old=np.full((3,4),.05);new=old+np.array([[0,0,0,0],[1,3,5,7],[2,4,6,8]])*1e-6
  r=depth_comparison([0,.01,.02],old,new,.05)
  self.assertAlmostEqual(r['maximum_absolute_native_column_difference_m'],8e-6)
  self.assertAlmostEqual(r['right_strip_difference_m'][-1],7e-6)
  self.assertEqual(r['right_strip_definition'],'Mean of last two fine native columns; same physical 2.5mm strip')
 def test_bad_geometry_or_nonfinite_rejected(self):
  for a,b in [(np.zeros((3,3)),np.zeros((3,3))),(np.zeros((3,4)),np.full((3,4),np.nan)),(np.zeros((3,4)),np.zeros((2,4)))]:
   with self.assertRaises(ValueError):depth_comparison([0,.01,.02],a,b,.05)
if __name__=='__main__':unittest.main()
