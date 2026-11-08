import unittest,numpy as np
from wing_cell_identity_v001 import centroid_bijection
class CellIdentity(unittest.TestCase):
 def test_permutation_and_rounding_boundary(self):
  a=np.array([[5e-9,0,0],[1,2,3],[2,3,4.]])
  b=a[[2,0,1]].copy();b[1,0]-=1e-15
  i,d=centroid_bijection(a,b);self.assertEqual(i.tolist(),[2,0,1]);self.assertLess(float(d.max()),2e-10)
 def test_duplicate_refused(self):
  with self.assertRaises(ValueError):centroid_bijection(np.array([[0.,0,0],[0.,0,0]]),np.array([[0.,0,0],[0.,0,0]]))
 def test_missing_cell_refused(self):
  with self.assertRaises(ValueError):centroid_bijection(np.array([[0.,0,0],[1.,0,0]]),np.array([[0.,0,0],[2.,0,0]]))
if __name__=='__main__':unittest.main()
