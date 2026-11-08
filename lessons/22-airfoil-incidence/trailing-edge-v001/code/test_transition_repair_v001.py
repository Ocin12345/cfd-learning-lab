import unittest
import numpy as np
from local_refine_triangles_v002 import refine_mesh,edge
from triangle_quality_v002 import cell_quality
from quality_aware_refinement_v001 import repair_refinement
class TransitionRepair(unittest.TestCase):
 def data(self):
  xy=np.array([[0.,0.],[1,0],[1,.25],[0,.25]]);tri=np.array([[0,1,2],[0,2,3]]);b={edge(a,c):4 for a,c in [(0,1),(1,2),(2,3),(3,0)]};return xy,tri,b,np.array([True,False])
 def test_actual_transition_pattern_is_repaired(self):
  xy,t,b,m=self.data();q,tt,bb,p=refine_mesh(xy,t,b,m);self.assertLess(cell_quality(q,tt).min(),.22)
  q,tt,bb,p=repair_refinement(xy,t,b,m,np.array([.5,.125]),1.);self.assertGreaterEqual(cell_quality(q,tt).min(),.22);self.assertTrue(np.array_equal(q[:len(xy)],xy));self.assertAlmostEqual(p['area_after'],.25,places=12)
 def test_excessive_extent_rejected(self):
  xy,t,b,m=self.data()
  with self.assertRaises(ValueError):repair_refinement(xy,t,b,m,np.array([.5,.125]),.1)
if __name__=='__main__':unittest.main()