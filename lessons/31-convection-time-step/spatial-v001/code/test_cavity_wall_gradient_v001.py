import unittest,numpy as np
from cavity_wall_gradient_v001 import wall_gradient_heat
class Heat(unittest.TestCase):
 def test_linear_conduction_flux_and_refinement(self):
  for n in (6,12,48,96):
   x=(np.arange(n)+.5)/n;T=np.tile(300.4-.8*x,(n,1));q=wall_gradient_heat(T,.025478873239436622)
   self.assertAlmostEqual(q['left'],.025478873239436622*.8,places=12);self.assertAlmostEqual(q['right'],-q['left'],places=12);self.assertEqual(q['top'],0);self.assertEqual(q['bottom'],0)
 def test_invalid_geometry_or_material(self):
  for T,k in [(np.ones((3,4)),1),(np.full((4,4),np.nan),1),(np.ones((4,4)),0)]:
   with self.assertRaises(ValueError):wall_gradient_heat(T,k)
if __name__=='__main__':unittest.main()