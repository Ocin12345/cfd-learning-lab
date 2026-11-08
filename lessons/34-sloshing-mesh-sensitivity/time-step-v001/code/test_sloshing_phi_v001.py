import unittest
from sloshing_phi_v001 import phi_wall_flux
class Flux(unittest.TestCase):
 def test_zero_wall_flux_and_internal_count(self):
  text='dimensions [0 3 -1 0 0 0 0]; internalField nonuniform List<scalar> 2 (1 -1); boundaryField { tankWalls { type calculated; value uniform 0; } }'
  self.assertEqual(phi_wall_flux(text,4,2),0.)
 def test_reject_wrong_units_count_and_nonfinite(self):
  text='dimensions [0 3 -1 0 0 0 0]; internalField nonuniform List<scalar> 2 (1 -1); boundaryField { tankWalls { value uniform 0; } }'
  for t in [text.replace('[0 3 -1','[1 0 -1'),text.replace('2 (1 -1)','3 (1 -1)'),text.replace('1 -1','nan -1')]:
   with self.assertRaises(AssertionError):phi_wall_flux(t,4,2)
if __name__=='__main__':unittest.main()
