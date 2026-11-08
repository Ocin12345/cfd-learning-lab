import unittest,sys
from pathlib import Path
R=Path('E:/CFD-Learning-Lab');W=Path(__file__).parent
sys.path[:0]=[str(W),str(R/'repo/cfd-learning-lab/scripts'),str(R/'local/finite-wing-continuation-v001'),str(R/'local/finite-wing-mesh-sensitivity-v001'),str(R/'local/showcase-improvements-2026-11-08-v001/sloshing-time-step-v001')]
from wing_repair_native_v003 import expected_settings
from wing_repair_native_v002 import expected_settings as previous
class Controls(unittest.TestCase):
 def test_only_declared_pressure_momentum_damping_changes(self):
  a=previous(.001,10.);b=expected_settings(.001,10.,.25)
  for k in ['explicit_pressure_under_relaxation','explicit_momentum_under_relaxation']:
   self.assertEqual(b['controls']['p_v_controls'][k],.25);a['controls']['p_v_controls'][k]=.25
  self.assertEqual(a,b)
 def test_comparator_settings_preserved(self):
  self.assertEqual(previous(.001,10.),expected_settings(.001,10.,.5))
 def test_reject_outside_declared_damping(self):
  for v in [0.,1.,float('nan')]:
   with self.assertRaises(ValueError):expected_settings(.001,10.,v)
if __name__=='__main__':unittest.main()
