import unittest,numpy as np
from sloshing_cycles_v001 import cycle_evidence
class Cycles(unittest.TestCase):
 def test_exact_bracketed_triangle_and_partial_ends(self):
  t=np.array([0,.25,.5,.75,1,1.25,1.5,1.75,2]);y=np.array([-1,0,1,0,-1,0,1,0,-1],dtype=float)
  r=cycle_evidence(t,y,0);self.assertEqual(r['positive_crossing_times_s'],[.25,1.25]);self.assertEqual(len(r['complete_cycles']),1);self.assertEqual(r['complete_cycles'][0]['period_s'],1);self.assertEqual(r['complete_cycles'][0]['peak_to_peak'],2)
 def test_no_complete_cycle_and_no_extrapolation(self):
  r=cycle_evidence([0,1,2],[-1,1,2],0);self.assertEqual(r['positive_crossing_times_s'],[.5]);self.assertEqual(r['complete_cycles'],[])
 def test_bad_time_nonfinite_or_ambiguous_plateau(self):
  for t,y in [([0,0,1],[-1,0,1]),([0,1,2],[-1,float('nan'),1]),([0,1,2,3],[-1,0,0,1])]:
   with self.assertRaises(ValueError):cycle_evidence(t,y,0)
if __name__=='__main__':unittest.main()
