import unittest
from pathlib import Path
import audit_jet_energy_qualification_v002 as a
class Clock(unittest.TestCase):
 def test_actual_current_clock_ignores_saved_history(self):
  root=Path('E:/CFD-Learning-Lab/native/runs/showcase-heated-jet-energy-v001/fluent-qualification-v001')
  for s in range(21):a.clock(root/f'heat-{s:04d}.dat.h5',s)
 def test_wrong_expected_current_time_still_rejected(self):
  p=Path('E:/CFD-Learning-Lab/native/runs/showcase-heated-jet-energy-v001/fluent-qualification-v001/heat-0001.dat.h5')
  with self.assertRaises(AssertionError):a.clock(p,0)
if __name__=='__main__':unittest.main()
