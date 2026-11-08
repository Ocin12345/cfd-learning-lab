import unittest
from pathlib import Path
import audit_jet_energy_qualification_v001 as old
class Clock(unittest.TestCase):
 def test_actual_current_clock_ignores_saved_history(self):
  p=Path('E:/CFD-Learning-Lab/native/runs/showcase-heated-jet-energy-v001/fluent-qualification-v001/heat-0001.dat.h5')
  old.clock(p,1)
if __name__=='__main__':unittest.main()
