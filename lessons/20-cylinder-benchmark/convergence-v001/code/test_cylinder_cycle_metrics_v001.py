import unittest
import numpy as np
from cylinder_cycle_metrics_v001 import cycle_metrics
class CycleMetrics(unittest.TestCase):
 def data(self):
  t=np.arange(0,4.0001,.005);return t,3.2+.03*np.cos(12*np.pi*t),np.sin(6*np.pi*t)
 def test_known_period_and_mean(self):
  t,cd,cl=self.data();m=cycle_metrics(t,cd,cl);self.assertTrue(m['last_six_cycles_stable_under_declared_criterion']);self.assertAlmostEqual(m['last_six_mean_period_s'],1/3,places=4);self.assertAlmostEqual(m['last_six_mean_Cd'],3.2,places=4)
 def test_growing_amplitude_not_stable(self):
  t,cd,cl=self.data();m=cycle_metrics(t,cd,cl*(1+.05*t));self.assertFalse(m['last_six_cycles_stable_under_declared_criterion'])
 def test_duplicate_time_rejected(self):
  t,cd,cl=self.data();t[2]=t[1]
  with self.assertRaises(ValueError):cycle_metrics(t,cd,cl)
 def test_insufficient_cycles_not_claimed(self):
  t,cd,cl=self.data();self.assertFalse(cycle_metrics(t[:100],cd[:100],cl[:100])['last_six_cycles_stable_under_declared_criterion'])
if __name__=='__main__':unittest.main()