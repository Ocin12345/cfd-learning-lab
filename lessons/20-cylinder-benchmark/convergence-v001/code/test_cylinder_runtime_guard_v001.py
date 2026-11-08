from pathlib import Path
import unittest
from cylinder_runtime_guard_v001 import check_native_step
N=Path('E:/CFD-Learning-Lab/native/runs/showcase-cylinder-convergence-v001')
class NativeRuntimeGuard(unittest.TestCase):
 def paths(self,name):
  p=N/name;return p/'vortex-final.dat.h5',p/'force-reports/mass-000810.txt'
 def test_real_pass(self):
  d,m=self.paths('fluent-coupled-pilot-v001');r=check_native_step(d,m,4.05);self.assertTrue(r['passed'])
 def test_real_residual_miss_rejected(self):
  d,m=self.paths('fluent-pressure-baseline-v001')
  with self.assertRaises(ValueError):check_native_step(d,m,4.05)
 def test_wrong_native_time_rejected(self):
  d,m=self.paths('fluent-coupled-pilot-v001')
  with self.assertRaises(ValueError):check_native_step(d,m,4.055)
if __name__=='__main__':unittest.main()