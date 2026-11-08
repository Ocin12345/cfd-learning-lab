from pathlib import Path
import unittest
from analyse_vortex_openfoam_v002 import residual_steps
from audit_cylinder_openfoam_v001 import verify_residual_window
class NativeOpenFOAMCoverage(unittest.TestCase):
 def test_actual_old_window_and_missing_step(self):
  p=Path('E:/CFD-Learning-Lab/native/runs/2026-11-01-vortex-shedding-v001/openfoam-extended/solver-continuation-to-4s.log');text=p.read_text();rows=[r for r in residual_steps(text,p.name) if r['time_s']>3.95+1e-8]
  self.assertEqual(len(verify_residual_window(rows,3.95,10)),10)
  with self.assertRaises(ValueError):verify_residual_window(rows[:-1],3.95,10)
 def test_missed_terminal_residual_rejected(self):
  p=Path('E:/CFD-Learning-Lab/native/runs/2026-11-01-vortex-shedding-v001/openfoam-extended/solver-continuation-to-4s.log');rows=[r for r in residual_steps(p.read_text(),p.name) if r['time_s']>3.95+1e-8];rows[-1]['variables']['p']['terminal_initial']=1.01e-5
  with self.assertRaises(ValueError):verify_residual_window(rows,3.95,10)
if __name__=='__main__':unittest.main()