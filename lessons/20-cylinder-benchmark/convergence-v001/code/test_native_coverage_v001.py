import unittest
from pathlib import Path
from analyse_vortex_fluent_v003 import terminal_rows
from audit_cylinder_window_v001 import verify_step_rows
class NativeCoverage(unittest.TestCase):
 def test_actual_parent_window_and_missing_step(self):
  p=Path('E:/CFD-Learning-Lab/native/runs/2026-11-01-vortex-shedding-v001/fluent-corrected-v003/solver-transcript.txt');r=[a for a in terminal_rows(p.read_text()) if a['time_s']>3.95+1e-8]
  self.assertEqual(len(verify_step_rows(r,3.95,10)),10)
  with self.assertRaises(ValueError):verify_step_rows(r[:-1],3.95,10)
if __name__=='__main__':unittest.main()