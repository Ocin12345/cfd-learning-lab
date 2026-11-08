import unittest,json
from pathlib import Path
from sloshing_openfoam_steps_v001 import parse_steps
R=Path('E:/CFD-Learning-Lab')
class Steps(unittest.TestCase):
 def test_actual_baseline_residuals_and_alpha(self):
  r=parse_steps((R/'native/runs/sloshing-refinement-v001/production/openfoam/solver.log').read_text(),0,2400,.0005)
  q=json.loads((R/'local/sloshing-refinement-v001/openfoam-production-audit.json').read_text())['step_records'];self.assertEqual(r,q)
 def test_global_time_missing_corrector_and_alpha_miss(self):
  t='Time = .02025s\nPIMPLE: Iteration 8\n'
  for k in ['Ux','Uy','p_rgh']:t+=f'Solving for {k}, Initial residual = 1e-6, Final residual = 1e-9\n'
  t+='Min(alpha.water) = 0 Max(alpha.water) = 1\n'
  self.assertTrue(parse_steps(t,80,81,.00025)[0]['passed'])
  for v in [t.replace('.02025s','.02s'),t.replace('Iteration 8','Iteration 7')]:
   with self.assertRaises(AssertionError):parse_steps(v,80,81,.00025)
  self.assertFalse(parse_steps(t.replace('= 1\n','= 1.01\n'),80,81,.00025)[0]['passed'])
if __name__=='__main__':unittest.main()
