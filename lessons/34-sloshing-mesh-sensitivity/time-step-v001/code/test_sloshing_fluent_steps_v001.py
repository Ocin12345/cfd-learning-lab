import unittest,json
from pathlib import Path
from sloshing_fluent_steps_v001 import parse_steps
R=Path('E:/CFD-Learning-Lab')
class Steps(unittest.TestCase):
 def test_actual_baseline_inner_counts(self):
  r=parse_steps((R/'native/runs/sloshing-refinement-v001/production/fluent/solver-transcript.txt').read_text(errors='replace'),0,2400,.0005)
  expected=json.loads((R/'local/sloshing-refinement-v001/fluent-production-audit.json').read_text())['step_records']
  self.assertEqual([x['inner_iterations'] for x in r],[x['inner_iterations'] for x in expected]);self.assertTrue(all(x['passed'] for x in r))
 def test_global_continuation_and_rejections(self):
  t='   900  1e-5 1e-5 1e-5  0:00:01 400\n   906  1e-5 1e-5 1e-5  0:00:01 394\nFlow time = .02025s, time step = 81\n'
  r=parse_steps(t,80,81,.00025);self.assertEqual(r[0]['inner_iterations'],6)
  for v in [t.replace('step = 81','step = 1'),t.replace('.02025s','.02s'),t.replace('394','0').replace('906','1301')]:
   with self.assertRaises(AssertionError):parse_steps(v,80,81,.00025)
if __name__=='__main__':unittest.main()
