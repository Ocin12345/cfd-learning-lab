import unittest
from cavity_step_audit_v001 import audit_steps
class Steps(unittest.TestCase):
 def block(self,n,t,res='1e-7 2e-7 3e-7 1e-10'):
  return f'Updating solution at time level N\n 3 {res} 0:00:01 10\nFlow time = {t}s, time step = {n}\n'
 def test_complete_and_terminal_limits(self):
  a=audit_steps(self.block(1,.01)+self.block(2,.02),2,.01);self.assertTrue(all(r['passed'] for r in a));self.assertEqual(a[-1]['step'],2)
 def test_actual_miss_preserved(self):
  a=audit_steps(self.block(1,.01,'2e-6 2e-7 3e-7 1e-10'),1,.01);self.assertFalse(a[0]['passed'])
 def test_missing_or_wrong_time_rejected(self):
  for text in [self.block(2,.02),self.block(1,.02)]:
   with self.assertRaises(ValueError):audit_steps(text,1,.01)
if __name__=='__main__':unittest.main()