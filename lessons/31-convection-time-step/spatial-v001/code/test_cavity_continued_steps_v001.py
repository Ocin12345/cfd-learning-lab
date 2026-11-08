import unittest
from cavity_continued_steps_v001 import audit_continued_steps
class Continued(unittest.TestCase):
 def block(self,n,t,r='1e-7 2e-7 3e-7 1e-10'):return f'Updating solution at time level N\n 3001 {r} 0:00:01 199\nFlow time = {t}s, time step = {n}\n'
 def test_global_steps_not_reset(self):
  a=audit_continued_steps(self.block(101,1.01)+self.block(102,1.02),100,102,.01);self.assertEqual([r['step'] for r in a],[101,102]);self.assertTrue(all(r['passed'] for r in a))
 def test_no_duplicate_prefix_or_gap(self):
  for text in [self.block(100,1),self.block(102,1.02),self.block(101,1.02)]:
   with self.assertRaises(ValueError):audit_continued_steps(text,100,101,.01)
 def test_miss_is_not_pass(self):
  a=audit_continued_steps(self.block(101,1.01,'1e-7 2e-6 3e-7 1e-10'),100,101,.01);self.assertFalse(a[0]['passed'])
if __name__=='__main__':unittest.main()