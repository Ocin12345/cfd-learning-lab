import unittest
from prepare_airfoil_continuation_v001 import rewrite_control
class Restart(unittest.TestCase):
 def test_only_top_level_restart_limits_change(self):
  raw='startFrom startTime; startTime 0; endTime 3000; writeInterval 3000; functions { f { writeInterval 1; type forceCoeffs; } }'
  new=rewrite_control(raw,3000,9000)
  self.assertIn('endTime 9000;',new);self.assertIn('writeInterval 1;',new);self.assertIn('startFrom latestTime;',new)
 def test_bad_prior_or_limit_rejected(self):
  raw='startFrom startTime; startTime 0; endTime 3000; writeInterval 3000; functions { }'
  for prior,end in ((2000,9000),(3000,3000)):
   with self.assertRaises(ValueError):rewrite_control(raw,prior,end)
if __name__=='__main__':unittest.main()
