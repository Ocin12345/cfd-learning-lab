import unittest
from prepare_finer_extension_v001 import validate_miss,rewrite_control
class Extension(unittest.TestCase):
 def test_24000_residual_only_miss_admitted(self):
  validate_miss(dict(iteration=24000,source_case='case',numerical_checks_passed=False,failed_checks=['openfoam.residuals']),'case',24000)
 def test_mass_miss_or_wrong_checkpoint_not_admitted(self):
  for failures,iteration in [(['openfoam.mass_balance'],24000),(['openfoam.residuals'],3000)]:
   with self.assertRaises(ValueError):validate_miss(dict(iteration=iteration,source_case='case',numerical_checks_passed=False,failed_checks=failures),'case',24000)
 def test_only_execution_controls_change(self):
  old='startFrom latestTime; startTime 3000; endTime 24000; writeInterval 24000; functions { forceCoefficients { writeInterval 1; type forceCoeffs; } }'
  s=rewrite_control(old,24000,32000)
  self.assertIn('startTime 24000;',s);self.assertIn('endTime 32000;',s);self.assertIn('writeInterval 1;',s)
if __name__=='__main__':unittest.main()
