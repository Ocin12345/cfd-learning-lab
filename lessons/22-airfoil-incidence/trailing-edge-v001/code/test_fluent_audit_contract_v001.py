import unittest
from audit_refined_fluent_v001 import applicable_checks, check_targets
class Contract(unittest.TestCase):
 def test_five_degree_symmetry_is_inapplicable_without_hiding_residual_miss(self):
  a=applicable_checks({'fluent.residuals':False,'fluent.Cd_stability':True,'fluent.zero_angle_lift_diagnostic':False},5)
  self.assertEqual(a,{'fluent.residuals':False,'fluent.Cd_stability':True})
  self.assertFalse(all(a.values()))
 def test_wrong_angle_rejected(self):
  with self.assertRaises(ValueError):applicable_checks({'fluent.residuals':True},0)
 def test_weaker_targets_rejected(self):
  t=dict(residual_target=1e-7,relative_mass_imbalance_max=1e-6,drag_relative_range_max=1e-4,lift_absolute_range_max=1e-5,minimum_orthogonal_quality=.2,last_window_samples=20,additional_iterations=100)
  check_targets(t)
  t['residual_target']=1e-6
  with self.assertRaises(ValueError):check_targets(t)
if __name__=='__main__':unittest.main()
