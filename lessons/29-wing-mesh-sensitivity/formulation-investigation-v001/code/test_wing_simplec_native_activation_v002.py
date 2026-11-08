import wing_import_paths_v001
import unittest,json,copy
from wing_repair_common_v001 import N
from wing_simplec_capability_helpers_v002 import validate_state
class NativeActivation(unittest.TestCase):
 def pair(self):
  p=N/'fluent-simplec-capability-v001';return json.loads((p/'loaded-source-settings.json').read_text()),json.loads((p/'candidate-settings.json').read_text())
 def test_actual_native_default_activation(self):
  a,b=self.pair();self.assertTrue(validate_state(a,b))
 def test_enhanced_correction_cannot_be_enabled(self):
  a,b=self.pair();b['controls']['advanced']['enhanced_skewness_correction']=True
  with self.assertRaises(ValueError):validate_state(a,b)
 def test_native_momentum_default_cannot_change(self):
  a,b=self.pair();b['controls']['advanced']['multi_grid']['mg_controls']['x-mom']['termination_criteria']=.2
  with self.assertRaises(ValueError):validate_state(a,b)
 def test_extra_equation_not_admitted(self):
  a,b=self.pair();b['controls']['advanced']['multi_grid']['mg_controls']['energy']={'termination_criteria':.1}
  with self.assertRaises(ValueError):validate_state(a,b)
if __name__=='__main__':unittest.main()
