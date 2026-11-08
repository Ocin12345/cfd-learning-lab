import wing_import_paths_v001
import unittest,numpy as np
from wing_retained_residual_helpers_v001 import expected_probe_settings,classify_stage
from wing_repair_native_v003 import expected_settings
class Retention(unittest.TestCase):
 def test_only_storage_flag_added(self):
  a=expected_settings(.1,20.,.5);a['controls']['advanced']['expert']['save_cell_residuals']=True;self.assertEqual(a,expected_probe_settings())
 def test_assembled_stage(self):self.assertEqual(classify_stage(np.array([2e-5,-3e-5]),5e-5,1e-8),'assembled_initial')
 def test_corrected_stage(self):self.assertEqual(classify_stage(np.array([2e-9,-3e-9]),5e-5,5e-9),'corrected_final')
 def test_unknown_stage_refused(self):
  with self.assertRaises(ValueError):classify_stage(np.array([1.]),5e-5,1e-8)
 def test_ambiguous_stage_refused(self):
  with self.assertRaises(ValueError):classify_stage(np.array([5e-5]),5e-5,5e-5)
 def test_nonfinite_refused(self):
  with self.assertRaises(ValueError):classify_stage(np.array([np.nan]),5e-5,1e-8)
if __name__=='__main__':unittest.main()
