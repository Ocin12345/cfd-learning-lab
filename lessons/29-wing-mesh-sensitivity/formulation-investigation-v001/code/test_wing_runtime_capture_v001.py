import unittest
from types import SimpleNamespace
from wing_runtime_residual_export_v001 import scalar_names,select_svars
class Capture(unittest.TestCase):
 def test_public_runtime_catalogue(self):
  s=SimpleNamespace(fields=SimpleNamespace(field_data=SimpleNamespace(scalar_fields=SimpleNamespace(allowed_values=lambda:['continuity-residual']))));self.assertEqual(scalar_names(s),['continuity-residual'])
 def test_capture_base_fields_and_actual_residuals(self):
  base=['SV_CENTROID','SV_U','SV_V','SV_W','SV_P'];self.assertEqual(select_svars(base+['SV_P_R','SV_U_R','SV_IGNORED']),base+['SV_P_R','SV_U_R'])
 def test_missing_residual_is_not_success(self):
  with self.assertRaises(ValueError):select_svars(['SV_CENTROID','SV_U','SV_V','SV_W','SV_P'])
if __name__=='__main__':unittest.main()
