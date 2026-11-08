import wing_import_paths_v001
import unittest,numpy as np
from wing_simplec_runtime_v001 import verify_simplec_runtime
class Runtime(unittest.TestCase):
 def fixture(self):
  c=np.array([[0.,0,0],[1.,0,0]]);u=np.zeros((2,3));p=np.ones(2);a={'SV_CENTROID':c.ravel(),'SV_U':u[:,0],'SV_V':u[:,1],'SV_W':u[:,2],'SV_P':p}
  for k in ['SV_P_R','SV_U_R','SV_V_R','SV_W_R']:a[k]=np.array([.2,.3])
  return c,u,p,a
 def test_unknown_stage_is_recorded_not_called_coupled(self):
  c,u,p,a=self.fixture();r=verify_simplec_runtime(c,u,p,a,dict.fromkeys(['continuity','x-velocity','y-velocity','z-velocity'],3.),4.)
  self.assertEqual(r['pressure_stage'],'unidentified');self.assertFalse(r['raw_numerator_matches']['SV_P_R'])
 def test_raw_match_only_labels_observed_numerator(self):
  c,u,p,a=self.fixture();r=verify_simplec_runtime(c,u,p,a,dict.fromkeys(['continuity','x-velocity','y-velocity','z-velocity'],.5),4.)
  self.assertEqual(r['pressure_stage'],'matches_HDF_numerator_only');self.assertTrue(all(r['raw_numerator_matches'].values()))
 def test_runtime_field_or_identity_mismatch_rejected(self):
  c,u,p,a=self.fixture();a['SV_P']=np.array([8.,1.])
  with self.assertRaises(ValueError):verify_simplec_runtime(c,u,p,a,dict.fromkeys(['continuity','x-velocity','y-velocity','z-velocity'],.5),4.)
if __name__=='__main__':unittest.main()
