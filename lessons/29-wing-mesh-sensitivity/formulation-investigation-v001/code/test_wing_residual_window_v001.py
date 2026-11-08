import wing_import_paths_v001
import unittest,numpy as np
from wing_residual_window_helpers_v001 import declared_iterations,verify_runtime
class Window(unittest.TestCase):
 def sample(self):
  centres=np.array([[0.,0,0],[1.,0,0],[2.,0,0]]);u=np.array([[1.,0,0],[1.,.1,0],[1.,0,.1]]);p=np.array([.1,.2,.3]);order=np.array([2,0,1]);r=np.array([1e-5,-2e-5,3e-5]);arr={'SV_CENTROID':centres[order].ravel(),'SV_U':u[order,0],'SV_V':u[order,1],'SV_W':u[order,2],'SV_P':p[order]}
  for k in ['SV_P_R','SV_U_R','SV_V_R','SV_W_R']:arr[k]=r[order]
  return centres,u,p,arr,{k:6e-5 for k in ['continuity','x-velocity','y-velocity','z-velocity']}
 def test_finite_global_window(self):
  self.assertEqual(declared_iterations(1502,1510),list(range(1503,1511)))
  with self.assertRaises(ValueError):declared_iterations(1502,1511)
 def test_exact_runtime_identity(self):
  c,u,p,a,r=self.sample();v=verify_runtime(c,u,p,a,r,1e-8);np.testing.assert_equal(v['residuals']['SV_P_R'],[1e-5,-2e-5,3e-5]);self.assertEqual(v['stage'],'assembled_initial')
 def test_tampered_pressure_refused(self):
  c,u,p,a,r=self.sample();a['SV_P'][0]+=.1
  with self.assertRaises(ValueError):verify_runtime(c,u,p,a,r,1e-8)
 def test_wrong_residual_sum_refused(self):
  c,u,p,a,r=self.sample();a['SV_P_R'][0]*=2
  with self.assertRaises(ValueError):verify_runtime(c,u,p,a,r,1e-8)
 def test_wrong_stage_not_substituted(self):
  c,u,p,a,r=self.sample();a['SV_P_R']*=1e-3;r['continuity']=6e-8
  with self.assertRaises(ValueError):verify_runtime(c,u,p,a,r,6e-8)
if __name__=='__main__':unittest.main()
