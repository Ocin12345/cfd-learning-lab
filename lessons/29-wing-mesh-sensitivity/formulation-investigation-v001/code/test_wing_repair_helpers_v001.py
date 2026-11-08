import unittest
import numpy as np
import wing_repair_helpers_v001 as h
class Helpers(unittest.TestCase):
 def api(self,n):self.assertTrue(hasattr(h,n),n+' missing');return getattr(h,n)
 def test_exact_monitor_echo_and100newrows(self):
  rows=['1500 .23 .30']+[f'{i} .23 .30' for i in range(1501,1601)];r=self.api('force_window')('header\n'+'\n'.join(rows),1501,1600,dict(iteration=1500,Cd=.23,Cl=.30));self.assertEqual(len(r),100);self.assertEqual(r[0]['iteration'],1501)
 def test_conflicting_or_repeated_rows_rejected(self):
  f=self.api('force_window')
  for text in ['1500 .24 .30\n1501 .23 .30','1501 .23 .30\n1501 .23 .30']:
   with self.assertRaises(ValueError):f(text,1501,1501,dict(iteration=1500,Cd=.23,Cl=.30))
 def test_missing_and_nonfinite_monitor_rejected(self):
  f=self.api('force_window')
  for text in ['1501 .23 .30','1501 nan .30\n1502 .23 .30']:
   with self.assertRaises(ValueError):f(text,1501,1502,None)
 def test_global_residual_window_and_four_equations(self):
  q={k:np.ones(3)*1e-7 for k in ['continuity','x-velocity','y-velocity','z-velocity']};r=self.api('residual_window')([1500,1501,1502],q,1501,1502);self.assertEqual([v['iteration'] for v in r],[1501,1502])
 def test_missing_z_or_wrong_global_coverage_rejected(self):
  f=self.api('residual_window');q={k:np.ones(2)*1e-7 for k in ['continuity','x-velocity','y-velocity','z-velocity']}
  with self.assertRaises(ValueError):f([1,2],q,1501,1502)
  del q['z-velocity']
  with self.assertRaises(ValueError):f([1501,1502],q,1501,1502)
 def test_residual_prefix_numerator_and_scale_preserved(self):
  f=self.api('same_residual_prefix');p={'continuity':np.array([[1.,2.,0.,0.],[.1,2.,0.,0.]])};n={'continuity':np.vstack([p['continuity'],[.01,2.,0.,0.]])};self.assertTrue(f([1,2],p,[1,2,3],n,2));n['continuity'][0,1]=3
  with self.assertRaises(ValueError):f([1,2],p,[1,2,3],n,2)
 def test_original_threshold_read_without_changing_flag(self):
  class Toggle:
   def __init__(self):self.v=False
   def __call__(self):return self.v
   def set_state(self,x):self.v=x
  class Eq:
   def __init__(self):self.check_convergence=Toggle()
   def absolute_criteria(self):
    if not self.check_convergence():raise RuntimeError('inactive')
    return 1e-5
  e=Eq();self.assertEqual(self.api('original_target')(e),1e-5);self.assertFalse(e.check_convergence())
 def test_target_getter_restores_flag_on_failure(self):
  class Toggle:
   def __init__(self):self.v=False
   def __call__(self):return self.v
   def set_state(self,x):self.v=x
  class Eq:
   def __init__(self):self.check_convergence=Toggle()
   def absolute_criteria(self):raise RuntimeError('read error')
  e=Eq()
  with self.assertRaisesRegex(RuntimeError,'read error'):self.api('original_target')(e)
  self.assertFalse(e.check_convergence())
if __name__=='__main__':unittest.main()
