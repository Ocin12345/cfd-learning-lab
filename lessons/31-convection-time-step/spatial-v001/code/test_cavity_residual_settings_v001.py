import unittest
from cavity_residual_settings_v001 import configure_residuals
class Equation:
 def __init__(self):self.check_convergence=False;self._criteria=None
 @property
 def absolute_criteria(self):return self._criteria
 @absolute_criteria.setter
 def absolute_criteria(self,v):
  if not self.check_convergence:raise RuntimeError('Native setting inactive')
  self._criteria=v
class Container(dict):
 def get_object_names(self):return list(self)
class Residual:
 def __init__(self):self.equations=Container({k:Equation() for k in ('continuity','x-velocity','y-velocity','energy')})
class Restart(unittest.TestCase):
 def test_loaded_disabled_control_restored_after_setting(self):
  r=Residual();configure_residuals(r,False)
  self.assertTrue(all(not e.check_convergence for e in r.equations.values()));self.assertEqual(r.equations['energy'].absolute_criteria,1e-9);self.assertEqual(r.equations['continuity'].absolute_criteria,1e-6)
 def test_heat_checks_enabled(self):
  r=Residual();configure_residuals(r,True);self.assertTrue(all(e.check_convergence for e in r.equations.values()))
if __name__=='__main__':unittest.main()