import unittest,copy
from cavity_restart_settings_v001 import compare_restart_settings
KEYS=['general','energy','viscous','material','boundaries','references','methods','controls']
class Reference(unittest.TestCase):
 def fixture(self):
  a={k:{} for k in KEYS};a['methods']={'p_v_coupling':{'flow_scheme':'PISO'}};a['references']={'zone':'fluid','depth':1.};a['material']={'conductivity':.025478873239436622};return a
 def test_checkpoint_zone_is_enforced(self):
  a=self.fixture();old=copy.deepcopy(a);old['references']['zone']=''
  with self.assertRaises(AssertionError):compare_restart_settings(a,old)
  self.assertTrue(compare_restart_settings(a,a)['passed'])
 def test_only_original_serialization_roundoff(self):
  a=self.fixture();b=copy.deepcopy(a);b['material']['conductivity']=.02547887323943662
  r=compare_restart_settings(a,b);self.assertEqual(len(r['roundoff_only']),1)
 def test_physical_reference_change_rejected(self):
  a=self.fixture();b=copy.deepcopy(a);b['references']['depth']=1.001
  with self.assertRaises(AssertionError):compare_restart_settings(a,b)
if __name__=='__main__':unittest.main()