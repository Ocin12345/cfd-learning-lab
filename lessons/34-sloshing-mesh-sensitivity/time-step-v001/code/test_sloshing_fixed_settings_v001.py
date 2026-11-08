import unittest,copy
from sloshing_fixed_settings_v001 import compare_fixed
class Settings(unittest.TestCase):
 def test_same_and_serialization_roundoff(self):
  a={'methods':{'piso':True,'tol':1e-10},'material':[998.2,.001003]};b=copy.deepcopy(a);b['material'][0]*=1+1e-15;self.assertTrue(compare_fixed(a,b)['passed'])
 def test_reject_method_material_or_extra_key(self):
  a={'method':'PISO','rho':998.2,'flag':False}
  for k,v in [('method','SIMPLE'),('rho',998.3),('flag',True),('extra',1)]:
   b=dict(a);b[k]=v
   with self.assertRaises(AssertionError):compare_fixed(a,b)
if __name__=='__main__':unittest.main()
