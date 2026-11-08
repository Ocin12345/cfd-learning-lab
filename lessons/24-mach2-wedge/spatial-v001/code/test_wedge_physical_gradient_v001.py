import unittest,numpy as np
from wedge_physical_gradient_v001 import gradient_magnitude
class Gradient(unittest.TestCase):
 def test_linear_physical_field_on_skew_grid(self):
  i,j=np.meshgrid(np.linspace(0,1,37),np.linspace(0,1,29));x=i+.1*j*j;y=j+.2*i
  np.testing.assert_allclose(gradient_magnitude(2*x-3*y,x,y),np.sqrt(13),atol=1e-11,rtol=0)
 def test_degenerate_grid_rejected(self):
  i,j=np.meshgrid(np.arange(5.),np.arange(5.))
  with self.assertRaises(ValueError):gradient_magnitude(i,i,2*i)
if __name__=='__main__':unittest.main()