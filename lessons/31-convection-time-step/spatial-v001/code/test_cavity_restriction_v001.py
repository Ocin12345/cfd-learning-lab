import unittest,numpy as np
from cavity_restriction_v001 import restrict_square,compare_fields
class Restriction(unittest.TestCase):
 def test_affine_scalar_and_vector_locations(self):
  n=7;x=(np.arange(2*n)+.5)/(2*n);xx,yy=np.meshgrid(x,x);fine=(2*xx+3*yy).reshape(-1);u=np.stack([-yy,xx],axis=-1).reshape(-1,2)
  c=(np.arange(n)+.5)/n;xc,yc=np.meshgrid(c,c)
  np.testing.assert_allclose(restrict_square(fine,n),(2*xc+3*yc).reshape(-1),atol=1e-14)
  np.testing.assert_allclose(restrict_square(u,n),np.stack([-yc,xc],axis=-1).reshape(-1,2),atol=1e-14)
 def test_conserves_integral_for_arbitrary_field(self):
  q=np.sin(np.arange(144.));coarse=restrict_square(q,6);self.assertAlmostEqual(q.mean(),coarse.mean(),places=14)
 def test_invalid_shape_and_nonfinite_rejected(self):
  for q,n in [(np.ones(15),2),(np.ones(16),0),(np.full(16,np.nan),2)]:
   with self.assertRaises(ValueError):restrict_square(q,n)
 def test_zero_start_absolute_only(self):
  c=dict(T=np.ones(4)*300,U=np.zeros((4,2)));f=dict(T=np.ones(16)*300,U=np.zeros((16,2)))
  a=compare_fields(c,f,2,relative=False);self.assertEqual(a['max_temperature_difference_K'],0);self.assertIsNone(a['velocity_relative_L2'])
if __name__=='__main__':unittest.main()