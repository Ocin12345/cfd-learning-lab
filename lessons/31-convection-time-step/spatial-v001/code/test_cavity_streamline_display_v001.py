import unittest
import numpy as np
from cavity_streamline_display_v001 import streamline_display
class Display(unittest.TestCase):
 def test_stagnation_mask_and_unchanged_native_rotation(self):
  c=np.arange(5,dtype=float);x,y=np.meshgrid(c,c);u=np.stack((-(y-2),x-2),axis=-1);old=u.copy();seed=np.array([[2.,2.],[1.,2.],[4.,2.]])
  a,b,s=streamline_display(u,c,seed,.1)
  self.assertTrue(a.mask[2,2]);self.assertTrue(b.mask[2,2]);np.testing.assert_array_equal(s,seed[1:]);np.testing.assert_array_equal(u,old)
 def test_quiet_field_has_no_seeds(self):
  a,b,s=streamline_display(np.zeros((3,3,2)),np.arange(3),np.array([[1.,1.]]),.1)
  self.assertTrue(a.mask.all());self.assertEqual(s.shape,(0,2))
if __name__=='__main__':unittest.main()
