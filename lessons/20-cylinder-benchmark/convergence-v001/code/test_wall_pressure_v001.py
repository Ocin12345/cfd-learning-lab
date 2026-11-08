import unittest
import numpy as np
from cylinder_wall_pressure_v001 import angular_stencil
class WallPressureStencil(unittest.TestCase):
 def test_linear_angular_pressure_at_rear(self):
  theta=np.array([-.02,.03]);xy=np.c_[.2+.05*np.cos(theta),.2+.05*np.sin(theta)];s=angular_stencil(xy,0.);p=2+3*theta;self.assertAlmostEqual(float(np.dot(p[s['indices']],s['weights'])),2,places=12)
 def test_front_wrap_around_pi(self):
  delta=np.array([-.02,.03]);theta=np.pi+delta;xy=np.c_[.2+.05*np.cos(theta),.2+.05*np.sin(theta)];s=angular_stencil(xy,np.pi);p=5+2*delta;self.assertAlmostEqual(float(np.dot(p[s['indices']],s['weights'])),5,places=12)
 def test_one_sided_stencil_rejected(self):
  theta=np.array([.02,.03]);xy=np.c_[.2+.05*np.cos(theta),.2+.05*np.sin(theta)]
  with self.assertRaises(ValueError):angular_stencil(xy,0.)
if __name__=='__main__':unittest.main()