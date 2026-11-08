import unittest
import numpy as np
from shock_width_v001 import connected_width
class Width(unittest.TestCase):
 def test_known_linear_transition(self):
  d=np.linspace(-.1,.1,41);q=np.clip((d+.02)/.04,0,1);r=connected_width(d,q);self.assertAlmostEqual(r['width_m'],.032,places=12);self.assertAlmostEqual(r['midpoint_m'],0,places=12)
 def test_descending_input_is_equivalent(self):
  d=np.linspace(-.1,.1,41);q=np.clip((d+.02)/.04,0,1);a=connected_width(d,q);b=connected_width(d[::-1],q[::-1]);self.assertAlmostEqual(a['width_m'],b['width_m'])
 def test_missing_bracket_rejected(self):
  with self.assertRaises(ValueError):connected_width(np.linspace(0,1,20),np.linspace(.2,.8,20))
 def test_does_not_bridge_unrelated_shocks(self):
  d=np.arange(9,dtype=float);q=np.array([0,.1,.25,.75,.95,.9,.3,.4,.5]);r=connected_width(d,q);self.assertLess(r['width_m'],5);self.assertLess(r['upper_bracket_indices'][1],6)
 def test_disconnected_threshold_rejected(self):
  with self.assertRaises(ValueError):connected_width(np.arange(8,dtype=float),[0,.1,.4,.8,.2,.4,.6,.95])
if __name__=='__main__':unittest.main()
