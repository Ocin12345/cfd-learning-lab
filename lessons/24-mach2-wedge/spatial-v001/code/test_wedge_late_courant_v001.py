import unittest,numpy as np
from wedge_late_courant_v001 import record_courant
class CFL(unittest.TestCase):
 def test_inherited_startup_and_scaled_first_advance(self):
  a=record_courant([.3,.3,.15,.15],[.5e-6]*3,[.0070005,.007001,.0070015],1e-6)
  self.assertEqual(a['raw_startup_courant'],[.3,.3]);self.assertAlmostEqual(a['maximum_advance_scaled_courant'],.15)
  self.assertEqual(a['records'][2]['pre_advance_time_s'],.0070005)
 def test_actual_large_advance_rejected(self):
  with self.assertRaises(ValueError):record_courant([.3,.3,.15],[1e-6,.5e-6],[.007001,.0070015],1e-6)
 def test_large_later_native_record_rejected(self):
  with self.assertRaises(ValueError):record_courant([.3,.3,.16],[.5e-6]*2,[.0070005,.007001],1e-6)
if __name__=='__main__':unittest.main()