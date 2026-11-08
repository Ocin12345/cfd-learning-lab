import unittest
import numpy as np
import jet_energy_restart_v001 as r
class Restart(unittest.TestCase):
 def api(self,n):self.assertTrue(hasattr(r,n),n+' missing');return getattr(r,n)
 def rows(self,start,end):return [dict(step=i,time_s=i*.005,residual=1e-10) for i in range(start,end+1)]
 def test_exact_two_segment_coverage(self):
  a=self.api('merge_steps')(self.rows(1,20),self.rows(21,600));self.assertEqual(len(a),600);self.assertEqual(a[-1]['step'],600)
 def test_restart_duplicate_is_rejected(self):
  with self.assertRaises(ValueError):self.api('merge_steps')(self.rows(1,20),self.rows(20,600))
 def test_missing_terminal_step_rejected(self):
  with self.assertRaises(ValueError):self.api('merge_steps')(self.rows(1,20),self.rows(21,599))
 def test_clock_reset_or_wrong_split_rejected(self):
  a=self.rows(21,600);a[0]['time_s']=.005
  with self.assertRaises(ValueError):self.api('merge_steps')(self.rows(1,20),a)
  with self.assertRaises(ValueError):self.api('merge_steps')(self.rows(1,19),self.rows(20,600))
 def test_exact_restart_arrays(self):
  a={'T':np.array([300.,310.]),'phi':np.array([-.1,.1])};self.assertTrue(self.api('same_restart')(a,{k:v.copy() for k,v in a.items()}))
 def test_any_restart_field_change_rejected(self):
  f=self.api('same_restart');a={'T':np.array([300.,310.]),'phi':np.array([-.1,.1])};b={k:v.copy() for k,v in a.items()};b['T'][0]+=1e-10
  with self.assertRaises(ValueError):f(a,b)
  with self.assertRaises(ValueError):f(a,{'T':a['T']})
  with self.assertRaises(ValueError):f({'T':[float('nan')]},{'T':[float('nan')]})
if __name__=='__main__':unittest.main()
