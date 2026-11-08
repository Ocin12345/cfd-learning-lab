import unittest
from unittest.mock import patch
import sloshing_wait_admission_v001 as g
class Wait(unittest.TestCase):
 def test_sources_rechecked_after_resource_wait(self):
  with patch.object(g,'verify',side_effect=[1,AssertionError('changed source')]) as v,patch.object(g,'available_memory',return_value=0),patch.object(g.time,'sleep'),patch.object(g,'prelaunch') as p:
   with self.assertRaisesRegex(AssertionError,'changed source'):g.wait_admission('x.json')
   self.assertEqual(v.call_count,2);p.assert_not_called()
 def test_transient_gate_resource_race_retries_all_checks(self):
  with patch.object(g,'verify',return_value=1) as v,patch.object(g,'available_memory',return_value=8*2**30),patch.object(g.shutil,'disk_usage',return_value=type('D',(),{'free':20*2**30})()),patch.object(g.time,'sleep'),patch.object(g,'prelaunch',side_effect=[AssertionError(),{'verified_inputs':1}]) as p:
   self.assertEqual(g.wait_admission('x.json'),{'verified_inputs':1});self.assertEqual(v.call_count,2);self.assertEqual(p.call_count,2)
if __name__=='__main__':unittest.main()
