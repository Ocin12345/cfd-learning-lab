import unittest
from unittest.mock import Mock
import jet_energy_parent_gate_v001 as g
class Gate(unittest.TestCase):
 def api(self):self.assertTrue(hasattr(g,'wait_checked'));return g.wait_checked
 def test_parent_rechecked_after_wait_and_before_return(self):
  events=[];parent=lambda:events.append('parent');n=[0]
  def gate():
   events.append('gate');n[0]+=1
   if n[0]==1:raise BlockingIOError('RAM')
   return {'accepted':True}
  result=self.api()(parent,gate,lambda x:events.append('wait'))
  self.assertTrue(result['accepted']);self.assertEqual(events,['parent','gate','wait','parent','gate','parent','gate'])
 def test_parent_change_while_waiting_is_fatal(self):
  parent=Mock(side_effect=[None,RuntimeError('parent changed')]);gate=Mock(side_effect=BlockingIOError('RAM'))
  with self.assertRaisesRegex(RuntimeError,'parent changed'):self.api()(parent,gate,lambda x:None)
  self.assertEqual(gate.call_count,1)
 def test_resource_race_rechecks_parent(self):
  parent=Mock();gate=Mock(side_effect=[{},BlockingIOError('RAM race'),{},{}]);self.api()(parent,gate,lambda x:None);self.assertEqual(parent.call_count,4)
if __name__=='__main__':unittest.main()
