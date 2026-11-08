import wing_import_paths_v001
import unittest,numpy as np,tempfile,json
from pathlib import Path
from types import SimpleNamespace
from wing_simplec_qualification_helpers_v001 import checkpoints,scale_guard,query_iteration
from wing_repair_common_v001 import R
class Expressions(dict):
 def get_state(self):return {k:v.definition for k,v in self.items()}
 def get_object_names(self):return list(self)
 def create(self,name):self[name]=SimpleNamespace(definition='',get_value=lambda:1520.)
class Qualification(unittest.TestCase):
 def test_finite100_schedule(self):
  self.assertEqual(checkpoints(),[1520,1540,1560,1580,1600]);self.assertEqual(sum(np.diff([1500,*checkpoints()])),100)
 def test_continuity_scale_preserved(self):
  i=np.arange(1,1601);v=np.zeros((1600,4));v[:,1]=2.2344914086503613;self.assertTrue(scale_guard(i,v,2.2344914086503613,1600));v[1590,1]*=2
  with self.assertRaises(ValueError):scale_guard(i,v,2.2344914086503613,1600)
 def test_missing_iteration_refused(self):
  with self.assertRaises(ValueError):scale_guard(np.array([1501,1503]),np.ones((2,4)),1.,1503)
 def test_live_query_writes_only_new_directory_and_restores(self):
  with tempfile.TemporaryDirectory(dir=R/'temp') as d:
   ex=Expressions();st=SimpleNamespace(setup=SimpleNamespace(named_expressions=ex));self.assertEqual(query_iteration(st,Path(d),'actual1520',1520),1520);self.assertEqual(ex,{});self.assertEqual(json.loads((Path(d)/'actual1520-live-iteration.json').read_text())['query_value'],1520.)
   with self.assertRaises(FileExistsError):query_iteration(st,Path(d),'actual1520',1520)
 def test_bad_live_counter_restores_temporary_expression(self):
  with tempfile.TemporaryDirectory(dir=R/'temp') as d:
   ex=Expressions();st=SimpleNamespace(setup=SimpleNamespace(named_expressions=ex))
   with self.assertRaises(ValueError):query_iteration(st,Path(d),'wrong',1540)
   self.assertEqual(ex,{})
if __name__=='__main__':unittest.main()
