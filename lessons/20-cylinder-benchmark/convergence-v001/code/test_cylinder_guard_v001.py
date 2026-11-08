import hashlib,json,tempfile,unittest
from pathlib import Path
from cylinder_guard_v001 import verify_inputs
class GuardTests(unittest.TestCase):
 def test_fresh_inputs_and_fresh_output(self):
  with tempfile.TemporaryDirectory(prefix='cylinder-guard-') as d:
   p=Path(d);f=p/'source.py';f.write_text('unchanged');m=p/'inputs.json';m.write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest()}));self.assertEqual(verify_inputs(m,p/'new'),1)
 def test_changed_source_rejected(self):
  with tempfile.TemporaryDirectory(prefix='cylinder-guard-') as d:
   p=Path(d);f=p/'source.py';f.write_text('old');m=p/'inputs.json';m.write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest()}));f.write_text('changed')
   with self.assertRaises(ValueError):verify_inputs(m,p/'new')
 def test_reused_output_rejected(self):
  with tempfile.TemporaryDirectory(prefix='cylinder-guard-') as d:
   p=Path(d);m=p/'inputs.json';m.write_text('{}');out=p/'existing';out.mkdir()
   with self.assertRaises(ValueError):verify_inputs(m,out)
if __name__=='__main__':unittest.main()