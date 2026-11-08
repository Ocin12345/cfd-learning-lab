import os,subprocess,sys,unittest
from pathlib import Path
W=Path(__file__).parent
class Deployment(unittest.TestCase):
 def test_native_import_with_no_inherited_pythonpath(self):
  e=os.environ.copy();e.pop('PYTHONPATH',None)
  p=subprocess.run([sys.executable,'-B','-c','import wing_import_paths_v001; import wing_repair_native_v003; print("native_import_ok")'],cwd=W,env=e,capture_output=True,text=True)
  self.assertEqual(p.returncode,0,p.stderr);self.assertIn('native_import_ok',p.stdout)
if __name__=='__main__':unittest.main()
