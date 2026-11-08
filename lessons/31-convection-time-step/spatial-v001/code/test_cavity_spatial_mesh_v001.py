import unittest,numpy as np,uuid,hashlib
from pathlib import Path
from cavity_spatial_mesh_v001 import generate,refinement_proof
class Geometry(unittest.TestCase):
 def path(self):return Path('E:/CFD-Learning-Lab/temp')/('cavity-spatial-test-'+uuid.uuid4().hex)
 def test_parent_partition_and_nodes(self):
  a,b=self.path(),self.path();generate(a,4);generate(b,8)
  with np.load(a/'geometry.npz') as c,np.load(b/'geometry.npz') as f:
   q=refinement_proof(c,f);self.assertEqual(q['old_node_max_error_m'],0);self.assertLess(q['maximum_parent_area_error_m2'],1e-17)
 def test_exact_original48_mesh(self):
  p=self.path();generate(p,48);old=Path('E:/CFD-Learning-Lab/native/meshes/natural-convection-cavity-48-v001/cavity.msh')
  self.assertEqual(hashlib.sha256((p/'cavity.msh').read_bytes()).hexdigest(),hashlib.sha256(old.read_bytes()).hexdigest())
 def test_refuse_existing_directory(self):
  p=self.path();p.mkdir()
  with self.assertRaises(AssertionError):generate(p,4)
if __name__=='__main__':unittest.main()