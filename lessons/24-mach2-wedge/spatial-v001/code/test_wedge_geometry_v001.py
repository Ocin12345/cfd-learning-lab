import hashlib,tempfile,unittest,uuid
from pathlib import Path
import numpy as np
from wedge_mesh_v001 import make_mesh,compare_refinement
ROOT=Path('E:/CFD-Learning-Lab');TEMP=ROOT/'temp'
def fresh():
 p=TEMP/('wedge-geometry-test-'+uuid.uuid4().hex);p.mkdir();return p
class Geometry(unittest.TestCase):
 def test_positive_area_and_same_domain(self):
  a,b=fresh(),fresh();make_mesh(a,10,6);make_mesh(b,20,12)
  r=compare_refinement(a/'geometry.npz',b/'geometry.npz');self.assertTrue(r['same_boundary_polygon']);self.assertTrue(r['old_nodes_preserved']);self.assertLess(r['maximum_parent_area_partition_error_m2'],1e-14)
 def test_corner_must_be_a_node(self):
  with self.assertRaises(ValueError):make_mesh(fresh(),11,6)
 def test_baseline_ascii_mesh_reproduced(self):
  p=fresh();make_mesh(p,200,112);old=ROOT/'native/runs/mach2-wedge-fluent-v001/wedge.msh'
  self.assertEqual(hashlib.sha256((p/'wedge.msh').read_bytes()).hexdigest(),hashlib.sha256(old.read_bytes()).hexdigest())
if __name__=='__main__':unittest.main()
