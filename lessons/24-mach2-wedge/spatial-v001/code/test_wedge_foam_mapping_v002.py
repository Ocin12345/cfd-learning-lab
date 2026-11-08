import unittest
from pathlib import Path
import numpy as np
from wedge_foam_mapping_v002 import native_mapping
from analyse_nozzle_openfoam_v002 import field,FIELDS
R=Path('E:/CFD-Learning-Lab')
class Mapping(unittest.TestCase):
 def setUp(self):
  with np.load(R/'native/runs/mach2-wedge-fluent-v001/geometry.npz') as g:self.xy=g['centres'].copy()
  self.C=field(R/'native/runs/mach2-wedge-openfoam-v001/0.00025/C',22400,*FIELDS['C'])
 def test_actual_baseline_bijection(self):
  order,m=native_mapping(self.C,self.xy,.01);self.assertEqual(len(np.unique(order)),22400);self.assertLess(m['maximum_euclidean_error_m'],2e-10)
 def test_wrong_extrusion_rejected(self):
  a=self.C.copy();a[0,2]=.02
  with self.assertRaises(ValueError):native_mapping(a,self.xy,.01)
 def test_duplicate_native_cell_rejected(self):
  a=self.C.copy();a[-1]=a[0]
  with self.assertRaises(ValueError):native_mapping(a,self.xy,.01)
if __name__=='__main__':unittest.main()
