from pathlib import Path
import unittest
import numpy as np
from triangle_quality_v002 import cell_quality
class NativeQualityPredictor(unittest.TestCase):
 def test_equilateral(self):
  xy=np.array([[0.,0.],[1,0],[.5,np.sqrt(3)/2]]);q=cell_quality(xy,np.array([[0,1,2]]));self.assertAlmostEqual(float(q[0]),1,places=12)
 def test_known_skinny_right_angle(self):
  xy=np.array([[0.,0.],[1,0],[0,.1]]);q=cell_quality(xy,np.array([[0,1,2]]));self.assertAlmostEqual(float(q[0]),float(.008/(.0101*3*np.sqrt(3))),places=12)
 def test_actual_rejected_native_mesh(self):
  p=Path('E:/CFD-Learning-Lab/native/meshes/showcase-airfoil-trailing-edge-v001/level1-v001/mesh-arrays.npz');d=np.load(p);q=cell_quality(d['xy'],d['triangles']);self.assertAlmostEqual(float(q.min()),.115009,delta=1e-6);self.assertLess(float(q.min()),.2)
if __name__=='__main__':unittest.main()