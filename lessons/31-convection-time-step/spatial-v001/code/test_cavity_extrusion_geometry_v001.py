import unittest,numpy as np
from pathlib import Path
from cavity_extrusion_geometry_v001 import verify_extrusion
from analyse_nozzle_openfoam_v002 import field,FIELDS
from plot_nozzle_native_v004 import foam_list
class Extrusion(unittest.TestCase):
 def data(self,fine):
  r=Path('E:/CFD-Learning-Lab');c=r/('native/runs/showcase-convection-spatial-v001/openfoam-control-v001' if fine else 'native/runs/natural-convection-cavity-48-v001/openfoam-control');t='0' if fine else '600';n=96 if fine else 48
  C=field(c/t/'C',n*n,*FIELDS['C']);mesh=c/'0/polyMesh';mesh=mesh if mesh.exists() else c/'constant/polyMesh';_,text=foam_list(mesh/'points');points=np.fromstring(text.replace('(',' ').replace(')',' '),sep=' ').reshape(-1,3)
  x=(np.arange(n)+.5)*.05/n;xx,yy=np.meshgrid(x,x);reference=np.stack([xx,yy],axis=-1).reshape(-1,2);return C,points,reference
 def test_actual_coarse_and_fine(self):
  for flag in (False,True):
   c,p,g=self.data(flag);a=verify_extrusion(c,p,g);self.assertAlmostEqual(a['native_depth_m'],.005);self.assertEqual(a['expected_z_bounds_m'],[-.0025,.0025])
 def test_wrong_centroid_plane_rejected(self):
  c,p,g=self.data(True);c[:,2]+=.0025
  with self.assertRaises(ValueError):verify_extrusion(c,p,g)
 def test_wrong_depth_rejected(self):
  c,p,g=self.data(True);p[:,2]*=2
  with self.assertRaises(ValueError):verify_extrusion(c,p,g)
if __name__=='__main__':unittest.main()