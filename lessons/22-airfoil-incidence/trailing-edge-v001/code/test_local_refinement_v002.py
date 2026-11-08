import unittest
import numpy as np
from local_refine_triangles_v002 import refine_mesh

def boundary_of(tri):
 counts={}
 for a,b,c in tri:
  for u,v in [(a,b),(b,c),(c,a)]:e=tuple(sorted((u,v)));counts[e]=counts.get(e,0)+1
 return {e:4 for e,n in counts.items() if n==1}
class LocalRefinement(unittest.TestCase):
 def verify(self,xy,tri,mark,expected):
  b=boundary_of(tri);q,t,bb,p=refine_mesh(xy,tri,b,np.array(mark,dtype=bool));self.assertEqual(len(t),expected);self.assertTrue(np.array_equal(q[:len(xy)],xy));self.assertAlmostEqual(p['area_before'],p['area_after'],places=12);self.assertEqual(set(boundary_of(t)),set(bb));self.assertTrue(p['conforming']);return q,t,bb,p
 def test_one_red_and_one_green(self):
  xy=np.array([[0.,0.],[1,0],[1,1],[0,1]]);tri=np.array([[0,1,2],[0,2,3]]);q,t,b,p=self.verify(xy,tri,[True,False],6);self.assertEqual(len(q),7)
 def test_two_edge_green_closure(self):
  xy=np.array([[0.,0.],[1,0],[.5,.8],[.5,-.8],[1.5,.8]]);tri=np.array([[0,1,2],[1,0,3],[2,1,4]]);q,t,b,p=self.verify(xy,tri,[False,True,True],11);self.assertEqual(sum(np.array(p['parent_ids'])==0),3)
 def test_empty_mark_keeps_mesh(self):
  xy=np.array([[0.,0.],[1,0],[0,1]]);tri=np.array([[0,1,2]]);q,t,b,p=self.verify(xy,tri,[False],1);self.assertTrue(np.array_equal(q,xy));self.assertTrue(np.array_equal(t,tri))
 def test_degenerate_triangle_rejected(self):
  xy=np.array([[0.,0.],[1,0],[2,0]]);tri=np.array([[0,1,2]])
  with self.assertRaises(ValueError):refine_mesh(xy,tri,boundary_of(tri),np.array([True]))
if __name__=='__main__':unittest.main()