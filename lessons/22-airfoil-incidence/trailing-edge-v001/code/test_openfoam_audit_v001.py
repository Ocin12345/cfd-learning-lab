import unittest
from audit_refined_openfoam_v001 import outcome
class Contract(unittest.TestCase):
 def rows(self):return [dict(Cd=.42,Cl=.26) for _ in range(100)]
 def test_no_zero_angle_gate(self):
  a=outcome(dict(Ux=1e-8,Uy=1e-8,p=1e-8),dict(inlet=-1,outlet=1,airfoil=0),self.rows())
  self.assertTrue(all(a['checks'].values()))
 def test_worst_pressure_residual_miss_preserved(self):
  a=outcome(dict(Ux=1e-8,Uy=1e-8,p=1.1e-7),dict(inlet=-1,outlet=1,airfoil=0),self.rows())
  self.assertFalse(a['checks']['openfoam.residuals'])
 def test_incomplete_window_rejected(self):
  with self.assertRaises(ValueError):outcome(dict(Ux=1e-8,Uy=1e-8,p=1e-8),dict(inlet=-1,outlet=1,airfoil=0),self.rows()[:-1])
 def test_mass_and_force_misses_preserved(self):
  r=self.rows();r[-1]['Cl']+=2e-5
  a=outcome(dict(Ux=1e-8,Uy=1e-8,p=1e-8),dict(inlet=-1,outlet=1.00002,airfoil=0),r)
  self.assertFalse(a['checks']['openfoam.mass_balance']);self.assertFalse(a['checks']['openfoam.Cl_stability'])
if __name__=='__main__':unittest.main()
