import unittest
from wedge_endpoint_v001 import accepted
class Endpoints(unittest.TestCase):
 def values(self):return dict(relative_field_change_l2=dict(p=.001,mach=.001),outlet_minimum_cell_mach=1.64,maximum_ideal_gas_relative_error=1e-6,analytical_angle_check=True,analytical_pressure_check=True,shock_fit_columns=192)
 def residuals(self):return dict(continuity=1e-7,**{'x-velocity':1e-7,'y-velocity':1e-7},energy=1e-8)
 def test_steady_native_endpoint_admitted(self):self.assertTrue(all(accepted(self.residuals(),self.values(),1e-7).values()))
 def test_residual_pass_cannot_hide_field_drift(self):
  m=self.values();m['relative_field_change_l2']['p']=.003;self.assertFalse(accepted(self.residuals(),m,1e-7)['field_stability'])
 def test_native_mass_and_outlet_limits(self):
  m=self.values();m['outlet_minimum_cell_mach']=.99;a=accepted(self.residuals(),m,.002);self.assertFalse(a['mass_balance']);self.assertFalse(a['supersonic_outlet'])
 def test_energy_or_coverage_miss_retained(self):
  r=self.residuals();r['energy']=2e-6;m=self.values();m['shock_fit_columns']=159;a=accepted(r,m,1e-7);self.assertFalse(a['residuals']);self.assertFalse(a['column_coverage'])
if __name__=='__main__':unittest.main()
