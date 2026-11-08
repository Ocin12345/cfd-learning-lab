import unittest
import numpy as np
import jet_energy_ledger_v001 as e
class Energy(unittest.TestCase):
 def api(self,name):
  self.assertTrue(hasattr(e,name),name+' not implemented');return getattr(e,name)
 def test_boundary_datum_shift_sign(self):
  self.assertEqual(self.api('fluent_heat300')({'inlet':20.,'outlet':-40.},{'inlet':-.01,'outlet':.01},1000.),{'inlet':10.,'outlet':-30.})
 def test_unequal_patch_keys_rejected(self):
  with self.assertRaises(ValueError):self.api('fluent_heat300')({'a':1.},{'b':0.},0.)
 def test_discrete_conservation_cumulative_and_signs(self):
  r=self.api('ledger')([1,2],[.5,1.],[0.,1.,2.],{'hot':[4.,6.],'out':[-2.,-4.]},.5)
  self.assertTrue(r['passed']);np.testing.assert_allclose(r['cumulative_net_input_J_per_m'],[1.,2.]);np.testing.assert_allclose(r['cumulative_patch_absolute_J_per_m'],[3.,8.])
 def test_cancellation_does_not_mask_local_misses(self):
  r=self.api('ledger')([1,2],[.5,1.],[0.,1.,2.],{'q':[4.,0.]},.5)
  self.assertFalse(r['passed']);self.assertEqual(r['step_relative_error'],[.5,1.]);self.assertEqual(r['cumulative_relative_error'][-1],0.)
 def test_right_endpoint_not_trapezoid(self):
  r=self.api('ledger')([1,2],[.1,.2],[0.,.1,.4],{'q':[1.,3.]},.1)
  self.assertTrue(r['passed']);self.assertAlmostEqual(r['cumulative_net_input_J_per_m'][-1],.4)
 def test_zero_energy_exact(self):
  r=self.api('ledger')([1],[.1],[0.,0.],{'q':[0.]},.1);self.assertTrue(r['passed']);self.assertEqual(r['step_relative_error'],[0.])
 def test_missing_duplicate_and_wrong_clock_rejected(self):
  f=self.api('ledger')
  for steps,times in [([1,3],[.1,.3]),([1,1],[.1,.1]),([1,2],[.1,.3])]:
   with self.assertRaises(ValueError):f(steps,times,[0.,1.,2.],{'q':[10.,10.]},.1)
 def test_nonfinite_or_wrong_shape_rejected(self):
  f=self.api('ledger')
  for energy,heat in [([0.,float('nan')],[1.]),([0.,1.],[float('inf')]),([0.],[1.])]:
   with self.assertRaises(ValueError):f([1],[.1],energy,{'q':heat},.1)
 def test_original_limit_unchanged(self):
  f=self.api('ledger');self.assertFalse(f([1],[1.],[0.,1.],{'q':[1.01]},1.)['passed']);self.assertTrue(f([1],[1.],[0.,1.],{'q':[1.004]},1.)['passed'])
 def test_outlet_upwind_and_depth_normalization(self):
  f=self.api('foam_patch_heat300')
  # No diffusion on outflow, -rho*cp*phi*(T-300)/depth=-100 W/m
  self.assertAlmostEqual(f('right',[310.],[.001],dx=.01,depth=.01,rho=1.,cp=100.,k=.1),-100.)
 def test_plate_diffusion_and_inlet_backflow(self):
  f=self.api('foam_patch_heat300')
  self.assertAlmostEqual(f('plate',[310.],[0.],dx=.01,depth=.01,rho=1.,cp=100.,k=.1),2.)
  self.assertAlmostEqual(f('right',[310.],[-.001],dx=.01,depth=.01,rho=1.,cp=100.,k=.1),-2.)
 def test_fluent_cap_initial_echo_not_terminal_iteration(self):
  f=self.api('fluent_steps')
  t=' 1750 0.0e+00 0:00:00 30\n 1751 1.0e-09 0:00:00 29\nFlow time = 0.005s, time step = 1\n'
  r=f(t,1,1,.005,30);self.assertEqual(r[0]['inner_iterations'],1);self.assertEqual(r[0]['residual'],1e-9)
 def test_fluent_no_new_iteration_or_wrong_step_rejected(self):
  f=self.api('fluent_steps')
  for t in [' 1 1e-9 0:00:00 30\nFlow time = 0.005s, time step = 1',' 1 1e-9 0:00:00 29\nFlow time = 0.005s, time step = 2']:
   with self.assertRaises(ValueError):f(t,1,1,.005,30)
 def test_continuation_global_steps_and_missing_report_rejected(self):
  f=self.api('fluent_steps');t=' 501 1e-9 0:00:00 29\nFlow time = 0.105s, time step = 21'
  self.assertEqual(f(t,21,21,.005,30)[0]['step'],21)
  with self.assertRaises(ValueError):f(t,21,22,.005,30)
if __name__=='__main__':unittest.main()
