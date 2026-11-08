import wing_import_paths_v001
import unittest,copy,ast
from wing_simplec_capability_helpers_v001 import validate_state,diff_inventory,assert_zero_advance_source
class Capability(unittest.TestCase):
 def pair(self):
  src=dict(general={'steady':True},material={'rho':1.,'mu':.0025},references={'area':4.},methods={'p_v_coupling':{'flow_scheme':'Coupled','flux_type':'Rhie-Chow: momentum based','flux_auto_select':True},'pseudo_time_method':{'formulation':{'coupled_solver':'off'}},'spatial_discretization':{'momentum':'second-order-upwind'}},controls={'p_v_controls':{'flow_courant_number':20.},'under_relaxation':{'body-force':1.,'density':1},'advanced':{'multi_grid':{'mg_controls':{'pressure-coupled':{'termination_criteria':.1}},'amg_controls':{'options':{'verbosity':0},'scalar_parameters':{'max_cycle':30},'coupled_parameters':{'max_cycle':10}}},'expert':{'save_cell_residuals':False}}},residual={'scale':True,'normalize':False})
  dst=copy.deepcopy(src);dst['methods']['p_v_coupling'].update(flow_scheme='SIMPLEC',skewness_correction_itr_count=1);dst['methods']['pseudo_time_method']['formulation']={'segregated_solver':'off'};dst['controls'].pop('p_v_controls');dst['controls']['under_relaxation'].update(pressure=.7,mom=.7);dst['controls']['advanced']['multi_grid']['mg_controls']={'pressure':{'termination_criteria':.001},'momentum':{'termination_criteria':.1}};dst['controls']['advanced']['multi_grid']['amg_controls']['options']['verbosity']=1;dst['controls']['advanced']['expert']['save_cell_residuals']=True
  return src,dst
 def test_declared_transition_only(self):
  src,dst=self.pair();self.assertTrue(validate_state(src,dst));self.assertTrue(diff_inventory(src,dst))
 def test_fixed_physics_change_refused(self):
  src,dst=self.pair();dst['material']['rho']=1.01
  with self.assertRaises(ValueError):validate_state(src,dst)
 def test_spatial_scheme_change_refused(self):
  src,dst=self.pair();dst['methods']['spatial_discretization']['momentum']='first-order-upwind'
  with self.assertRaises(ValueError):validate_state(src,dst)
 def test_other_scalar_amg_change_refused(self):
  src,dst=self.pair();dst['controls']['advanced']['multi_grid']['amg_controls']['scalar_parameters']['max_cycle']=99
  with self.assertRaises(ValueError):validate_state(src,dst)
 def test_wrong_declared_controls_refused(self):
  src,dst=self.pair();dst['controls']['under_relaxation']['pressure']=.3
  with self.assertRaises(ValueError):validate_state(src,dst)
 def test_active_pseudo_time_refused(self):
  src,dst=self.pair();dst['methods']['pseudo_time_method']['formulation']['segregated_solver']='global-time-step'
  with self.assertRaises(ValueError):validate_state(src,dst)
 def test_exact_leaf_inventory(self):
  r=diff_inventory({'a':{'x':1},'b':{}},{'a':{'y':2},'c':False});self.assertEqual({v['path'] for v in r},{'a/x','a/y','b','c'})
 def test_zero_advance_source_contract(self):
  self.assertTrue(assert_zero_advance_source('session.settings.file.write_case_data(file_name="x")'))
  for code in ['st.solution.run_calculation.iterate(iter_count=1)','session.execute_tui("/solve/iterate 1")','st.solution.initialization.hybrid_initialize()']:
   with self.assertRaises(ValueError):assert_zero_advance_source(code)
if __name__=='__main__':unittest.main()
