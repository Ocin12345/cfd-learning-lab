import copy,ast
from sloshing_fixed_settings_v001 import compare_fixed

def _remove(d,path):
 node=d
 for key in path[:-1]:
  if key not in node:return
  node=node[key]
 node.pop(path[-1],None)

def fixed_projection(settings):
 d=copy.deepcopy(settings);coupling=d['methods']['p_v_coupling'];d['methods']['p_v_coupling']={k:coupling[k] for k in ['flux_type','flux_auto_select']}
 for path in [('methods','pseudo_time_method','formulation'),('controls','p_v_controls'),('controls','under_relaxation','pressure'),('controls','under_relaxation','mom'),('controls','advanced','multi_grid','mg_controls'),('controls','advanced','multi_grid','amg_controls','coupled_parameters'),('controls','advanced','multi_grid','amg_controls','options','verbosity'),('controls','advanced','expert','save_cell_residuals')]:_remove(d,path)
 return d

def validate_state(source,candidate):
 try:compare_fixed(fixed_projection(source),fixed_projection(candidate))
 except Exception as e:raise ValueError('Unexpected fixed physics/spatial/control change: '+str(e)) from e
 p=candidate['methods']['p_v_coupling'];allowed={'flow_scheme','skewness_correction_itr_count','neighbor_correction_itr_count','skewness_neighbor_coupling','coupled_form','solve_n_phase','flux_type','flux_auto_select'}
 if not set(p).issubset(allowed) or p.get('flow_scheme')!='SIMPLEC' or p.get('skewness_correction_itr_count')!=1:raise ValueError('Undeclared coupling method/keys/corrector')
 if p.get('solve_n_phase',False):raise ValueError('Unexpected multiphase solve')
 formulation=candidate['methods']['pseudo_time_method']['formulation']
 if formulation!={'segregated_solver':'off'}:raise ValueError('Pseudo time must remain off')
 controls=candidate['controls'];u=controls['under_relaxation']
 if u.get('pressure')!=.7 or u.get('mom')!=.7:raise ValueError('Wrong prescribed implicit relaxation')
 mg=controls['advanced']['multi_grid']['mg_controls']
 if 'pressure' not in mg or mg['pressure']['termination_criteria']!=.001 or not set(mg).issubset({'pressure','momentum','x-velocity','y-velocity','z-velocity'}):raise ValueError('Wrong active pressure AMG control')
 if not controls['advanced']['expert']['save_cell_residuals'] or controls['advanced']['multi_grid']['amg_controls']['options']['verbosity']!=1:raise ValueError('Diagnostic outputs not active')
 return True

def _flatten(value,path=()):
 if isinstance(value,dict) and value:
  result={}
  for k,v in value.items():result.update(_flatten(v,path+(k,)))
  return result
 return {path:value}

def diff_inventory(source,candidate):
 a=_flatten(source);b=_flatten(candidate);result=[]
 for path in sorted(set(a)|set(b)):
  if path not in a or path not in b or a[path]!=b[path] or type(a[path]) is not type(b[path]):result.append(dict(path='/'.join(path),before=a.get(path,{'__missing__':True}),after=b.get(path,{'__missing__':True})))
 return result

def assert_zero_advance_source(text):
 forbidden={'iterate','dual_time_iterate','initialize','hybrid_initialize','standard_initialize','patch','renormalize','reset','reset_residuals','execute_tui','remesh','smooth','improve'}
 for n in ast.walk(ast.parse(text)):
  if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in forbidden:raise ValueError('Forbidden advance/reset/mesh call: '+n.func.attr)
 return True
