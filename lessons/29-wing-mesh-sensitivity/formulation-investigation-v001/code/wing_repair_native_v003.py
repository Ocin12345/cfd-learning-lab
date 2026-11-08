import h5py,numpy as np,json,re
from wing_repair_common_v001 import *
from wing_repair_helpers_v001 import same_residual_prefix,residual_window,force_window
from plot_finite_wing_native_v002 import read_fluent_geometry,indexed
from audit_wing_continuation_v001 import native
from wing_fluent_v001 import force_stability
from sloshing_fixed_settings_v001 import compare_fixed
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'

def state(st):
 return dict(general=st.setup.general.get_state(),viscous=st.setup.models.viscous.get_state(),material=st.setup.materials.fluid['air'].get_state(),cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),inlet=st.setup.boundary_conditions.velocity_inlet['inlet'].get_state(),farfield=st.setup.boundary_conditions.velocity_inlet['farfield'].get_state(),outlet=st.setup.boundary_conditions.pressure_outlet['outlet'].get_state(),walls=dict(wing=st.setup.boundary_conditions.wall['wing'].get_state()),references=st.setup.reference_values.get_state(),methods=st.solution.methods.get_state(),controls=st.solution.controls.get_state(),residual=st.solution.monitor.residual.get_state(),report_definitions=st.solution.report_definitions.get_state())
def expected_settings(criterion,courant=20.,relaxation=.5):
 if relaxation not in (.25,.5):raise ValueError("Undeclared coupled relaxation")
 a=json.loads((F/'final-settings.json').read_text());a={k:a[k] for k in ['general','viscous','material','cell_zones','inlet','farfield','outlet','walls','references','methods','controls','residual','report_definitions']}
 a['controls']['p_v_controls']['explicit_pressure_under_relaxation']=relaxation;a['controls']['p_v_controls']['explicit_momentum_under_relaxation']=relaxation;a['controls']['p_v_controls']['flow_courant_number']=courant;a['controls']['advanced']['multi_grid']['mg_controls']['pressure-coupled']['termination_criteria']=criterion;a['controls']['advanced']['multi_grid']['amg_controls']['options']['verbosity']=1
 return a
def raw_residual(data):
 with h5py.File(data) as h:
  it=None;raw={}
  for n in EQ:
   g=h['results/residuals/phase-1/'+n];i=g['iterations'][...].reshape(-1);v=g['data'][...]
   assert np.isfinite(v[:,:2]).all() and (v[:,1]>0).all()
   if it is not None:assert np.array_equal(i,it)
   it=i;raw[n]=v
 return it,raw

def same_native(case,data,source_case,source_data):
 a=read_fluent_geometry(case);b=read_fluent_geometry(source_case)
 for n in ['points','tets','faces','owner','neighbour']:assert np.array_equal(a[n],b[n]),'Native geometry/order changed:'+n
 n=len(a['tets']);nf=len(a['faces']);assert n==224972
 with h5py.File(data) as x,h5py.File(source_data) as y:
  for name,kind,num in [('SV_U','cells',n),('SV_V','cells',n),('SV_W','cells',n),('SV_P','cells',n),('SV_DENSITY','cells',n),('SV_MU_LAM','cells',n),('SV_FLUX','faces',nf)]:
   p='results/1/phase-1/'+kind+'/'+name;v=indexed(x[p],num);w=indexed(y[p],num);assert np.isfinite(v).all() and np.array_equal(v,w),'Changed restartfield:'+name
 return dict(passed=True,cells=n,geometry_order_exact=True,initial_UVW_p_rho_mu_phi_exact=True)

def audit(label):
 o=N/('fluent-'+label+'-v001');criterion=.1 if label=='baseline' else .001;courant=10. if label in ('low-cfl10','low-relax25') else 20.;relaxation=.25 if label=='low-relax25' else .5;status=json.loads((o/'status.json').read_text());assert status['stage']=='closed'
 assert sha(o/'source.cas.h5')==sha(F/'wing-final.cas.h5') and sha(o/'source.dat.h5')==sha(F/'wing-final.dat.h5')
 compare_fixed(json.loads((o/'loaded-source-settings.json').read_text()),{k:v for k,v in json.loads((F/'final-settings.json').read_text()).items() if k in state_keys()})
 compare_fixed(json.loads((o/'candidate-settings.json').read_text()),expected_settings(criterion,courant,relaxation));compare_fixed(json.loads((o/'final-settings.json').read_text()),expected_settings(criterion,courant,relaxation))
 init=same_native(o/'loaded-restart.cas.h5',o/'loaded-restart.dat.h5',o/'source.cas.h5',o/'source.dat.h5');pi,pr=raw_residual(o/'source.dat.h5');ni,nr=raw_residual(o/'wing-final.dat.h5');same_residual_prefix(pi,pr,ni,nr,1500)
 scaled={n:nr[n][:,0]/nr[n][:,1] for n in EQ};rows=residual_window(ni,scaled,1501,1600)
 parent=json.loads((W/'native-evidence-inventory-v001.json').read_text())['Fluent']['forces'];forces=force_window((o/'forces-history.out').read_text(),1501,1600,parent);fs=force_stability(forces,1600)
 # Actual final native mesh/fields/material/flux plus cell and global mass closure.
 result,it,res,g=native(o/'wing-final.cas.h5',o/'wing-final.dat.h5');assert result['iteration']==1600 and result['cells']==224972
 original=read_fluent_geometry(o/'source.cas.h5')
 for n in ['points','tets','faces','owner','neighbour']:assert np.array_equal(g[n],original[n]),n
 quality=json.loads((o/'geometry-quality.json').read_text());assert quality['minimum_orthogonal_quality']>=.10 and quality['minimum_cell_volume_m3']>0
 for label_force,k in [('drag','Cd'),('lift','Cl')]:
  text=(o/'force-reports'/(label_force+'-1600.txt')).read_text();found=re.findall(r'^Net\s+((?:'+NUM+r'\s+){5}'+NUM+r')\s*$',text,re.M);assert len(found)==1;v=list(map(float,found[0].split()));assert np.isclose(v[-1],forces[-1][k],rtol=1e-6,atol=1e-10) and np.isclose(v[2]/2,v[-1],rtol=1e-6,atol=1e-10)
 mass=(o/'force-reports/mass-1600.txt').read_text();net=re.findall(r'^\s*Net\s+('+NUM+r')\s*$',mass,re.M);assert len(net)==1 and np.isclose(-float(net[0]),result['net_outward_mass_kg_s'],rtol=1e-4,atol=1e-12)
 checks=dict(structure_and_fixed_settings=True,initial_and_geometry_exact=True,residual_endpoint=max(result['final_residuals'].values())<=1e-5,all_final20_residuals=max(max(r[n] for n in EQ) for r in rows[-20:])<=1e-5,mass=result['relative_mass_imbalance']<=1e-4,force_stability=fs['passed'],cell_mass_not_worse=result['cell_mass_imbalance_abs_sum_kg_s']<=6.536745231970203e-7)
 paths=[o/n for n in ['source.cas.h5','source.dat.h5','loaded-restart.cas.h5','loaded-restart.dat.h5','wing-final.cas.h5','wing-final.dat.h5','forces-history.out','solver-transcript.txt','loaded-source-settings.json','candidate-settings.json','final-settings.json','initial-field-check.json','geometry-quality.json','declared-test.json','solver-version.json']]+list((o/'force-reports').glob('*.txt'))+list(o.glob('iter-*.h5'))
 report=dict(label=label,numerical_pass=all(checks.values()),structural_pass=True,checks=checks,parent_iteration=1500,new_iterations=100,endpoint_iteration=1600,criterion=criterion,courant=courant,relaxation=relaxation,native=result,forces=forces[-1],force_stability=fs,new_residual_history=rows,new_force_history=forces,source_sha256={str(p):sha(p) for p in paths},physical_validation=False,mesh_independence=False,steady_iterations_not_physical_time=True)
 dump(o/'qualification-audit-v001.json',report);return report

def state_keys():return ['general','viscous','material','cell_zones','inlet','farfield','outlet','walls','references','methods','controls','residual','report_definitions']
