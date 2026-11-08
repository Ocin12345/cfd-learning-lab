"""Closed finite SIMPLEC qualification audit; no Coupled-stage assumption."""
import wing_import_paths_v001
from wing_repair_common_v001 import *
from wing_repair_native_v003 import raw_residual,same_native,state_keys,NUM
from wing_repair_helpers_v001 import same_residual_prefix,residual_window,force_window
from wing_simplec_qualification_helpers_v001 import checkpoints,scale_guard
from wing_simplec_runtime_v001 import verify_simplec_runtime
from wing_simplec_capability_helpers_v002 import validate_state
from wing_simplec_qualification_prelaunch_v001 import verify
from wing_fluent_v001 import force_stability
from audit_wing_continuation_v001 import native
from plot_finite_wing_native_v002 import read_fluent_geometry,indexed
from sloshing_fixed_settings_v001 import compare_fixed
import h5py,numpy as np,re
O=N/'fluent-simplec-qualification-v001'

def audit():
 verify();assert json.loads((O/'status.json').read_text())['stage']=='closed'
 for x,y in [('wing-final.cas.h5','source.cas.h5'),('wing-final.dat.h5','source.dat.h5'),('final-settings.json','expected-source-settings.json')]:assert sha(F/x)==sha(O/y)
 source={k:v for k,v in json.loads((F/'final-settings.json').read_text()).items() if k in state_keys()};compare_fixed(json.loads((O/'loaded-source-settings.json').read_text()),source)
 expected=json.loads((W/'accepted-native-SIMPLEC-settings-v001.json').read_text());validate_state(source,expected)
 for name in ['candidate-settings.json','final-settings.json']:compare_fixed(json.loads((O/name).read_text()),expected)
 for name in ['loaded-targets.json','candidate-targets.json']:assert json.loads((O/name).read_text())==dict.fromkeys(EQ,1e-5)
 for name in ['loaded-restart','configured-start']:same_native(O/(name+'.cas.h5'),O/(name+'.dat.h5'),O/'source.cas.h5',O/'source.dat.h5')
 for label in ['loaded','configured']:
  d=json.loads((O/(label+'-live-iteration.json')).read_text());assert d['query_value']==d['expected']==1500
 pi,pr=raw_residual(O/'source.dat.h5');assert np.array_equal(pi,np.arange(1,1501));scale=float(pr['continuity'][-1,1]);assert scale==2.2344914086503613
 for stem in ['loaded-restart','configured-start']:
  ci,cr=raw_residual(O/(stem+'.dat.h5'));same_residual_prefix(pi,pr,ci,cr,1500);assert ci[-1]==1500 and len(ci)==1500
 fi,fr=raw_residual(O/'wing-final.dat.h5');same_residual_prefix(pi,pr,fi,fr,1500);scale_guard(fi,fr['continuity'],scale,1600);scaled={k:fr[k][:,0]/fr[k][:,1] for k in EQ};rows=residual_window(fi,scaled,1501,1600)
 parent=json.loads((W/'native-evidence-inventory-v001.json').read_text())['Fluent']['forces'];forces=force_window((O/'forces-history.out').read_text(),1501,1600,parent);fs=force_stability(forces,1600)
 quality=json.loads((O/'geometry-quality.json').read_text());assert quality['minimum_orthogonal_quality']>=.10 and quality['minimum_cell_volume_m3']>0
 g0=read_fluent_geometry(O/'source.cas.h5');n=len(g0['tets']);centres=g0['points'][g0['tets']].mean(axis=1);states=[]
 paths=[O/x for x in ['source.cas.h5','source.dat.h5','expected-source-settings.json','loaded-restart.cas.h5','loaded-restart.dat.h5','configured-start.cas.h5','configured-start.dat.h5','wing-final.cas.h5','wing-final.dat.h5','solver-transcript.txt','forces-history.out','loaded-source-settings.json','candidate-settings.json','final-settings.json','loaded-field-proof.json','configured-field-proof.json','loaded-targets.json','candidate-targets.json','geometry-quality.json','solver-version.json','declared-test.json','launch-admission.json','status.json','loaded-live-iteration.json','configured-live-iteration.json','native-state-differences.json','active-AMG-names.json','active-under-relaxation-names.json','allowed-flow-schemes.json','after-switch-defaults-settings.json']]
 for end in checkpoints():
  case=O/f'iter-{end}.cas.h5';data=O/f'iter-{end}.dat.h5';r,it,res,g=native(case,data);assert r['iteration']==end and r['cells']==224972 and r['relative_mass_imbalance']<=1e-4 and max(r['final_residuals'].values())<.1
  for k in ['points','tets','faces','owner','neighbour']:assert np.array_equal(g[k],g0[k]),(end,k)
  ci,cr=raw_residual(data);same_residual_prefix(pi,pr,ci,cr,1500);scale_guard(ci,cr['continuity'],scale,end);assert all(np.isfinite(res[k]).all() for k in EQ)
  compare_fixed(json.loads((O/f'settings-{end}.json').read_text()),expected);assert json.loads((O/f'native-guard-{end}.json').read_text())==r
  live=json.loads((O/f'actual{end}-live-iteration.json').read_text());assert live['query_value']==live['expected']==end
  for label,k in [('drag','Cd'),('lift','Cl')]:
   text=(O/'force-reports'/f'{label}-{end}.txt').read_text();a=re.findall(r'^Net\s+((?:'+NUM+r'\s+){5}'+NUM+r')\s*$',text,re.M);assert len(a)==1;v=list(map(float,a[0].split()));assert np.isclose(v[-1],forces[end-1501][k],rtol=1e-6,atol=1e-10) and np.isclose(v[2]/2,v[-1],rtol=1e-6,atol=1e-10)
  a=re.findall(r'^\s*Net\s+('+NUM+r')\s*$',(O/'force-reports'/f'mass-{end}.txt').read_text(),re.M);assert len(a)==1 and np.isclose(-float(a[0]),r['net_outward_mass_kg_s'],rtol=1e-4,atol=1e-12)
  with h5py.File(data) as h:
   c=h['results/1/phase-1/cells'];u=np.column_stack([indexed(c[k],n) for k in ['SV_U','SV_V','SV_W']]);p=indexed(c['SV_P'],n)
  capture=O/'runtime'/str(end)
  with np.load(capture/'runtime-residual-arrays.npz') as a:proof=verify_simplec_runtime(centres,u,p,a,{k:float(cr[k][-1,0]) for k in EQ},r['cell_mass_imbalance_abs_sum_kg_s'])
  states.append(dict(iteration=end,native=r,runtime=proof,all_native_settings_match=True,force_reports_match_2N=True))
  paths += [case,data,O/f'settings-{end}.json',O/f'native-guard-{end}.json',O/f'actual{end}-live-iteration.json',capture/'runtime-residual-arrays.npz',capture/'runtime-residual-catalogue.json']+[O/'force-reports'/f'{s}-{end}.txt' for s in ['mass','drag','lift']]
 same_native(O/'wing-final.cas.h5',O/'wing-final.dat.h5',O/'iter-1600.cas.h5',O/'iter-1600.dat.h5');r,it,res,g=native(O/'wing-final.cas.h5',O/'wing-final.dat.h5');assert r['iteration']==1600
 checks=dict(structure_and_admitted_native_settings=True,initial_fields_topology_history_exact=True,continuity_scale_exact=True,full100_actual_residual_and_force_coverage=True,endpoint_residuals=max(r['final_residuals'].values())<=1e-5,all_final20_residuals=max(max(v[k] for k in EQ) for v in rows[-20:])<=1e-5,mass=r['relative_mass_imbalance']<=1e-4,force_stability=bool(fs['passed']),cell_mass_not_worse=r['cell_mass_imbalance_abs_sum_kg_s']<=6.536745231970203e-7)
 report=dict(structural_pass=True,numerical_pass=all(checks.values()),checks=checks,algorithm='SIMPLEC declared bundle',parent_iteration=1500,new_iterations=100,endpoint_iteration=1600,continuity_scale_exact=scale,native=r,forces=forces[-1],force_stability=fs,new_residual_history=rows,new_force_history=forces,states=states,momentum_denominator_history={k:fr[k][fi>1500,1].tolist() for k in EQ if k!='continuity'},native_field_and_mass_evidence_scope='five actual saved states, not all100',source_sha256={str(p):sha(p) for p in paths},physical_validation=False,mesh_independence=False,steady_iterations_not_physical_time=True,confirmation_not_admitted_until_actual_visual_inspection=True)
 verify();assert not (O/'qualification-audit-v001.json').exists();dump(O/'qualification-audit-v001.json',report);return report
if __name__=='__main__':print(json.dumps(audit(),indent=2))
