"""Closed, source-bound audit of all eight actual diagnostic states."""
import wing_import_paths_v001
from wing_repair_common_v001 import *
from wing_repair_native_v003 import raw_residual,same_native,state_keys
from wing_retained_residual_helpers_v001 import expected_probe_settings
from wing_residual_window_helpers_v001 import declared_iterations,verify_runtime
from wing_residual_window_prelaunch_v001 import verify
from wing_repair_helpers_v001 import same_residual_prefix,force_window
from wing_fluent_v001 import force_rows
from audit_wing_continuation_v001 import native
from plot_finite_wing_native_v002 import read_fluent,read_fluent_geometry,indexed,TetGeometry
from sloshing_fixed_settings_v001 import compare_fixed
from wing_amg_transcript_v001 import parse_blocks
import h5py,numpy as np,re
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
P=N/'fluent-retained-residual-probe-v002';O=N/'fluent-retained-residual-window-v001'

def audit_window():
 verify();assert json.loads((O/'status.json').read_text())['stage']=='closed';ends=declared_iterations(1502,1510)
 for x,y in [('wing-final.cas.h5','source.cas.h5'),('wing-final.dat.h5','source.dat.h5'),('final-settings.json','expected-source-settings.json')]:assert sha(P/x)==sha(O/y)
 compare_fixed(json.loads((O/'loaded-source-settings.json').read_text()),{k:v for k,v in json.loads((P/'final-settings.json').read_text()).items() if k in state_keys()})
 for name in ['candidate-settings.json','final-settings.json']:compare_fixed(json.loads((O/name).read_text()),expected_probe_settings())
 same_native(O/'loaded-restart.cas.h5',O/'loaded-restart.dat.h5',O/'source.cas.h5',O/'source.dat.h5')
 pi,pr=raw_residual(O/'source.dat.h5');assert pi[-1]==1502;fi,fr=raw_residual(O/'wing-final.dat.h5');assert fi[-1]==1510;same_residual_prefix(pi,pr,fi,fr,1502)
 q=json.loads((O/'geometry-quality.json').read_text());assert q['minimum_orthogonal_quality']>=.10 and q['minimum_cell_volume_m3']>0
 parent=force_rows((P/'forces-history.out').read_text())[-1];forces=force_window((O/'forces-history.out').read_text(),1503,1510,parent);blocks=parse_blocks((O/'solver-transcript.txt').read_text(),1503,1510)
 original=json.loads((N/'fluent-baseline-v001/qualification-audit-v001.json').read_text());g0=read_fluent_geometry(O/'source.cas.h5');n=len(g0['tets']);geom=TetGeometry(M/'geometry.npz');_,wing_ids=geom.map_native(read_fluent(O/'source.cas.h5',O/'source.dat.h5'));wing=np.zeros(n,bool);wing[g0['owner'][wing_ids]]=True;outer=np.zeros(n,bool);bf=np.flatnonzero(g0['neighbour']<0);outer[g0['owner'][np.setdiff1d(bf,wing_ids)]]=True;assert not np.any(wing&outer);wing_max_x=float(geom.wing_vertices[...,0].max())
 with np.load(W/'retained-probe-analysis-v002/mapped-retained-residuals.npz') as a:initial_res=a['SV_P_R'];initial_top=set(np.argsort(np.abs(initial_res))[::-1][:int(np.ceil(.01*n))].tolist())
 steps=[];paths=[O/x for x in ['source.cas.h5','source.dat.h5','expected-source-settings.json','loaded-restart.cas.h5','loaded-restart.dat.h5','wing-final.cas.h5','wing-final.dat.h5','solver-transcript.txt','forces-history.out','loaded-source-settings.json','candidate-settings.json','final-settings.json','initial-field-check.json','geometry-quality.json','solver-version.json','declared-test.json','status.json']]
 for j,end in enumerate(ends):
  case=O/f'iter-{end}.cas.h5';data=O/f'iter-{end}.dat.h5';capture=O/'runtime'/str(end);r,_,_,g=native(case,data);assert r['iteration']==end and r['cells']==224972 and r['relative_mass_imbalance']<=1e-4 and max(r['final_residuals'].values())<.1
  for k in ['points','tets','faces','owner','neighbour']:assert np.array_equal(g[k],g0[k]),(end,k)
  it,raw=raw_residual(data);same_residual_prefix(pi,pr,it,raw,1502);assert it[-1]==end and np.array_equal(it[-j-1:],np.arange(1503,end+1))
  values=np.array([raw[k][-1,0] for k in EQ]);scaled=np.array([raw[k][-1,0]/raw[k][-1,1] for k in EQ]);assert np.allclose(blocks[j]['initial'],values,rtol=6e-5,atol=1e-12) and np.allclose(blocks[j]['outer_printed'],scaled,rtol=6e-5,atol=1e-14)
  for k in EQ:assert np.isclose(r['final_residuals'][k],original['new_residual_history'][end-1501][k],rtol=1e-8,atol=1e-12),(end,k,'baseline replay')
  with h5py.File(data) as h:
   c=h['results/1/phase-1/cells'];u=np.column_stack([indexed(c[k],n) for k in ['SV_U','SV_V','SV_W']]);pressure=indexed(c['SV_P'],n)
  centres=g['points'][g['tets']].mean(axis=1)
  with np.load(capture/'runtime-residual-arrays.npz') as a:proof=verify_runtime(centres,u,pressure,a,{k:float(raw[k][-1,0]) for k in EQ},r['cell_mass_imbalance_abs_sum_kg_s'])
  guard=json.loads((capture/'native-runtime-guard.json').read_text());assert guard['iteration']==end and guard['capture_pass'] and guard['stage']=='assembled_initial' and guard['native']==r
  for p,h in guard['files_sha256'].items():assert sha(p)==h,p
  assert np.isclose(blocks[j]['terminal'][0],r['cell_mass_imbalance_abs_sum_kg_s'],rtol=6e-5,atol=1e-12)
  for label,k in [('drag','Cd'),('lift','Cl')]:
   assert np.isclose(forces[j][k],original['new_force_history'][end-1501][k],rtol=0,atol=1e-10)
   text=(O/'force-reports'/f'{label}-{end}.txt').read_text();found=re.findall(r'^Net\s+((?:'+NUM+r'\s+){5}'+NUM+r')\s*$',text,re.M);assert len(found)==1;v=list(map(float,found[0].split()));assert np.isclose(v[-1],forces[j][k],rtol=1e-6,atol=1e-10) and np.isclose(v[2]/2,v[-1],rtol=1e-6,atol=1e-10)
  text=(O/'force-reports'/f'mass-{end}.txt').read_text();found=re.findall(r'^\s*Net\s+('+NUM+r')\s*$',text,re.M);assert len(found)==1 and np.isclose(-float(found[0]),r['net_outward_mass_kg_s'],rtol=1e-4,atol=1e-12)
  res=proof['residuals']['SV_P_R'];mag=np.abs(res);order=np.argsort(mag)[::-1];sel=order[:int(np.ceil(.01*n))];total=float(mag.sum());top=set(sel.tolist());top_share=float(mag[sel].sum()/total);wing_share=float(mag[wing].sum()/total);minx=float(centres[sel,0].min())
  steps.append(dict(iteration=end,stage=proof['stage'],native=r,forces=forces[j],capture_pass=True,top1_cells=len(sel),top1_absolute_residual_fraction=top_share,top1_min_x_m=minx,top1_max_x_m=float(centres[sel,0].max()),top1_wing_adjacent_count=int(wing[sel].sum()),wing_adjacent_residual_fraction=wing_share,outer_boundary_adjacent_residual_fraction=float(mag[outer].sum()/total),downstream_of_wing_residual_fraction=float(mag[centres[:,0]>wing_max_x].sum()/total),top1_jaccard_with1502=len(top&initial_top)/len(top|initial_top),descriptive_persistence_pass=top_share>=.90 and minx>wing_max_x and wing_share<=.001,original_residual_pass=max(r['final_residuals'].values())<=1e-5,original_cell_mass_not_worse=r['cell_mass_imbalance_abs_sum_kg_s']<=6.536745231970203e-7,maximum_centroid_distance_m=proof['maximum_centroid_distance_m']))
  paths += [case,data,capture/'runtime-residual-arrays.npz',capture/'runtime-residual-catalogue.json',capture/'native-runtime-guard.json']+[O/'force-reports'/f'{label}-{end}.txt' for label in ['mass','drag','lift']]
 result=dict(diagnostic_capture_pass=True,additional_iterations=8,parent_iteration=1502,endpoint_iteration=1510,total_actual_steps_from_original1500=10,retained_map_iterations=[1502,*ends],first1501_map_unavailable_preserved=True,original_global_residual_limit=1e-5,original_full_qualification_not_claimed=True,original_qualification_pass=False,force100_window_not_available=True,steps=steps,wing_max_x_m=wing_max_x,descriptive_persistence_all_eight=all(s['descriptive_persistence_pass'] for s in steps),nine_map_conclusion_requires_parent1502_audit=True,no_new_algorithm_or_openfoam_solve=True,source_sha256={str(p):sha(p) for p in paths},physical_validation=False,full_cycle_or_time_independence=False)
 verify();dump(O/'retained-window-audit-v001.json',result);return result
if __name__=='__main__':print(json.dumps(audit_window(),indent=2))
