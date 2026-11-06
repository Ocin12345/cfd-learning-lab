"""Read-only audit of101 actual Fluent solid-slab snapshots; no solver launch."""
from pathlib import Path
import argparse,csv,json,re,shutil,math
import h5py,numpy as np
from thermal_slab_fluent_common_v003 import (sha,write_json,source_files,verify_sources,geometry_check,physical_check,topology,native_time,temperatures,state_metrics,step_metrics,energy_residuals,TARGETS,WALLS,NUMBER)

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 run=a.run.resolve();out=a.out.resolve();assert run.is_dir() and out.drive.upper()=='E:' and not out.exists(),'Fresh E: analysis directory required'
 out.mkdir(parents=True);checks={};errors={};evidence={};before={}
 def remember(path):
  p=Path(path).resolve();assert p.is_file(),str(p);before.setdefault(p,sha(p));return p
 def load(name):return json.loads(remember(run/name).read_text(encoding='utf-8-sig'))
 def attempt(name,fn):
  try:v=fn();checks[name]=True;return v
  except Exception as e:checks[name]=False;errors[name]=type(e).__name__+': '+str(e);return None
 status=attempt('execution_status_read',lambda:load('status.json')) or {};summary=attempt('summary_read',lambda:load('summary.json')) or {}
 auth=attempt('authorization_read',lambda:load('authorization.json')) or {};source=attempt('source_provenance_read',lambda:load('code-provenance.json')) or {}
 targets=attempt('declared_targets_read',lambda:load('declared-targets.json')) or {};quality=attempt('native_quality_read',lambda:load('geometry-quality.json')) or {}
 mesh=Path(source.get('mesh',str(run.parent/'mesh/slab.msh'))).resolve();geometry_path=mesh.parent/'geometry.npz';setup_path=mesh.parent/'setup.json'
 def provenance():
  assert auth['action']=='run_thermal_slab_fluent_v003' and auth['root']=='root' and auth['approved'] is True and auth['allowed_mode']=='Run' and Path(auth['out']).resolve()==run
  assert source['source_sha256']==sha(remember(run/'executed-script.py'))==auth['source_sha256']
  assert source['mesh_sha256']==sha(remember(mesh))==auth['mesh_sha256']
  assert source['geometry_npz_sha256']==sha(remember(geometry_path)) and source['mesh_setup_sha256']==sha(remember(setup_path))
  assert source['helper_sha256']==auth['helper_sha256'];verify_sources(auth['helper_sha256'])
  for name,digest in auth['helper_sha256'].items():assert sha(remember(run/'sources'/name))==digest
  assert source['snapshot_before_solver_launch'] is True and source['available_ram_before_launch_bytes']>=4*1024**3 and source['free_disk_before_launch_bytes']>=10*1024**3
  assert targets==TARGETS
  return dict(source_sha256=source['source_sha256'],helper_sha256=source['helper_sha256'],mesh_sha256=source['mesh_sha256'],geometry_npz_sha256=source['geometry_npz_sha256'])
 evidence['provenance']=attempt('all_reviewed_sources_mesh_and_targets_bound',provenance)
 checks['native_completed100steps200s101states']=status.get('stage') in ('completed','completed_with_export_errors') and summary.get('actual_steps')==100 and summary.get('actual_time_s')==200 and summary.get('native_states')==101 and summary.get('native_final_saved') is True
 checks['no_export_errors']=summary.get('export_errors')==[]
 checks['native_positive_volume0p002_and_quality']=quality.get('minimum_orthogonal_quality',-1)>=.15 and quality.get('minimum_cell_volume_m3',0)>0 and abs(quality.get('total_cell_volume_m3',0)-.002)<1e-10 and quality.get('checked_before_energy_advance') is True
 evidence['native_quality']=quality
 for name in ('loaded-settings.json','initialized-settings.json','first-step-settings.json','final-settings.json'):
  attempt('actual_solid_material_thermal_Euler_settings_'+name,lambda name=name:physical_check(load(name)))
 def gate():
  value=load('first-step-gate.json');assert value['passed'] is True and value['actual_step']==1 and value['actual_time_s']==2
  for key in ('fluid_zones_absent','constant_material_verified','native_temperature_time_verified','all240_native_boundary_heat_powers_verified','report_HDF_integrated_W_inward_sign_verified','global_Euler_energy_balance_verified'):assert value[key] is True
  return value
 evidence['runtime_first_step_gate']=attempt('first_step_native_proof_gate_passed_before_continuation',gate)
 geometry=attempt('shared_geometry_read',lambda:dict(np.load(remember(geometry_path),allow_pickle=False)))
 case=attempt('shared_native_case_read',lambda:remember(run/'slab.cas.h5'));final_case=attempt('final_native_case_read',lambda:remember(run/'slab-final.cas.h5'));final_data=attempt('final_native_data_read',lambda:remember(run/'slab-final.dat.h5'))
 fields=[];records=[];native_iterations=[];native_scaled=[]
 if geometry is not None and case:
  attempt('exact_canonical2000quad_geometry',lambda:geometry_check(geometry))
  mapped=attempt('native_topology_allfacegroups_and_cell_mapping',lambda:topology(case,geometry))
  if mapped:
   mapping,zones=mapped;evidence['native_face_groups']={name:b-a for name,(a,b) in zones.items()}
   def initial():
    path=remember(run/'snapshots/step-000000.dat.h5');t=native_time(path);assert t['time-step']==0 and t['flow-time']==0 and t['physical-time-step']==2
    T=temperatures(path,mapping);assert np.max(abs(T-300))<1e-10
    return T,dict(step=0,time_s=0.,native_dt_s=t['physical-time-step'],native_time_index=0,native_iteration=0,initial_unsolved=True,**state_metrics(T,geometry,0))
   value=attempt('native_initial_time0_uniform300K',initial)
   previous_energy=0.;previous_iteration=0
   if value:fields.append(value[0]);records.append(value[1]);previous_energy=value[1]['energy_J']
   for step in range(1,101):
    def audit_step(step=step):
     data=remember(run/f'snapshots/step-{step:06d}.dat.h5');report=remember(run/f'heat-reports/step-{step:06d}.txt')
     T,row=step_metrics(data,report,mapping,zones,geometry,step,previous_energy,previous_iteration)
     index,scaled=energy_residuals(data);assert np.array_equal(index,np.arange(1,len(index)+1))
     if native_iterations:assert np.array_equal(index[:len(native_iterations)],native_iterations) and np.array_equal(scaled[:len(native_scaled)],native_scaled)
     return T,row,index,scaled
    value=attempt(f'native_step{step:03d}_time_T_heatW_energy_residual',audit_step)
    if value:
     T,row,index,scaled=value;fields.append(T);records.append(row);previous_energy=row['energy_J'];previous_iteration=row['native_iteration'];native_iterations=index.tolist();native_scaled=scaled.tolist()
    else:break
   checks['all101_actual_states_audited']=len(fields)==len(records)==101
   if len(records)==101:
    evidence['maximum_global_Euler_energy_relative']=max(row['energy_balance_relative'] for row in records[1:]);evidence['maximum_step_terminal_energy_residual']=max(row['energy_residual'] for row in records[1:]);evidence['max_y_spread_K']=max(row['max_y_spread_K'] for row in records)
    evidence['final_metrics']=records[-1];evidence['early_transient_analytic_errors_reported_without_alltime_accuracy_claim']=True
    checks['final_analytic_max_error_le0p5K']=records[-1]['analytic_max_error_K']<=.5
    checks['final_analytic_weighted_RMS_error_le0p25K']=records[-1]['analytic_rms_error_K']<=.25
    checks['allsteps_energy_residual_le1e_10']=evidence['maximum_step_terminal_energy_residual']<=1e-10
    checks['allsteps_global_Euler_energy_balance_le1e_5']=evidence['maximum_global_Euler_energy_relative']<=1e-5
    def terminal_match():
     last=run/'snapshots/step-000100.dat.h5';assert native_time(final_data)==native_time(last)
     assert np.array_equal(temperatures(final_data,mapping),fields[-1])
     f_i,f_r=energy_residuals(final_data);assert np.array_equal(f_i,native_iterations) and np.array_equal(f_r,native_scaled)
     final_mapping,final_zones=topology(final_case,geometry);assert np.array_equal(mapping,final_mapping) and zones==final_zones
     return True
    attempt('final_saved_case_data_match_actual_step100',terminal_match)
    def history_match():
     history=load('time-history.json');assert history==records and len(history)==101;return True
    attempt('driver_history_recomputed_from_all_native_states',history_match)
    times=np.array([row['time_s'] for row in records]);Q=np.zeros((101,4));Q[1:]=[[row['wall_heat_in_W'][name] for name in WALLS] for row in records[1:]]
    np.savez_compressed(out/'fluent-fields-common.npz',**geometry,times_s=times,temperature_K=np.array(fields),native_step=np.arange(101),native_dt_s=np.array([row['native_dt_s'] for row in records]),wall_names=np.array(WALLS),wall_heat_in_W=Q,wall_heat_solved=np.arange(101)>0,energy_J=np.array([row['energy_J'] for row in records]),energy_balance_relative=np.array([row['energy_balance_relative'] for row in records[1:]]),energy_residual=np.array([row['energy_residual'] for row in records[1:]]),analytic_max_error_K=np.array([row['analytic_max_error_K'] for row in records]),analytic_rms_error_K=np.array([row['analytic_rms_error_K'] for row in records]))
    write_json(out/'fluent-time-metrics.json',records)
    csvrows=[]
    for row in records:
     flat={k:v for k,v in row.items() if not isinstance(v,dict)}
     if row['step']>0:flat.update({'Qin_'+name+'_W':row['wall_heat_in_W'][name] for name in WALLS})
     csvrows.append(flat)
    with (out/'fluent-time-metrics.csv').open('x',newline='',encoding='utf-8') as f:
     writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in csvrows for k in row)));writer.writeheader();writer.writerows(csvrows)
    with (out/'fluent-energy-residuals.csv').open('x',newline='',encoding='utf-8') as f:
     writer=csv.writer(f);writer.writerow(['native_iteration','scaled_energy_residual']);writer.writerows(zip(native_iterations,native_scaled))
 def transcript_check():
  logs=list(run.glob('fluent-*.trn'));assert len(logs)==1
  startup=remember(logs[0]).read_text(encoding='utf-8',errors='replace');text=remember(run/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace')
  assert 'student version' in startup.lower() and len(re.findall(r'^n\d+\*?\s+',startup,re.M))==2
  assert not re.search(r'Divergence detected|floating point exception|fatal error',text,re.I)
  # The only solved equation is energy; match each saved endpoint to actual printed residual.
  printed={int(m.group(1)):float(m.group(2)) for m in re.finditer(r'^\s*(\d+)\s+('+NUMBER+r')\s+\d+:\d+:\d+\s+\d+\s*$',text,re.M)}
  assert len(records)==101
  for row in records[1:]:
   value=printed[row['native_iteration']];assert abs(value-row['energy_residual'])<=6e-5*max(abs(row['energy_residual']),1e-300),'Printed/native terminal thermal residual mismatch'
  return dict(student_version_verified=True,native_compute_processes=2,actual_terminal_iterations=[row['native_iteration'] for row in records[1:]])
 evidence['native_licence_and_residual_transcript']=attempt('native_two_process_student_no_fatal_and_all100_printed_residuals',transcript_check)
 attempt('native_temperature_png_present',lambda:remember(run/'temperature.png'))
 checks['all_hashed_native_inputs_unchanged']=all(sha(path)==digest for path,digest in before.items())
 shutil.copyfile(__file__,out/'executed-analyser.py')
 for name,p in source_files().items():shutil.copyfile(p,out/name)
 report=dict(kind='Saved-native Fluent transient solid slab audit',checks=checks,failed_checks=[k for k,v in checks.items() if not v],all_checks_passed=all(checks.values()),errors=errors,evidence=evidence,physical_validation=False,mesh_independence=False,time_step_independence=False,all_cell_thermal_flux_conservation_claim=False,heat_power_note='Native integrated W inward sign crosschecked report/HDF/storage, no additional area multiplication',initial_heat_note='Initial state has no solved heat balance; NPZ initial Q placeholder0 is explicitly masked wall_heat_solved=False',analyser_sha256=sha(__file__),helper_sha256={name:sha(p) for name,p in source_files().items()})
 write_json(out/'fluent-checks.json',report);write_json(out/'fluent-source-manifest.json',dict(inputs=[dict(path=str(p),sha256=h) for p,h in before.items()],all_inputs_unchanged=checks['all_hashed_native_inputs_unchanged'],outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}))
 print(json.dumps(dict(all_checks_passed=report['all_checks_passed'],failed_checks=report['failed_checks'],errors=errors,evidence=evidence),allow_nan=False));return 0 if report['all_checks_passed'] else 1
if __name__=='__main__':raise SystemExit(main())
