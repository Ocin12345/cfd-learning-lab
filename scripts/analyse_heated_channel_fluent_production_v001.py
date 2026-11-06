"""Read-only all600-step native heated-channel production audit. No solver launch."""
from pathlib import Path
import argparse,json,re,shutil
import numpy as np
from heated_channel_fluent_production_common_v001 import *

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();run=a.run.resolve();out=a.out.resolve()
 assert run.is_dir() and out.drive.upper()=='E:' and not out.exists() and not out.is_relative_to(run);out.mkdir(parents=True);checks={};errors={};evidence={};before={}
 def remember(path):
  path=Path(path).resolve();assert path.is_file(),str(path);before.setdefault(path,sha(path));return path
 def load(name):return json.loads(remember(run/name).read_text(encoding='utf-8-sig'))
 def attempt(name,fn):
  try:value=fn();checks[name]=True;return value
  except Exception as e:checks[name]=False;errors[name]=type(e).__name__+': '+str(e);return None
 status=attempt('native_status',lambda:load('status.json')) or {};summary=attempt('native_summary',lambda:load('summary.json')) or {};source=attempt('native_provenance',lambda:load('code-provenance.json')) or {};auth=attempt('root_authorization',lambda:load('authorization.json')) or {};declared=attempt('declared_targets',lambda:load('declared-targets.json')) or {}
 mesh=Path(source.get('mesh',str(run.parent/'mesh/heated-channel.msh')));gp=mesh.parent/'geometry.npz'
 def provenance():
  assert auth['action']=='continue_heated_channel_fluent_production_v001' and auth['root']=='root' and auth['approved'] is True and auth['allowed_mode']=='Run' and Path(auth['out']).resolve()==run
  assert auth['total_thermal_steps']==source['total_thermal_steps']==600 and auth['additional_thermal_steps']==source['additional_thermal_steps']==599 and auth['additional_flow_iterations']==source['additional_flow_iterations']==0 and auth['monitor_history_retention_n_save']==source['monitor_history_retention_n_save']==40000
  assert auth['source_sha256']==source['source_sha256']==sha(remember(run/'executed-script.py')) and auth['helper_sha256']==source['helper_sha256']
  for n,h in source['helper_sha256'].items():assert sha(remember(run/'sources'/n))==h==sha(remember(Path(__file__).with_name(n)))
  assert auth['targets_sha256']==source['targets_sha256']==TARGETS_SHA==sha(remember(run/'TARGETS.md')) and declared==TARGETS
  assert auth['mesh_sha256']==source['mesh_sha256']==sha(remember(mesh)) and source['geometry_npz_sha256']==sha(remember(gp)) and source['mesh_setup_sha256']==sha(remember(mesh.parent/'setup.json'))
  manifest=load('parent-native-inputs.json');parent=Path(manifest['parent']);assert parent.resolve()==Path(auth['parent']).resolve()==Path(source['continuation_parent']).resolve()
  assert auth['parent_manifest_sha256']==source['parent_manifest_sha256']==sha(remember(run/'parent-native-inputs.json'))
  for n,h in manifest['files'].items():assert sha(remember(parent/n))==h
  assert auth['parent_case_sha256']==source['parent_case_sha256']==sha(remember(parent/'thermal-step000001.cas.h5')) and auth['parent_data_sha256']==source['parent_data_sha256']==sha(remember(parent/'thermal-step000001.dat.h5'))
  assert auth['selected_pilot_audit_sha256']==source['selected_pilot_audit_sha256']==manifest['selected_pilot_audit_sha256']==sha(remember(run/'accepted-pilot-audit.json'));accepted=load('accepted-pilot-audit.json');assert accepted['all_checks_passed'] and accepted['errors']=={}
  for n,row in load('inherited-native-files.json').items():assert sha(remember(run/n))==row['sha256']==sha(remember(row['source']))
  assert source['snapshot_before_solver_launch'] and source['available_ram_before_launch_bytes']>=4*1024**3 and source['free_disk_before_launch_bytes']>=10*1024**3
  return dict(source_sha256=source['source_sha256'],helper_sha256=source['helper_sha256'],parent_manifest_sha256=source['parent_manifest_sha256'],selected_pilot_audit_sha256=source['selected_pilot_audit_sha256'],parent=str(parent))
 evidence['provenance']=attempt('source_parent_pilot_receipt_mesh_targets_bound',provenance)
 checks['completed600_total599_additional_no_flow']=status.get('stage') in ('completed','completed_with_export_errors') and summary.get('actual_flow_iterations')==2000 and summary.get('additional_flow_iterations')==0 and summary.get('actual_thermal_steps')==600 and summary.get('additional_thermal_steps_executed')==599 and summary.get('native_saved_states')==601 and abs(summary.get('actual_time_s',-1)-3.)<=1e-10
 checks['native_graphics_and_other_exports_no_errors']=summary.get('export_errors')==[]
 q=attempt('native_quality',lambda:load('geometry-quality.json')) or {};checks['positive_native_mesh_volume0p1_quality0p15']=q.get('minimum_cell_volume_m3',0)>0 and abs(q.get('total_cell_volume_m3',0)-.1)<1e-9 and q.get('minimum_orthogonal_quality',0)>=.15
 def settings(name):
  thermal=True
  s=load(name);solver=s['general']['solver'];assert solver['type']=='pressure-based' and solver['two_dim_space']=='planar' and solver['time']==('transient' if thermal else 'steady')
  assert s['models']['viscous']['model']=='laminar' and s['models']['energy']['enabled'] is thermal;assert s['general']['operating_conditions']['gravity']['enable'] is False and s['general']['operating_conditions']['operating_pressure']==0
  assert s['cell_zone_conditions']['fluid']['fluid']['general']['material']=='air'
  for key,value in (('density',1.),('viscosity',.01)):
   assert s['material'][key]['option']=='value' and s['material'][key]['value']==value
  assert s['references']['depth']==1.;assert s['inlet']['momentum']['velocity_specification_method']=='Components' and [i['value'] for i in s['inlet']['momentum']['velocity_components']]==[1,0] and s['outlet']['momentum']['gauge_pressure']['value']==0
  for n in WALLS:assert s['walls'][n]['momentum']['wall_motion']=='Stationary Wall' and s['walls'][n]['momentum']['shear_condition']=='No Slip'
  if thermal:
   assert s['work_control_proof']['radiation_model']['verified_off_or_inactive'] is True
   if s['work_control_proof']['radiation_model']['active']:assert s['models']['radiation']['model']=='none'
   for key,value in (('specific_heat',1000.),('thermal_conductivity',10.)):assert s['material'][key]['option']=='value' and s['material'][key]['value']==value
   assert s['models']['energy']['inlet_diffusion'] is True
   for key in ('pressure_work','kinetic_energy','viscous_dissipation','viscous_heating'):assert s['work_control_proof'][key]['verified_off_or_inactive'] is True
   assert s['inlet']['thermal']['temperature']['value']==300. and s['outlet']['thermal']['backflow_total_temperature']['value']==300.
   for n in WALLS:assert s['walls'][n]['thermal']['thermal_condition']=='Temperature' and s['walls'][n]['thermal']['temperature']['value']==350.
   eq=s['controls']['equations'];assert any(v for v in eq.values()) and all(('temperature' in str(k).lower() or 'energy' in str(k).lower()) for k,v in eq.items() if v)
   assert s['methods']['transient_formulation']=='unsteady-1st-order' and s['methods']['spatial_discretization']['discretization_scheme']['temperature']=='first-order-upwind';run_state=s['run_calculation'];assert run_state['parameters']['time_step_size']==.005 and run_state['parameters']['max_iter_per_time_step']==60 and run_state['transient_controls']['type']=='Fixed' and run_state['transient_controls']['method']=='User-Specified'
   assert s['residual']['equations']['energy']['check_convergence'] is True and s['residual']['equations']['energy']['absolute_criteria']==1e-8
  else:
   assert s['methods']['p_v_coupling']['flow_scheme']=='SIMPLE' and s['methods']['spatial_discretization']['gradient_scheme']=='least-square-cell-based' and s['methods']['spatial_discretization']['discretization_scheme']['mom']=='second-order-upwind' and s['methods']['spatial_discretization']['discretization_scheme']['pressure']=='second-order'
   for n in EQUATIONS:assert s['residual']['equations'][n]['check_convergence'] is False
  return True

 for n in ('production-loaded-settings.json','production-start-settings.json','thermal-final-settings.json'):attempt('native_physical_settings_'+n,lambda n=n:settings(n))
 def retention():
  loaded=load('production-loaded-settings.json');start=load('production-start-settings.json');expected=json.loads(json.dumps(loaded['residual']));expected['options']['n_save']=40000;assert start['residual']==expected and load('thermal-final-settings.json')['residual']==expected
  assert load('monitor-retention-only-change.json')['all_other_residual_settings_identical'] is True and load('parent-physical-settings-comparison.json')['all_checks_passed'] is True
  return dict(before_n_save=loaded['residual']['options']['n_save'],after_n_save=40000,convergence_settings_unchanged=True)
 evidence['retention']=attempt('only_monitor_history_retention_changed',retention)
 g=attempt('source_geometry',lambda:dict(np.load(remember(gp),allow_pickle=False)));temperature=np.full((601,2304),np.nan);enthalpy=np.full((601,2304),np.nan);times=np.full(601,np.nan);steps=np.full(601,np.nan);dt=np.full(601,np.nan);frozen_delta=np.full((601,4),np.nan);metrics=[];residual_i=[];residual_v=[]
 if g is not None:
  attempt('canonical2304_quad_geometry',lambda:geometry_check(g));flow=snapshot(remember(run/'flow-final.cas.h5'),remember(run/'flow-final.dat.h5'),g);baseline=snapshot(remember(run/'flow-iter1800.cas.h5'),remember(run/'flow-iter1800.dat.h5'),g);fg=flow_gates(flow,baseline,g,run/'flow-final.dat.h5');checks.update(fg['checks']);evidence['flow']=fg['metrics'];assert fg==load('flow-gates.json')
  initial=snapshot(remember(run/'thermal-initial.cas.h5'),remember(run/'states/step-000000.dat.h5'),g,True);previous=None
  static=remember(run/'production-start.cas.h5');remember(run/'inherited-step000001.cas.h5')
  for step in range(601):
   data=remember(run/'states'/f'step-{step:06d}.dat.h5');case=run/'thermal-initial.cas.h5' if step==0 else (run/'inherited-step000001.cas.h5' if step==1 else static)
   def audit_state():
    current=snapshot(case,data,g,True);t=native_time(data);frozen=frozen_compare(flow,current)
    assert abs(t['flow-time']-step*.005)<=1e-10 and t['time-step']==step and abs(t['physical-time-step']-.005)<=1e-12
    temperature[step]=current['common']['T'];enthalpy[step]=current['common']['h'];times[step]=t['flow-time'];steps[step]=t['time-step'];dt[step]=t['physical-time-step'];frozen_delta[step]=[frozen['maximum_absolute_changes_SI'][n] for n in ('u','v','p','all_face_mass_flux')]
    checks[f'actual_state{step:06d}_frozen_flow_all_native_values']=frozen['passed']
    if step:
     heat=remember(run/'heat-reports'/f'step-{step:06d}.txt');mass=remember(run/'mass-reports'/f'step-{step:06d}.txt');transcript=remember(run/('pilot-solver-transcript.txt' if step==1 else 'solver-transcript.txt'))
     result=thermal_metrics(initial,current,g,run/'states/step-000000.dat.h5',data,heat,mass,initial_total_report=remember(run/'heat-reports/step-000000.txt'),previous=previous,previous_data=run/'states'/f'step-{step-1:06d}.dat.h5',expected_step=step,transcript_path=transcript)
     assert result==load('step-metrics/'+f'step-{step:06d}.json'),'Actual recomputed metrics differ from driver'
     for key,value in result['checks'].items():checks[f'step{step:06d}_'+key]=value
     m=result['metrics'];metrics.append(m);proof=m['native_actual_thermal_residual_proof'];residual_i.extend(proof['indices']);residual_v.extend(proof['scaled'])
    return current
   current=attempt(f'actual_state{step:06d}_native_time_fields_flux_reports_metrics',audit_state)
   if current is not None:previous=current
  checks['all601_native_states_and600_step_metric_records']=len(metrics)==600 and np.isfinite(temperature).all() and np.isfinite(enthalpy).all() and np.isfinite(frozen_delta).all()
  if len(metrics)==600:
   energy=np.sum(g['areas'][None,:]*1000*(temperature-300),axis=1);qkeys=('native_raw_Qin_W','native_reference_corrected_Qin_above300_W_per_m','reconstructed_PDE_Qin_above300_W_per_m');powers={key:np.array([[m[key][n] for n in BOUNDARIES] for m in metrics]) for key in qkeys}
   np.savez_compressed(out/'fluent-thermal-fields.npz',**g,temperature_K=temperature,native_enthalpy_J_kg=enthalpy,times_s=times,native_step=steps,native_dt_s=dt,frozen_u=flow['common']['u'],frozen_v=flow['common']['v'],frozen_p=flow['common']['p'],frozen_face_mass_flux_kg_s=flow['flux'],face_c0=flow['c0'],face_c1=flow['c1'],native_face_area=flow['face_area'],raw_native_edge_area=flow['raw_face_area'],native_face_mid=flow['face_mid'],native_cell_to_common=flow['mapping'],physical_boundary_names=np.array(BOUNDARIES),native_raw_boundary_Qin_W=powers[qkeys[0]],native_reference_corrected_boundary_Qin_W_per_m=powers[qkeys[1]],reconstructed_boundary_Qin_W_per_m=powers[qkeys[2]],energy_J_per_m=energy,storage_W_per_m=np.array([m['storage_W_per_m'] for m in metrics]),native_Euler_relative_defect=np.array([m['native_Euler_relative_defect'] for m in metrics]),reconstructed_PDE_relative_defect=np.array([m['reconstructed_PDE_relative_defect'] for m in metrics]),energy_final_scaled_residual=np.array([m['native_energy_final_scaled_residual'] for m in metrics]),actual_inner_iterations=np.array([m['native_actual_thermal_residual_proof']['actual_inner_iterations'] for m in metrics]),residual_iteration=np.array(residual_i),residual_energy_scaled=np.array(residual_v),frozen_field_maximum_absolute_changes_SI=frozen_delta,outlet_bulk_temperature_K=np.array([300.]+[m['actual_outward_flux_weighted_outlet_bulk_temperature_K'] for m in metrics]))
   evidence['thermal_summary']=dict(actual_steps=600,saved_states=601,first_time_s=times[0],last_time_s=times[-1],temperature_min_K=float(temperature.min()),temperature_max_K=float(temperature.max()),maximum_energy_final_scaled_residual=max(m['native_energy_final_scaled_residual'] for m in metrics),maximum_native_Euler_relative_defect=max(m['native_Euler_relative_defect'] for m in metrics),maximum_reconstructed_PDE_relative_defect=max(m['reconstructed_PDE_relative_defect'] for m in metrics),minimum_inner_iterations=min(m['native_actual_thermal_residual_proof']['actual_inner_iterations'] for m in metrics),maximum_inner_iterations=max(m['native_actual_thermal_residual_proof']['actual_inner_iterations'] for m in metrics),actual_total_thermal_inner_iterations=len(residual_i),native_h300_J_kg=metrics[0]['native_h300_J_kg'],maximum_native_h_cpDeltaT_error_J_kg=max(m['native_h_cpDeltaT_max_error_J_kg'] for m in metrics),initial_unsolved_heat_diagnostic=metrics[0]['native_initial_unsolved_heat_diagnostic'],maximum_frozen_field_changes_SI=dict(zip(('u','v','p','all_face_mass_flux'),np.max(frozen_delta,axis=0).tolist())),final_energy_J_per_m=float(energy[-1]),final_outlet_bulk_temperature_K=metrics[-1]['actual_outward_flux_weighted_outlet_bulk_temperature_K'],final_native_Qin_above300_W_per_m=metrics[-1]['native_reference_corrected_Qin_above300_W_per_m'])
  def final_native():
   final=snapshot(remember(run/'thermal-final.cas.h5'),remember(run/'thermal-final.dat.h5'),g,True);assert np.array_equal(final['common']['T'],temperature[-1]) and np.array_equal(final['common']['h'],enthalpy[-1]) and frozen_compare(flow,final)['passed'];assert native_time(run/'thermal-final.dat.h5')['time-step']==600;return True
  attempt('native_final_case_data_match_last_actual_saved_state',final_native)
 def logs():
  startup=list(run.glob('fluent-*.trn'));assert len(startup)==1;text=remember(startup[0]).read_text(encoding='utf-8',errors='replace');transcript=remember(run/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace');assert 'student version' in text.lower() and len(re.findall(r'^n\d+\*?\s+',text,re.M))==2;assert not re.search(r'api-set-var: the object is not active|Divergence detected|floating point exception|fatal error',transcript,re.I);return dict(student_version=True,actual_compute_processes=2)
 evidence['licence']=attempt('actual_hidden_two_process_student_and_no_fatal',logs)
 for x in run.rglob('*'):
  if x.is_file():remember(x)
 checks['all_hashed_inputs_unchanged']=all(sha(x)==h for x,h in before.items());shutil.copyfile(__file__,out/'executed-analyser.py')
 for n,x in source_files().items():shutil.copyfile(x,out/n)
 report=dict(kind='Heated-channel all600-step Fluent saved-native production audit',checks=checks,all_checks_passed=all(checks.values()),failed_checks=[k for k,v in checks.items() if not v],errors=errors,evidence=evidence,analyser_sha256=sha(__file__),physical_validation=False,mesh_independence=False,time_step_independence=False)
 write_json(out/'fluent-thermal-checks.json',report);write_json(out/'fluent-thermal-source-manifest.json',dict(inputs=[dict(path=str(x),sha256=h) for x,h in before.items()],all_inputs_unchanged=checks['all_hashed_inputs_unchanged'],outputs={x.name:sha(x) for x in out.iterdir() if x.is_file()}));print(json.dumps(dict(check_count=len(checks),all_checks_passed=report['all_checks_passed'],failed_checks=report['failed_checks'],errors=errors,evidence=evidence),allow_nan=False));return 0 if report['all_checks_passed'] else 1
if __name__=='__main__':raise SystemExit(main())
