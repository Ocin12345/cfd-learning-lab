"""Read-only saved-native audit of the single-step heated-channel Fluent pilot."""
from pathlib import Path
import argparse,csv,json,re,shutil,math
import h5py,numpy as np
from heated_channel_fluent_common_v001 import *
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();run=a.run.resolve();out=a.out.resolve()
 assert run.is_dir() and out.drive.upper()=='E:' and not out.exists() and not out.is_relative_to(run);out.mkdir(parents=True);checks={};errors={};evidence={};before={}
 def remember(p):p=Path(p).resolve();assert p.is_file(),str(p);before.setdefault(p,sha(p));return p
 def load(n):return json.loads(remember(run/n).read_text(encoding='utf-8-sig'))
 def attempt(name,fn):
  try:value=fn();checks[name]=True;return value
  except Exception as e:checks[name]=False;errors[name]=type(e).__name__+': '+str(e);return None
 status=attempt('native_status_read',lambda:load('status.json')) or {};summary=attempt('native_summary_read',lambda:load('summary.json')) or {};source=attempt('source_provenance_read',lambda:load('code-provenance.json')) or {};auth=attempt('root_authorization_read',lambda:load('authorization.json')) or {};declared=attempt('declared_targets_read',lambda:load('declared-targets.json')) or {};quality=attempt('native_quality_read',lambda:load('geometry-quality.json')) or {}
 mesh=Path(source.get('mesh',str(run.parent/'mesh/heated-channel.msh'))).resolve();geometry_path=mesh.parent/'geometry.npz'
 def provenance():
  assert auth['action']=='run_heated_channel_fluent_pilot_v001' and auth['root']=='root' and auth['approved'] is True and auth['allowed_mode']=='Run' and auth['pilot_thermal_steps']==1 and Path(auth['out']).resolve()==run
  assert source['source_sha256']==auth['source_sha256']==sha(remember(run/'executed-script.py'));assert source['mesh_sha256']==auth['mesh_sha256']==sha(remember(mesh))
  assert source['geometry_npz_sha256']==sha(remember(geometry_path)) and source['mesh_setup_sha256']==sha(remember(mesh.parent/'setup.json'));assert source['targets_sha256']==auth['targets_sha256']==TARGETS_SHA==sha(remember(run/'TARGETS.md'));assert declared==TARGETS
  assert source['helper_sha256']==auth['helper_sha256']
  for name,path in source_files().items():assert source['helper_sha256'][name]==sha(remember(run/'sources'/name))==sha(remember(path))
  assert source['snapshot_before_solver_launch'] and source['available_ram_before_launch_bytes']>=4*1024**3 and source['free_disk_before_launch_bytes']>=10*1024**3
  return dict(source_sha256=source['source_sha256'],helper_sha256=source['helper_sha256'],mesh_sha256=source['mesh_sha256'],targets_sha256=TARGETS_SHA)
 evidence['provenance']=attempt('source_mesh_one_step_authorization_and_targets_bound',provenance)
 checks['pilot_completed_exactly2000flow_plus1thermal_step']=status.get('stage')=='pilot_completed' and summary.get('actual_flow_iterations')==2000 and summary.get('actual_thermal_steps')==1 and summary.get('actual_time_s')==.005 and summary.get('flow_native_saved') is True and summary.get('pilot_native_saved') is True
 checks['production_not_authorized']=summary.get('production_authorized') is False;checks['no_export_errors']=summary.get('export_errors')==[]
 checks['native_mesh_positive_volume0p1_quality0p15']=quality.get('minimum_cell_volume_m3',0)>0 and abs(quality.get('total_cell_volume_m3',0)-.1)<1e-9 and quality.get('minimum_orthogonal_quality',-1)>=.15
 def settings(name,thermal):
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
 for name in ('flow-loaded-settings.json','flow-initialized-settings.json','flow-final-settings.json'):attempt('native_physical_settings_'+name,lambda name=name:settings(name,False))
 for name in ('thermal-before-patch-settings.json','thermal-initial-settings.json','thermal-step000001-settings.json'):attempt('native_physical_settings_'+name,lambda name=name:settings(name,True))
 geometry=attempt('source_geometry_read',lambda:dict(np.load(remember(geometry_path),allow_pickle=False)))
 allstates={}
 if geometry is not None:
  attempt('exact2304quad_rectangular_geometry',lambda:geometry_check(geometry))
  for stem,thermal in (('flow-initial',False),('flow-iter1800',False),('flow-final',False),('thermal-before-patch',True),('thermal-initial',True),('thermal-step000001',True)):
   allstates[stem]=attempt('native_mapped_fields_'+stem,lambda stem=stem,thermal=thermal:snapshot(remember(run/(stem+'.cas.h5')),remember(run/(stem+'.dat.h5')),geometry,thermal=thermal))
  if allstates.get('flow-initial') is not None:checks['native_flow_initial_U_V_p_from_rest']=all(float(abs(allstates['flow-initial']['common'][k]).max())<=1e-12 for k in ('u','v','p'))
  for endpoint in (500,1000,1500,2000):
   attempt('checkpoint'+str(endpoint)+'_case',lambda endpoint=endpoint:remember(run/f'flow-iter{endpoint:04d}.cas.h5'));attempt('checkpoint'+str(endpoint)+'_data',lambda endpoint=endpoint:remember(run/f'flow-iter{endpoint:04d}.dat.h5'))
  if allstates.get('flow-final') is not None and allstates.get('flow-iter1800') is not None:
   def flow_audit():
    value=flow_gates(allstates['flow-final'],allstates['flow-iter1800'],geometry,run/'flow-final.dat.h5');assert value==load('flow-gates.json');return value
   flow=attempt('native_flow_metrics_recomputed_match_driver',flow_audit)
   if flow:
    evidence['flow']=flow['metrics'];evidence['native_geometric_normal_orientation']=allstates['flow-final']['face_orientation_proof'];checks.update(flow['checks'])
   for stem,label in (('thermal-before-patch','frozen-before-patch.json'),('thermal-initial','frozen-initial.json'),('thermal-step000001','frozen-step000001.json')):
    if allstates.get(stem) is not None:
     value=attempt('native_frozen_comparison_recomputed_'+stem,lambda stem=stem:frozen_compare(allstates['flow-final'],allstates[stem]))
     if value is not None:checks['every_U_V_p_and_face_flux_unchanged_'+stem]=value['passed'];evidence[stem+'_frozen']=value;attempt('driver_frozen_record_matches_'+stem,lambda value=value,label=label:True if value==load(label) else (_ for _ in ()).throw(AssertionError('Frozen record differs')))
  if allstates.get('thermal-initial') is not None and allstates.get('thermal-step000001') is not None:
   def thermal_audit():
    value=thermal_metrics(allstates['thermal-initial'],allstates['thermal-step000001'],geometry,run/'thermal-initial.dat.h5',run/'thermal-step000001.dat.h5',remember(run/'reports/step000001-total-heat.txt'),remember(run/'reports/step000001-mass.txt'),initial_total_report=remember(run/'reports/initial-total-heat.txt'));assert value==load('pilot-thermal-gates.json');return value
   thermal=attempt('native_pilot_thermal_metrics_recomputed_match_driver',thermal_audit)
   if thermal:
    def sensible_diagnostic():
     report=run/'reports/step000001-sensible-heat.txt'
     if not report.is_file():return dict(native_sensible_report_available=False,note='Explicit command activation proof retained; no sensible field or power invented')
     actual=flux_report(remember(report),'W');result=dict(native_sensible_report_available=True,native_sensible_report_W=actual)
     with h5py.File(run/'thermal-step000001.dat.h5','r') as f:available='SV_HEAT_FLUX_SENSIBLE' in f['results/1/phase-1/faces']
     result['distinct_native_sensible_HDF_field_available']=available
     if available:
      raw,powers,coverage=boundary_array(run/'thermal-step000001.dat.h5','SV_HEAT_FLUX_SENSIBLE',allstates['thermal-step000001']['zones']);result.update(native_sensible_HDF_W=powers,coverage=coverage)
     return result
    evidence['native_sensible_heat_diagnostic']=attempt('native_sensible_report_availability_recorded_without_assumed_semantics',sensible_diagnostic)
    checks.update(thermal['checks']);evidence['thermal']=thermal['metrics'];initial=allstates['thermal-initial'];final=allstates['thermal-step000001'];reference=allstates['flow-final'];t0=native_time(run/'thermal-initial.dat.h5');t1=native_time(run/'thermal-step000001.dat.h5');m=thermal['metrics']
    np.savez_compressed(out/'fluent-pilot-fields.npz',**geometry,temperature_K=np.array([initial['common']['T'],final['common']['T']]),native_enthalpy_J_kg=np.array([initial['common']['h'],final['common']['h']]),times_s=np.array([t0['flow-time'],t1['flow-time']]),native_step=np.array([t0['time-step'],t1['time-step']]),native_dt_s=np.array([t0['physical-time-step'],t1['physical-time-step']]),frozen_u=reference['common']['u'],frozen_v=reference['common']['v'],frozen_p=reference['common']['p'],frozen_face_mass_flux_kg_s=reference['flux'],face_c0=reference['c0'],face_c1=reference['c1'],native_face_area=reference['face_area'],raw_native_edge_area=reference['raw_face_area'],native_face_mid=reference['face_mid'],native_cell_to_common=reference['mapping'],physical_boundary_names=np.array(BOUNDARIES),native_raw_boundary_Qin_W=np.array([[m['native_raw_Qin_W'][n] for n in BOUNDARIES]]),native_reference_corrected_boundary_Qin_W_per_m=np.array([[m['native_reference_corrected_Qin_above300_W_per_m'][n] for n in BOUNDARIES]]),reconstructed_boundary_Qin_W_per_m=np.array([[m['reconstructed_PDE_Qin_above300_W_per_m'][n] for n in BOUNDARIES]]),energy_J_per_m=np.array([m['energy_initial_J_per_m'],m['energy_step1_J_per_m']]),storage_W_per_m=np.array([m['storage_W_per_m']]),outlet_bulk_temperature_step1_K=np.array([m['actual_outward_flux_weighted_outlet_bulk_temperature_K']]))
 def flow_hashes():
  hashes=load('frozen-flow-source-hashes.json');assert all(sha(remember(run/n))==h for n,h in hashes.items());assert load('frozen-flow-original-files-check.json')['all_unchanged'] is True;return hashes
 evidence['original_flow_native_files']=attempt('frozen_original_flow_hashes_unchanged',flow_hashes)
 def native_logs():
  startup=list(run.glob('fluent-*.trn'));assert len(startup)==1;text=remember(startup[0]).read_text(encoding='utf-8',errors='replace');transcript=remember(run/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace');assert 'student version' in text.lower() and len(re.findall(r'^n\d+\*?\s+',text,re.M))==2
  assert not re.search(r'api-set-var: the object is not active|Divergence detected|floating point exception|fatal error',transcript,re.I)
  printed={int(m.group(1)):[float(m.group(i)) for i in (2,3,4)] for m in re.finditer(r'^\s*(\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+\d+:\d+:\d+\s+\d+\s*$',transcript,re.M)}
  final,_=terminal_residuals(run/'flow-final.dat.h5',EQUATIONS)
  for name,value in zip(EQUATIONS,printed[2000]):assert abs(value-final[name])<=6e-5*max(final[name],1e-300)
  energy,energy_indices=terminal_residuals(run/'thermal-step000001.dat.h5',('energy',));last=energy_indices['energy'][-1]
  printed_energy={int(m.group(1)):float(m.group(2)) for m in re.finditer(r'^\s*(\d+)\s+('+NUMBER+r')\s+\d+:\d+:\d+\s+\d+\s*$',transcript,re.M)}
  assert abs(printed_energy[last]-energy['energy'])<=6e-5*max(energy['energy'],1e-300)
  return dict(student_version=True,native_compute_processes=2,actual_terminal_energy_iteration=last)
 evidence['native_licence_and_logs']=attempt('native_two_process_student_and_printed_flow_residuals_no_fatal',native_logs)
 for p in run.rglob('*'):
  if p.is_file():remember(p)
 checks['all_hashed_native_inputs_unchanged']=all(sha(p)==h for p,h in before.items());shutil.copyfile(__file__,out/'executed-analyser.py')
 for n,p in source_files().items():shutil.copyfile(p,out/n)
 report=dict(kind='Heated-channel one-step Fluent saved-native pilot audit',checks=checks,all_checks_passed=all(checks.values()),failed_checks=[k for k,v in checks.items() if not v],errors=errors,evidence=evidence,production_authorized=False,physical_validation=False,mesh_independence=False,time_step_independence=False,analyser_sha256=sha(__file__))
 write_json(out/'fluent-pilot-checks.json',report);write_json(out/'fluent-pilot-source-manifest.json',dict(inputs=[dict(path=str(p),sha256=h) for p,h in before.items()],all_inputs_unchanged=checks['all_hashed_native_inputs_unchanged'],outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}));print(json.dumps(report,allow_nan=False));return 0 if report['all_checks_passed'] else 1
if __name__=='__main__':raise SystemExit(main())
