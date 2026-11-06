"""Read-only annular Couette Fluent PRODUCTION audit; fresh private outputs; no CFD."""
from pathlib import Path
import argparse,ctypes,datetime,json,math
import numpy as np
import annular_couette_fluent_common_v003 as c
def available_memory():
 class MS(ctypes.Structure):
  _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ['total','available','total_page','available_page','total_virtual','available_virtual','available_extended']]
 m=MS();m.length=ctypes.sizeof(m);c.require(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)),'Cannot query RAM');return m.available
def settings_gate(st,setup):
 c.require(st['general']['solver']['type']=='pressure-based' and st['general']['solver']['two_dim_space']=='planar' and st['general']['solver']['time']=='transient','Wrong native Fluent formulation')
 op=st['general']['operating_conditions'];c.require(op['operating_pressure']==0 and op['gravity']['enable'] is False,'Wrong native pressure/gravity')
 c.require(np.allclose(op['reference_pressure_location'],setup['pressure_reference_xy_m'],rtol=0,atol=1e-14),'Pressure reference is not sourcecell0 centre')
 c.require(st['viscous']['model']=='laminar' and st['energy']['enabled'] is False and st['dynamic_mesh']['enabled'] is False,'Unexpected physical model')
 c.require(st['material']['density']['value']==1 and st['material']['viscosity']['value']==.000125 and st['references']['depth']==1.,'Wrong material/depth')
 zone=st['cell_zones']['fluid'];c.require(zone['reference_frame']['frame_motion'] is False and zone['mesh_motion']['enable'] is False,'Rotating or moving fluid cell zone')
 inner=st['walls']['inner']['momentum'];outer=st['walls']['outer']['momentum'];c.require(inner['wall_motion']=='Moving Wall' and inner['relative']=='Absolute' and inner['velocity_spec']=='Rotational' and inner['shear_condition']==outer['shear_condition']=='No Slip' and outer['wall_motion']=='Stationary Wall','Wrong actual wall motion')
 methods=st['methods'];c.require(methods['p_v_coupling']['flow_scheme']=='SIMPLE' and methods['transient_formulation']=='unsteady-1st-order' and methods['spatial_discretization']['discretization_scheme']['mom']=='first-order-upwind','Wrong matched temporal/momentum order')
 c.require(all(st['residual']['equations'][name]['absolute_criteria']==1e-6 and st['residual']['equations'][name]['check_convergence'] is True for name in c.EQUATIONS),'Changed own solver residual target')
 return True
def analyse(run,out):
 run=Path(run).resolve();out=Path(out).resolve();c.require(not out.exists() and out.drive.upper()=='E:','Fresh private E audit required');c.require(available_memory()>=int(2.5*1024**3),'Read-only analysis requires unchanged2.5GiB RAM')
 out.mkdir(parents=True);before={};checks={};evidence={}
 def remember(p,expected=None):
  p=Path(p).resolve();digest=c.sha(p);c.require(expected is None or digest==expected,'Changed bound source '+str(p));before.setdefault(str(p),digest);return p
 def load(name):return c.read(remember(run/name))
 source=load('code-provenance.json');auth=load('authorization.json');status=load('status.json');rows=load('states.json');setup=load('mesh-setup.json');targets=load('declared-targets.json')
 c.require(auth['approved'] is True and auth['root']=='root' and auth['action']=='annular_couette_fluent_production_v001' and auth['phase']=='production-395-additional-steps' and Path(auth['out']).resolve()==run,'Wrong actual root pilot authorization')
 remember(run/'executed-production.py',auth['source_sha256']);remember(run/'executed-common.py',auth['common_sha256']);remember(run/'TARGETS.md',c.TARGETS_SHA)
 c.require(source['source_sha256']==auth['source_sha256'] and source['available_ram_bytes']>=4*1024**3 and source['free_E_bytes']>=10*1024**3 and source['processor_count']==2,'Launch source/resources/ranks proof differs')
 for row in source['inputs']:remember(row['path'],row['sha256'])
 c.require(targets==c.TARGETS,'Changed frozen numerical targets');mesh=next(Path(r['path']) for r in source['inputs'] if r['path'].endswith('.msh'))
 remember(mesh,c.MESH_SHA);remember(mesh.parent/'geometry.npz',c.GEOMETRY_SHA);remember(mesh.parent/'setup.json',c.SETUP_SHA)
 g=np.load(mesh.parent/'geometry.npz');c.geometry_check(g);remember(__file__);remember(c.__file__)
 resume=load('resume-preservation-proof.json');remember(run/resume['case'],resume['case_sha256']);remember(run/resume['data'],resume['data_sha256'])
 checks['actual_pilot_fields_preserved_before_continuation']=resume['actual_resumed_step']==5 and resume['time_s']==.25 and resume['no_physical_step_executed'] is True and resume['no_initialization_called'] is True and resume['all_fields_and_native_mapping_retained'] is True and all(v<=1e-12 for v in resume['maximum_absolute_field_changes_SI'].values())
 checks['source_and_authorization_bound']=True;checks['exact_pilot_state_coverage']=[r['step'] for r in rows]==list(range(401));checks['driver_completed']=status['stage']=='production_completed' and status['actual_additional_steps']==395 and status['actual_total_steps']==400 and status['actual_end_time_s']==20.
 quality=load('geometry-quality.json');checks['native_mesh_quality']=quality['checked_before_flow'] is True and quality['minimum_orthogonal_quality']>=.15 and quality['minimum_native_cell_volume']>0
 settings=load('final-settings.json');checks['actual_native_settings']=settings_gate(settings,setup);rotation=load('native-physical-readback.json');speed=rotation['rotation']['rotation_speed']
 checks['live_rotation_activity_and_units']=speed['active'] is True and abs(speed['rad_s']-1)<1e-12 and bool(speed['native_units']) and rotation['rotation']['origin']['active'] is True and np.max(abs(np.asarray(rotation['rotation']['origin']['state'],float)))<1e-12
 transcript=remember(run/'solver-transcript.txt');states=[];state_results=[];residuals=[];terminal_indices=[];paths=[];source_indices=None
 for item in rows:
  step=item['step'];case=remember(run/item['case'],item['case_sha256']);data=remember(run/item['data'],item['data_sha256']);s=c.read_state(case,data,g);gate=c.state_gates(s,step)
  c.require(gate==item['native_gates'] if step==0 else all(gate['checks'][k]==item['native_gates']['checks'][k] for k in gate['checks']),'Driver/native replay gate discrepancy')
  if states:c.require(np.array_equal(s['mapping'],states[0]['mapping']) and np.array_equal(s['face_mapping'],states[0]['face_mapping']) and np.array_equal(s['face_nodes'],states[0]['face_nodes']) and np.array_equal(s['native_cell_centres'],states[0]['native_cell_centres']),'Actual native fixed mesh/mapping changed')
  if step:
   terminal_log=remember(Path(source['accepted_pilot_run'])/'solver-transcript.txt') if step<=5 else transcript
   residual,index=c.terminal_residuals(data,terminal_log);c.require(residual==item['terminal_scaled_residuals'] and index==item['native_terminal_iteration'],'Stored native terminal residual mismatch')
   previous=terminal_indices[-1] if terminal_indices else 0;c.require(1<=index-previous<=100,'Actual inner iteration cap exceeded');terminal_indices.append(index);residuals.append([residual[k] for k in c.EQUATIONS]);gate['checks']['own_terminal_residuals']=all(v<=1e-6 for v in residual.values())
  gate['all_checks_passed']=all(gate['checks'].values());state_results.append(gate);states.append(s);paths.append(dict(index=step,time_s=s['native_time']['flow-time'],source_kind='native_fluent_hdf',case=str(case),data=str(data)))
  c.save_arrays(out/f'actual-native-state-{step:04d}.npz',s)
 c.require(states,'No actual native states');inverse=np.argsort(states[0]['mapping']);centres=states[0]['native_cell_centres'][inverse];areas=states[0]['native_cell_areas'][inverse]
 c.require(np.max(abs(centres-g['centres']))<1e-12 and np.max(abs(areas-g['areas']))<1e-14,'Common source-order geometry differs')
 checks['all_native_state_and_per_step_residual_gates']=bool(checks['exact_pilot_state_coverage'] and all(s['all_checks_passed'] for s in state_results))
 png=remember(run/'native-velocity.png');checks['actual_native_velocity_png_saved']=png.stat().st_size>1000
 times=np.array([s['native_time']['flow-time'] for s in states]);steps=np.array([r['step'] for r in rows],int)
 np.savez_compressed(out/'native-audit-history.npz',time_s=times,step_index=steps,U_source_history_m_s=np.stack([np.column_stack([s['common']['u'],s['common']['v']]) for s in states]),p_source_history_Pa=np.stack([s['common']['p'] for s in states]),
  source_centres_m=centres,source_areas_m2=areas,source_native_cell_index=inverse,dimensionless_divergence_max=np.array([s['metrics']['maximum_dimensionless_divergence'] for s in state_results]),
  dimensionless_divergence_RMS=np.array([s['metrics']['volume_weighted_RMS_dimensionless_divergence'] for s in state_results]),angular_momentum_per_depth=np.array([s['metrics']['angular_momentum_z_per_unit_depth_kg_m_s'] for s in state_results]),
  solved_time_s=times[1:],terminal_scaled_residuals=np.array(residuals),terminal_residual_equations=np.array(c.EQUATIONS),terminal_native_iteration=np.array(terminal_indices))
 assessment=dict(available=False,final_field_numerically_steady=False)
 if len(states)==401:
  final=np.column_stack([states[400]['common']['u'],states[400]['common']['v']]);old=np.column_stack([states[360]['common']['u'],states[360]['common']['v']]);den=float(np.sum(areas[:,None]*final*final));L=state_results[400]['metrics']['angular_momentum_z_per_unit_depth_kg_m_s'];oldL=state_results[360]['metrics']['angular_momentum_z_per_unit_depth_kg_m_s'];c.require(den>0 and L>0,'Invalid final steady-assessment denominator');change=float(np.sqrt(np.sum(areas[:,None]*(final-old)**2)/den));Lchange=float(abs(L-oldL)/abs(L));assessment=dict(available=True,velocity360to400_relative_volume_L2_final_denominator=change,angular_momentum360to400_relative_change=Lchange,threshold=1e-3,final_field_numerically_steady=bool(change<=1e-3 and Lchange<=1e-3),physical_validation=False)
 checks['source_bytes_unchanged_during_read']=all(c.sha(path)==digest for path,digest in before.items());all_primary=all(checks.values())
 record=dict(kind='Annular Couette actual Fluent400step production audit',checks={k:bool(v) for k,v in checks.items()},all_primary_checks_passed=bool(all_primary),state_results=state_results,native_state_paths=paths,
  steady_assessment=assessment,native_polygon_volume_definition=states[0]['volume_definition'],pressure_convention='Native gauge Pa; sourcecell0 centroid reference; cross-code visual comparison subtracts each volume-weighted mean',
  audited_fields_npz='native-audit-history.npz',audited_fields_npz_sha256=c.sha(out/'native-audit-history.npz'),source_native_cell_index_convention='source_to_native; inverse of native mapping',native_velocity_png=str(png),torque_reported=False,torque_reason='Native torque signs/units/depth integration not independently proved; no torque balance claim',physical_validation=False,mesh_independence=False,time_step_independence=False,production_authorized_by_separate_root_record=True)
 c.write_json(out/'checks.json',record)
 outputs={p.name:c.sha(p) for p in out.iterdir() if p.is_file()}
 c.write_json(out/'source-manifest.json',dict(analyser_sha256=c.sha(__file__),common_sha256=c.sha(c.__file__),inputs=[dict(path=p,sha256=h) for p,h in before.items()],outputs=outputs,all_primary_checks_passed=bool(all_primary),physical_validation=False))
 return record
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=analyse(a.run,a.out);print(json.dumps(r,indent=2,allow_nan=False));return 0 if r['all_primary_checks_passed'] else 1
if __name__=='__main__':raise SystemExit(main())
