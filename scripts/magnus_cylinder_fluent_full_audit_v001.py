"""Read-only streaming Magnus Fluent FULL501state audit; never launches a solver."""
from pathlib import Path
import argparse,ctypes,hashlib,importlib.util,json,sys

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()
def available_memory():
 class MS(ctypes.Structure):_fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ['total','available','total_page','available_page','total_virtual','available_virtual','available_extended']]
 m=MS();m.length=ctypes.sizeof(m)
 if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):raise ctypes.WinError()
 return m.available
def checked_common(path,digest):
 if sha(path)!=digest:raise ValueError('Common reader source changed BEFORE import')
 spec=importlib.util.spec_from_file_location('magnus_audit_common_checked',path);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c);return c
EXPECTED_PRODUCTION_SHA='af8bf603915005174fa390ce091dc8b3f61ca6777e25b75b88bf7eae0dfdbb65'
EXPECTED_COMMON_SHA='1a61e08a2f3624e6ccd36b3c5f244c67d3cbd89eee4865c2a67d030378905820'
def analyse(run,out,source_sha256,common_sha256):
 if sha(__file__)!=source_sha256:raise ValueError('Audit source differs from root-reviewed SHA')
 if common_sha256!=EXPECTED_COMMON_SHA:raise ValueError('Unreviewed Magnus common source')
 c=checked_common(Path(__file__).with_name('magnus_cylinder_fluent_common_v001.py'),common_sha256)
 import numpy as np
 run=Path(run).resolve();out=Path(out).resolve();c.require(run.parent==c.RUN.resolve() and run.is_dir(),'Wrong actual model run');c.require(not out.exists() and out.drive.upper()=='E:' and out.is_relative_to(c.LOCAL),'Fresh model-private E output required');ram=available_memory();c.require(ram>=int(2.5*1024**3),'Read-only audit requires2.5GiB RAM')
 before={}
 def remember(path,expected=None):
  path=Path(path).resolve();c.require(path.is_file(),'Missing bound source '+str(path));key=str(path)
  if key not in before:before[key]=c.sha(path)
  c.require(expected is None or before[key]==expected,'Changed bound source '+str(path));return path
 def load(name):return c.read(remember(run/name))
 source=load('code-provenance.json');auth=load('authorization.json');status=load('status.json');rows=load('states.json');setup=load('mesh-setup.json')
 c.require(auth['approved'] is True and auth['root']=='root' and Path(auth['out']).resolve()==run and auth['solver_processes']==2 and auth['dt_s']==.02 and auth['tandem_pair_and_local_record_complete'] is True,'Wrong root authorization/ranks/time')
 phase=auth['phase'];c.require(phase=='production-495-additional-steps','Unknown full-audit phase');endpoint=500
 expected_action='magnus_cylinder_fluent_production_v001';c.require(auth['action']==expected_action,'Wrong model source action')
 c.require(auth['common_sha256']==common_sha256,'Audited reader differs from actually executed reader');remember(run/'executed-common.py',common_sha256);remember(run/'executed-production.py',EXPECTED_PRODUCTION_SHA);c.require(auth['source_sha256']==EXPECTED_PRODUCTION_SHA,'Actual production source differs from reviewed continuation');remember(run/'TARGETS.md',c.TARGETS_SHA);remember(run/'EXECUTION-PLAN-v001.md',c.PLAN_SHA);remember(run/'mesh-setup.json',c.SETUP_SHA);remember(run/'authorization.json',source['root_authorization_sha256'])
 for r in source['inputs']:remember(r['path'],r['sha256'])
 producer=c.checked_module(remember(Path(__file__).with_name('magnus_cylinder_fluent_production_v001.py'),EXPECTED_PRODUCTION_SHA),EXPECTED_PRODUCTION_SHA);producer.validate_budget(auth);pilot_parent,accepted_rows,accepted_inputs=producer.verify_accepted_pilot(c,auth)
 for path,digest in accepted_inputs.items():remember(path,digest)
 c.require(rows[:6]==accepted_rows,'Authoritative original pilot0..5 ledger/settings changed')
 remember(__file__,source_sha256);remember(c.__file__,common_sha256);remember(c.MESH,c.MESH_SHA);remember(c.MESH.parent/'geometry.npz',c.GEOMETRY_SHA);remember(c.MESH.parent/'setup.json',c.SETUP_SHA);remember(c.LOCAL/'TARGETS.md',c.TARGETS_SHA);remember(c.LOCAL/'EXECUTION-PLAN-v001.md',c.PLAN_SHA)
 c.require(source['source_sha256']==auth['source_sha256'] and source['processor_count']==2 and source['available_ram_bytes']>=4*1024**3 and source['free_E_bytes']>=10*1024**3,'Missing actual launch source/resource proof')
 checks=dict(source_authorization_resources=True,exact_state_coverage=[r['step'] for r in rows]==list(range(endpoint+1)),driver_completed=status['stage']=='production_completed' and status['actual_total_steps']==500 and status['actual_additional_steps']==495 and status['actual_end_time_s']==10.)
 c.require(rows and [r['step'] for r in rows]==list(range(len(rows))),'Noncontiguous/duplicate native state ledger');quality=load('geometry-quality.json');checks['native_mesh_quality']=quality['checked_before_flow'] is True and quality['minimum_orthogonal_quality']>=.15 and quality['minimum_native_cell_volume']>0 and abs(quality['total_native_cell_volume']-c.POLYGON_AREA)<=quality['total_volume_print_rounding_tolerance_m3']
 for name in ['loaded-settings.json','initialized-settings.json']:c.settings_gate(load(name))
 checks['final_settings_saved']=(run/'final-settings.json').is_file()
 if checks['final_settings_saved']:c.settings_gate(load('final-settings.json'))
 resume=load('resume-preservation-proof.json');c.require(resume['actual_resumed_step']==5 and resume['time_s']==.1 and resume['no_initialization_called'] is True and resume['all_primary_fields_native_mesh_and_body_forces_retained'] is True and resume['diagnostic_duplicate_excluded_from501_state_history'] is True,'Native resume preservation proof failed');remember(run/resume['case'],resume['case_sha256']);remember(run/resume['data'],resume['data_sha256']);checks['native_resume_preservation_proof']=all(v<=1e-12 for v in resume['maximum_absolute_field_changes_SI'].values())
 checks['actual_native_physics_methods']=True;g=np.load(c.MESH.parent/'geometry.npz',allow_pickle=False);c.geometry_check(g);out.mkdir(parents=True)
 nstates=len(rows);nsolved=nstates-1;u=np.lib.format.open_memmap(out/'U-source-history.npy',mode='w+',dtype='float64',shape=(nstates,c.N,2));p=np.lib.format.open_memmap(out/'p-source-history.npy',mode='w+',dtype='float64',shape=(nstates,c.N))
 time=np.empty(nstates);steps=np.empty(nstates,int);dmax=np.empty(nstates);drms=np.empty(nstates);imbalance=np.empty(nstates);courant=np.empty(nstates);residuals=np.empty((nsolved,3));iterations=np.empty(nsolved,int);force={term:np.empty((nsolved,1,2)) for term in ['pressure','viscous','total']};paths=[];state_results=[];cache=None;previous=0;endpoint_reports={}
 for i,row in enumerate(rows):
  step=row['step'];case=remember(run/row['case'],row['case_sha256']);data=remember(run/row['data'],row['data_sha256']);c.require(case.is_relative_to(c.RUN) and data.is_relative_to(c.RUN),'Native state outside model')
  if cache is None:cache=c.map_native_mesh(case,g)
  s=c.read_state(case,data,g,cache);gate=c.state_gates(s,step);c.settings_gate(c.read(remember(run/row['settings'],row['settings_sha256'])))
  for key,value in gate['checks'].items():c.require(row['native_gates']['checks'][key] is value,'Saved driver/native replay gate differs '+key)
  if step:
   transcript=remember(case.parent.parent/'solver-transcript.txt');rr,index=c.terminal_residuals(data,transcript);c.require(rr==row['terminal_scaled_residuals'] and index==row['native_terminal_iteration'],'Driver/native residual differs');c.require(1<=index-previous<=100,'Observed native inner cap exceeded');previous=index;residuals[i-1]=[rr[k] for k in c.EQUATIONS];iterations[i-1]=index;gate['checks']['own_terminal_residuals']=bool(all(v<=1e-6 for v in rr.values()));forces=c.body_forces(s);c.require(forces==row['body_forces'],'Driver/native per-body force differs')
   for term in force:force[term][i-1]=[forces[name][term+'_coefficient'] for name in c.BODIES]
   if step==5:
    duplicate=c.read_state(run/resume['case'],run/resume['data'],g,cache);c.require(c.state_gates(duplicate,5)['all_checks_passed'],'Native resume duplicate gates fail');actual_resume=producer.resume_differences(dict(**s['native'],all_native_face_mass_flux=s['flux']),dict(**duplicate['native'],all_native_face_mass_flux=duplicate['flux']));c.require(actual_resume==resume['maximum_absolute_field_changes_SI'] and c.body_forces(duplicate)==forces,'Independent actual resume preservation differs from proof');del duplicate
   if step in [5,500]:
    reports={}
    for body in c.BODIES:
     for label,direction in [('x',[1,0,0]),('y',[0,1,0])]:
      r=row['native_force_reports'][body+'_'+label];path=remember(case.parent.parent/r['path'],r['sha256']);actual=c.force_report_agreement(forces,body,direction,path);c.require(actual==r['report'],'Driver/native ASCII report differs');reports[body+'_'+label]=actual
    endpoint_reports[str(step)]=reports
  gate['all_checks_passed']=bool(all(gate['checks'].values()));state_results.append(gate);time[i]=s['native_time']['flow-time'];steps[i]=step;u[i,:,0]=s['common']['u'];u[i,:,1]=s['common']['v'];p[i]=s['common']['p'];m=gate['metrics'];dmax[i]=m['maximum_dimensionless_divergence'];drms[i]=m['volume_weighted_RMS_dimensionless_divergence'];imbalance[i]=m['relative_mass_imbalance'] if m['relative_mass_imbalance'] is not None else np.nan;courant[i]=m['maximum_matched_cell_Courant'];paths.append(dict(index=step,time_s=time[i],source_kind='native_fluent_hdf',case=str(case),data=str(data)))
  if i%50==0:print(json.dumps(dict(read_only_native_state=i,total=nstates)),flush=True)
 u.flush();p.flush();inverse=np.argsort(cache['mapping']);centres=cache['native_cell_centres'][inverse];areas=cache['native_cell_areas'][inverse];c.require(np.max(abs(centres-g['centres']))<1e-12 and np.max(abs(areas-g['areas']))<1e-14,'Common source-order geometry differs')
 checks['all_native_solved_step_gates']=all(r['all_checks_passed'] for r in state_results);checks['separate_native_force_endpoint_reports']=str(endpoint) in endpoint_reports and '5' in endpoint_reports
 fields=out/'native-audit-history.npz';np.savez_compressed(fields,time_s=time,step_index=steps,U_source_history_m_s=u,p_source_history_Pa=p,source_centres_m=centres,source_areas_m2=areas,source_native_cell_index=inverse,solved_time_s=time[1:],dimensionless_divergence_max=dmax,dimensionless_divergence_RMS=drms,relative_boundary_imbalance=imbalance,native_flux_Courant_max=courant,body_pressure_force_coefficient_history=force['pressure'],body_viscous_force_coefficient_history=force['viscous'],body_total_force_coefficient_history=force['total'],**c.final_surface_pressure(s),body_names=np.array(c.BODIES),force_component_names=np.array(['x','y']),terminal_scaled_residuals=residuals,terminal_residual_equations=np.array(c.EQUATIONS),terminal_native_iteration=iterations)
 png=run/'native-velocity.png';png_saved=png.is_file()
 if png_saved:remember(png)
 checks['source_bytes_unchanged_during_read']=all(c.sha(path)==digest for path,digest in before.items());all_primary=bool(all(checks.values()));record=dict(kind='Actual native Fluent rotating-cylinder streaming full501state audit',phase=phase,checks={k:bool(v) for k,v in checks.items()},all_primary_checks_passed=all_primary,state_results=state_results,native_state_paths=paths,force_endpoint_report_proof=endpoint_reports,actual_last_saved_step=int(steps[-1]),actual_last_saved_time_s=float(time[-1]),surface_pressure_is_actual_last_saved_state=True,initial_forces_available=False,initial_force_reason='Uncomputed initialization; solved-only force arrays have no invented initial row',audited_fields_npz=fields.name,audited_fields_npz_sha256=c.sha(fields),source_native_cell_index_convention='source_to_native; inverse native mapping',native_velocity_png=str(png) if png_saved else None,native_velocity_png_available=png_saved,pressure_convention='Native gauge Pa; operating0Pa and fixed-pressure outlet; no imposed closed-domain reference-point constraint; fixed outlet supplies pressure gauge',native_polygon_volume_definition=cache['volume_definition'],vorticity_available=False,signed_vorticity_used=False,force_normalization_N_per_body=c.NORMALIZATION,wall_velocity_target='Omega_z=-2rad/s cross actual native facet centre, projected tangent; checked only solved steps1+',body_surface_native_face_ids_convention='Native Fluent1-based face ids; source ids also1-based, sorted source order; actual saved face SV_P, never cell interpolation',native_flux_Courant_definition='dt=.02s *sum(abs(nativeface massflux))/(2*rho*actualnativecellarea*depth1); descriptive',analysis_memory_policy='One native state and fixed mapped mesh; E-backed U/p histories; no all-state native-object list',available_ram_before_read_bytes=ram,physical_validation=False,mesh_independence=False,time_step_independence=False,periodicity_established=False)
 c.write_json(out/'checks.json',record);outputs={q.name:c.sha(q) for q in out.iterdir() if q.is_file()};c.write_json(out/'source-manifest.json',dict(analyser_sha256=source_sha256,common_sha256=common_sha256,inputs=[dict(path=p,sha256=h) for p,h in before.items()],outputs=outputs,all_primary_checks_passed=all_primary,physical_validation=False));return record

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--source-sha256',required=True);p.add_argument('--common-sha256',required=True);a=p.parse_args();r=analyse(a.run,a.out,a.source_sha256,a.common_sha256);print(json.dumps(dict(all_primary_checks_passed=r['all_primary_checks_passed'],native_states=len(r['native_state_paths']))));return 0 if r['all_primary_checks_passed'] else 1
if __name__=='__main__':raise SystemExit(main())
