"""Read-only native Couette/tandem reader fixtures; no solver, WSL or new case."""
from pathlib import Path
import argparse,copy,importlib.util,json,sys
import numpy as np
ROOT=Path('E:/CFD-Learning-Lab');LOCAL=ROOT/'local/magnus-cylinder-v001'

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args();out=a.out.resolve()
 if out.exists() or not out.is_relative_to(LOCAL):raise ValueError('Fresh Magnus-private E output required')
 source=Path(__file__).with_name('magnus_cylinder_fluent_common_v001.py');spec=importlib.util.spec_from_file_location('magnus_fixture_common',source);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
 before={};results=[]
 def remember(path):
  path=Path(path);before[str(path)]=c.sha(path);return path
 def geometry(path):
  with np.load(remember(path),allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
 def passed(name,detail):results.append(dict(name=name,passed=True,evidence=detail))
 def rejected(name,fn):
  try:fn()
  except (ValueError,KeyError,IndexError) as e:passed(name,dict(rejected=True,reason=str(e)));return
  raise AssertionError('Negative fixture accepted: '+name)
 remember(__file__);remember(source)
 cg=geometry(ROOT/'native/runs/2026-11-02-annular-couette-v001/mesh/geometry.npz');cr=ROOT/'native/runs/2026-11-02-annular-couette-v001/fluent-pilot-v002'
 cc=remember(cr/'states/state-0005.cas.h5');cd=remember(cr/'states/state-0005.dat.h5');cm=c.map_native_mesh(cc,cg,counts=dict(interior=3008,inner=64,outer=64),types=dict(interior=2,inner=3,outer=3));cs=c.read_fields(cd,cm)
 assert len(np.unique(cm['mapping']))==1536 and len(np.unique(cm['face_mapping']))==3136
 passed('actual_Couette_all_polygon_cells_faces_and_zones',dict(cells=1536,faces=3136,mesh_fingerprint=cm['mesh_fingerprint']))
 rm=c.moving_wall_metrics(cs,'inner',1.)
 assert rm['maximum_velocity_error_m_s']<=1e-10 and rm['maximum_normal_velocity_m_s']<=1e-10
 passed('actual_Couette_positive_rotation_native_face_UV_and_tangency',rm)
 wrong=c.moving_wall_metrics(cs,'inner',-1.)
 assert wrong['maximum_velocity_error_m_s']>.09
 passed('opposite_rotation_sign_is_detectable',wrong)
 rr=c.read(remember(cr/'native-physical-readback.json'))['rotation'];st=c.read(remember(cr/'loaded-settings.json'));c.rotation_settings_gate(rr,st['walls']['inner']['momentum'],1.)
 passed('actual_Couette_live_angular_units_and_settings',rr)
 surface=c.mapped_surface_pressure(cs,'inner',np.arange(3009,3073));assert surface['body_surface_pressure_final_Pa'].shape==(64,) and np.isfinite(surface['body_surface_pressure_final_Pa']).all()
 passed('actual_Couette_boundary_pressure_and_source_native_face_mapping',dict(faces=64,source_face_ids=surface['body_surface_source_face_ids'].tolist(),native_face_ids=surface['body_surface_native_face_ids'].tolist(),pressure_source='Actual saved face SV_P, not cell interpolation'))
 rejected('wrong_expected_surface_face_identity_rejected',lambda:c.mapped_surface_pressure(cs,'inner',np.arange(3010,3074)))
 badstate=dict(cs);badstate['pressure_coverage']=cs['pressure_coverage'].copy();badstate['pressure_coverage'][cs['zones']['inner'][0]]=False
 rejected('missing_actual_surface_pressure_rejected',lambda:c.mapped_surface_pressure(badstate,'inner',np.arange(3009,3073)))

 rejected('wrong_clockwise_setting_rejected',lambda:c.rotation_settings_gate(rr,st['walls']['inner']['momentum'],-2.))
 bad=copy.deepcopy(rr);bad['rotation_speed']['native_units']='rpm';rejected('ambiguous_native_angular_units_rejected',lambda:c.rotation_settings_gate(bad,st['walls']['inner']['momentum'],1.))
 tg=geometry(ROOT/'native/runs/2026-11-02-tandem-offset-wake-v001/mesh/geometry.npz');tr=ROOT/'native/runs/2026-11-02-tandem-offset-wake-v001/fluent-pilot-v001'
 tc=remember(tr/'states/state-0005.cas.h5');td=remember(tr/'states/state-0005.dat.h5');counts=dict(interior=18692,inlet=60,outlet=60,top=160,bottom=160,body1=32,body2=32);types=dict(interior=2,inlet=10,outlet=5,top=7,bottom=7,body1=3,body2=3)
 tm=c.map_native_mesh(tc,tg,counts,types);ts=c.read_fields(td,tm)
 assert len(np.unique(tm['mapping']))==9472 and len(np.unique(tm['face_mapping']))==19196 and np.isfinite(ts['flux']).all()
 passed('actual_tandem_all_cells_faces_native_phi_and_open_zones',dict(cells=9472,faces=19196,mesh_fingerprint=tm['mesh_fingerprint']))
 old_rows=c.read(remember(tr/'states.json'));old_co=old_rows[-1]['native_gates']['metrics']['maximum_matched_cell_Courant'];actual_co=c.flux_metrics(ts)['maximum_matched_cell_Courant'];assert abs(actual_co-old_co*.02/.05)<1e-13
 passed('actual_tandem_native_phi_Courant_scales_with_frozen_Magnus_dt',dict(old_dt_s=.05,new_dt_s=.02,old_native_Courant=old_co,new_same_field_Courant=actual_co,no_new_flow_computed=True))
 pilot_source=remember(Path(__file__).with_name('magnus_cylinder_fluent_pilot_v001.py'));pspec=importlib.util.spec_from_file_location('magnus_fixture_pilot',pilot_source);pilot=importlib.util.module_from_spec(pspec);pspec.loader.exec_module(pilot)
 from types import SimpleNamespace
 for native_run,area,label in [(cr,float(cg['areas'].sum()),'Couette'),(tr,float(tg['areas'].sum()),'tandem')]:
  transcript=remember(native_run/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace');q=pilot.quality(transcript,SimpleNamespace(require=c.require,POLYGON_AREA=area));assert q['minimum_orthogonal_quality']>=.15 and q['exact_native_saved_polygon_mapping_required_separately'] is True
  passed('actual_'+label+'_rounded_native_quality_print_parser',q)
  rejected('wrong_'+label+'_polygon_volume_rejected',lambda text=transcript:pilot.quality(text,SimpleNamespace(require=c.require,POLYGON_AREA=area*1.01)))

 forces={}
 for body in ['body1','body2']:
  lo,hi=ts['zones'][body];forces[body]=c.integrate_body(ts['face_pressure'][lo:hi],ts['wall_shear'][lo:hi],ts['face_area_vector'][lo:hi])
  for axis,direction in [('x',[1,0,0]),('y',[0,1,0])]:
   report=remember(tr/'force-reports'/f'step0005-{body}-{axis}.txt');proof=c.force_report_agreement(forces,body,direction,report);passed('actual_native_HDF_ASCII_force_'+body+'_'+axis,proof)
 lo,hi=ts['zones']['body1'];badforce=copy.deepcopy(forces);badforce['body1']=c.integrate_body(ts['face_pressure'][lo:hi],ts['wall_shear'][lo:hi]*np.linalg.norm(ts['face_area_vector'][lo:hi],axis=1)[:,None],ts['face_area_vector'][lo:hi])
 rejected('extra_area_multiplier_on_integrated_shear_rejected_by_native_report',lambda:c.force_report_agreement(badforce,'body1',[1,0,0],tr/'force-reports/step0005-body1-x.txt'))
 corrupt={k:v.copy() for k,v in tg.items()};corrupt['face_zone'][np.flatnonzero(corrupt['face_zone']==8)[0]]=9
 rejected('wrong_source_body_zone_rejected',lambda:c.map_native_mesh(tc,corrupt,counts,types))
 corrupt={k:v.copy() for k,v in cg.items()};corrupt['quads'][1]=corrupt['quads'][0]
 rejected('nonbijective_source_cell_geometry_rejected',lambda:c.map_native_mesh(cc,corrupt,dict(interior=3008,inner=64,outer=64),dict(interior=2,inner=3,outer=3)))
 mg=geometry(c.MESH.parent/'geometry.npz');c.geometry_check(mg);passed('selected_Magnus_geometry_matches_frozen_scope',dict(cells=len(mg['quads']),faces=len(mg['face_nodes']),polygon_area_m2=float(mg['areas'].sum()),no_native_Magnus_fields_used=True))
 badg={k:v.copy() for k,v in mg.items()};badg['nodes'][1]=badg['nodes'][0];rejected('duplicate_Magnus_seam_vertex_rejected',lambda:c.geometry_check(badg))
 assert all(c.sha(path)==digest for path,digest in before.items());passed('all_fixture_source_and_native_bytes_preserved',dict(files=len(before)))
 out.mkdir();record=dict(kind='Actual old native reader compatibility and negative fixtures; not a Magnus CFD result',passed=True,solver_launched=False,wsl_started=False,new_native_case_created=False,tests=results,inputs=[dict(path=k,sha256=v) for k,v in before.items()]);c.write_json(out/'checks.json',record);print(json.dumps(dict(passed=True,tests=len(results),out=str(out))));return 0
if __name__=='__main__':raise SystemExit(main())