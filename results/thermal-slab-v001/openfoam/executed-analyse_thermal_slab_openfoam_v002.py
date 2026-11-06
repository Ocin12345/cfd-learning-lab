"""Read-only source-bound audit of the actual Euler transient solid slab.

No CFD/WSL is started. Writes only a fresh external E: audit folder.
"""
from pathlib import Path
import argparse,csv,ctypes,hashlib,json,os,re,sys
import numpy as np
from prepare_openfoam_thermal_slab_v002 import (
 PATCHES,EMPTY,FLOAT,CELLS,INTERNAL_FACES,TOTAL_FACES,POINTS,REFERENCE_SHA,INSTALLED_SOURCES,
 require,sha,clean,dictionaries,native_poly,native_patches,native_geometry,mesh_quality,json_new,copied,checked_case,
)
from thermal_slab_reference_v001 import cell_average_temperature
class Evidence:
 def __init__(self):self.before={}
 def touch(self,path):
  path=Path(path).resolve();require(path.is_file(),'Missing actual native evidence '+str(path))
  if path not in self.before:self.before[path]=sha(path)
  return path
 def text(self,path):return self.touch(path).read_text(encoding='utf-8-sig')
 def read_json(self,path):return json.loads(self.text(path))
 def finish(self):
  rows=[dict(path=str(p),bytes=p.stat().st_size,sha256_before=s,sha256_after=sha(p)) for p,s in self.before.items()]
  require(all(r['sha256_before']==r['sha256_after'] for r in rows),'Actual input changed during read-only audit')
  return dict(inputs_unchanged_during_analysis=True,sources=rows,solver_started_or_modified=False)
def memory_floor():
 if os.name=='nt':
  class Status(ctypes.Structure):
   _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ('totalPhysical','availablePhysical','totalPage','availablePage','totalVirtual','availableVirtual','availableExtended')]
  s=Status();s.length=ctypes.sizeof(s);require(bool(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(s))),'Cannot inspect available RAM')
  require(s.availablePhysical>=int(2.5*1024**3),'Postprocess requires2.5GiB available RAM')
def finite_tokens(text):
 require(not re.search(r'(?i)(?<![A-Za-z0-9_])[-+]?(?:nan|inf(?:inity)?)(?![A-Za-z0-9_])',text),'Nonfinite native field token')
def value(text,count,components,label='internalField'):
 kind='vector' if components==3 else 'scalar'
 m=re.search(r'\b'+label+r'\s+uniform\s+([^;]+);',text)
 if m:
  row=np.fromstring(m[1].replace('(',' ').replace(')',' '),sep=' ')
  require(row.size==components and np.isfinite(row).all(),'Invalid native uniform '+label)
  return np.tile(row,(count,1)) if components>1 else np.full(count,row[0])
 m=re.search(r'\b'+label+r'\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
 require(m is not None and int(m[1])==count,'Native '+label+' count differs')
 data=np.fromstring(m[2].replace('(',' ').replace(')',' '),sep=' ')
 require(data.size==count*components and np.isfinite(data).all(),'Native '+label+' length/finite gate')
 return data.reshape(count,components) if components>1 else data
def field(raw,count,components,dimensions):
 text=clean(raw);finite_tokens(text);require(re.search(r'\bformat\s+ascii\s*;',text),'ASCII native field required')
 d=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',text)
 require(d is not None and list(map(int,d[1].split()))==dimensions,'Native dimensions differ')
 return value(text,count,components)
def boundaries(raw):
 text=clean(raw);finite_tokens(text);require('boundaryField' in text,'Missing actual native boundaries')
 rows=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',text.split('boundaryField',1)[1]))
 require(set(rows)==set(PATCHES)|{EMPTY},'Actual native patch set differs')
 return rows
def temperature(raw,patches):
 result=field(raw,CELLS,1,[0,0,0,1,0,0,0]);rows=boundaries(raw)
 for name,expected in (('hot',400),('cold',300)):
  require(re.search(r'\btype\s+fixedValue\s*;',rows[name]),'Native isothermal patch type differs')
  require(np.max(np.abs(value(rows[name],patches[name]['count'],1,'value')-expected))<=1e-12,'Actual native wall temperature differs')
 for name in ('insulated_bottom','insulated_top'):
  require(re.search(r'\btype\s+zeroGradient\s*;',rows[name]),'Native insulation temperature type differs')
 require(re.search(r'\btype\s+empty\s*;',rows[EMPTY]),'Native temperature extrusion planes not empty')
 return result
def mapping(centres,points,volumes,geometry_file):
 require(centres.shape==(CELLS,3) and np.max(np.abs(centres[:,2]))<=1e-10,'Native C count/extrusion centre differs')
 ix=np.rint(centres[:,0]/.001-.5).astype(int);iy=np.rint(centres[:,1]/.001-.5).astype(int)
 require(np.all((ix>=0)&(ix<100)&(iy>=0)&(iy<20)),'Native C outside slab')
 expected=np.column_stack(((ix+.5)*.001,(iy+.5)*.001))
 require(np.max(np.abs(centres[:,:2]-expected))<=1e-10,'Native C off canonical slab grid')
 keys=iy*100+ix;require(np.array_equal(np.sort(keys),np.arange(CELLS)),'Native C grid bijection failed')
 with np.load(geometry_file,allow_pickle=False) as g:
  nodes=np.asarray(g['nodes']);quads=np.asarray(g['quads']);source_c=np.asarray(g['centres']);areas=np.asarray(g['areas'])
  sx=np.asarray(g['ix']).astype(int);sy=np.asarray(g['iy']).astype(int)
  require(nodes.shape==(2121,2) and quads.shape==(CELLS,4) and source_c.shape==(CELLS,2),'Canonical source shapes differ')
  require(np.array_equal(g['grid_shape'],[20,100]),'Canonical source grid shape differs')
 skeys=sy*100+sx;require(np.array_equal(np.sort(skeys),np.arange(CELLS)),'Canonical source cells not bijective')
 source_by_key=np.argsort(skeys);native_to_source=source_by_key[keys]
 centre_error=float(np.max(np.abs(centres[:,:2]-source_c[native_to_source])))
 area_error=float(np.max(np.abs(volumes/.01-areas[native_to_source])))
 require(centre_error<=1e-10 and area_error<=1e-12,'Actual native/source centre or area differs')
 points=np.asarray(points);require(points.shape==(POINTS,3),'Native points count differs')
 pints=np.rint(points[:,:2]/.001).astype(int);nints=np.rint(nodes/.001).astype(int)
 require(np.max(np.abs(points[:,:2]-.001*pints))<=1e-10,'Native point off canonical lattice')
 require(np.all((pints[:,0]>=0)&(pints[:,0]<=100)&(pints[:,1]>=0)&(pints[:,1]<=20)),'Native vertex bounds differ')
 nkeys=nints[:,1]*101+nints[:,0];pkeys=pints[:,1]*101+pints[:,0]
 lookup={int(k):i for i,k in enumerate(nkeys)};require(len(lookup)==2121 and all(int(k) in lookup for k in pkeys),'Projected native vertex absent from source')
 point_map=np.array([lookup[int(k)] for k in pkeys],dtype=np.int64)
 low=np.abs(points[:,2]+.005)<=1e-10;high=np.abs(points[:,2]-.005)<=1e-10
 require(np.all(low^high) and low.sum()==high.sum()==2121,'Native points not exactly two extrusion planes')
 require(np.array_equal(np.sort(point_map[low]),np.arange(2121)) and np.array_equal(np.sort(point_map[high]),np.arange(2121)),'Projected native vertex bijection failed')
 vertex_error=float(np.max(np.abs(points[:,:2]-nodes[point_map])));require(vertex_error<=1e-10,'Actual native vertex differs from source')
 xleft=np.min(nodes[quads,:,][...,0],axis=1)[native_to_source]
 xright=np.max(nodes[quads,:,][...,0],axis=1)[native_to_source]
 return dict(native_cell_to_source_cell=native_to_source,native_cell_index_grid=np.argsort(keys).reshape(20,100),source_cell_index_grid=source_by_key.reshape(20,100),
  ix=ix,iy=iy,native_points=points,native_point_to_source_node=point_map,x_left_m=xleft,x_right_m=xright),dict(native_cells=CELLS,source_cells=CELLS,
  centre_max_error_m=centre_error,area_max_error_m2=area_error,vertex_max_error_m=vertex_error,source_vertices_per_plane=2121,native_vertices=POINTS,
  native_to_source_cells_bijective=True,native_to_grid_bijective=True,projected_vertex_bijection_each_plane=True,ordering_assumed=False),quads
def verify_cell_vertices(geometry,m,quads):
 vertices=[set() for _ in range(CELLS)]
 for i,face in enumerate(geometry['faces']):
  vertices[geometry['owner'][i]].update(face)
  if i<INTERNAL_FACES:vertices[geometry['neighbour'][i]].update(face)
 for cell,labels in enumerate(vertices):
  require(len(labels)==8,'Native thermal cell does not have8vertices')
  expected=set(quads[m['native_cell_to_source_cell'][cell]].tolist())
  for sign in (-1,1):
   projected={int(m['native_point_to_source_node'][v]) for v in labels if sign*m['native_points'][v,2]>0}
   require(projected==expected,'Actual native cell differs from canonical source quad')
def face_areas(geometry):
 pts=np.asarray(geometry['points']);result=[]
 for face in geometry['faces']:
  p=pts[face];result.append(.5*float(np.linalg.norm(np.cross(p,np.roll(p,-1,axis=0)).sum(axis=0))))
 require(all(np.isfinite(result)) and min(result)>0,'Invalid actual native face areas')
 return np.asarray(result)
def heat_field(raw,patches,areas):
 internal=field(raw,CELLS,1,[1,0,-3,0,0,0,0]);require(np.max(np.abs(internal))<=1e-14,'wallHeatFlux internal storage must be zero')
 rows=boundaries(raw);Q=[];stats=[]
 for name in PATCHES:
  patch=patches[name];q=value(rows[name],patch['count'],1,'value')
  a=areas[patch['start']:patch['start']+patch['count']]
  power=float(np.dot(q,a));Q.append(power);stats.append([float(q.min()),float(q.max()),power,power/float(a.sum())])
 require(re.search(r'\btype\s+empty\s*;',rows[EMPTY]),'Heat-flux field planes not empty')
 return np.asarray(Q),np.asarray(stats)
def heat_table(text):
 require(re.search(r'^#.*Time.*patch.*Q\s*\[W\]',text,re.M),'Native wallHeatFlux table header lacks exact Time/patch/QW semantics')
 rows={};duplicates=0
 for line in text.splitlines():
  if not line.strip() or line.lstrip().startswith('#'):continue
  tokens=line.split();require(len(tokens)==6 and tokens[1] in PATCHES,'Native heat table row/patch differs')
  time=float(tokens[0]);values=np.array(list(map(float,tokens[2:])))
  require(time.is_integer() and int(time) in range(0,201,2) and np.isfinite(values).all(),'Invalid native heat table time/value')
  key=(int(time),tokens[1])
  if key in rows:
   require(np.allclose(rows[key],values,rtol=1e-12,atol=1e-10),'Repeated heat table row changed at same actual time')
   duplicates+=1
  rows[key]=values
 require(all((t,n) in rows for t in range(2,201,2) for n in PATCHES),'Native heat table missing solved time/patch')
 return rows,duplicates
def residual_history(text):
 require(not re.search(r'FOAM FATAL|Segmentation fault|Floating point exception\s*$',text,re.M),'Actual thermal solver fatal')
 require(re.search(r'^End\s*$',text,re.M),'Actual thermal solver log missing End')
 sections=re.split(r'^Time = ('+FLOAT+r')s?\s*$',text,flags=re.M);rows=[]
 for i in range(1,len(sections),2):
  time=float(sections[i]);samples=re.findall(r'Solving for e, Initial residual = ('+FLOAT+r'), Final residual = ('+FLOAT+r'), No Iterations (\d+)',sections[i+1])
  require(len(samples)==1,'Require exactly one actual energy solve per physical step')
  initial,final,iterations=samples[0]
  require(time.is_integer() and np.isfinite(float(initial)) and np.isfinite(float(final)),'Invalid thermal residual')
  rows.append(dict(time_s=int(time),initial_energy_linear_residual=float(initial),final_energy_linear_residual=float(final),linear_iterations=int(iterations)))
 require([r['time_s'] for r in rows]==list(range(2,201,2)),'Require actual physical times2..200 by2')
 return rows
def time_metadata(raw,time,index):
 text=clean(raw);result={}
 for name in ('beginTime','value','index','deltaT','deltaT0'):
  m=re.search(r'\b'+name+r'\s+('+FLOAT+r')\s*;',text);require(m is not None,'Missing native time metadata '+name);result[name]=float(m[1])
 require(result==dict(beginTime=0,value=time,index=index,deltaT=2,deltaT0=2),'Actual native time/dt/index differs from Euler schedule')
 return result
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 case=checked_case(a.run);out=a.out.resolve();require(not out.exists() and out.drive.upper()=='E:' and case not in out.parents,'Fresh external E: audit output required')
 memory_floor();e=Evidence();checks={}
 source=Path(__file__);prep=source.with_name('prepare_openfoam_thermal_slab_v002.py');reference=source.with_name('thermal_slab_reference_v001.py')
 for p in (source,prep,reference):e.touch(p)
 require(sha(reference)==REFERENCE_SHA,'Frozen analytical reference changed')
 declaration=e.read_json(case/'declared-comparison.json');status=e.read_json(case/'run-status.json')
 require(status.get('stage')=='finished_execution' and status.get('final_time_s')==200 and status.get('time_steps')==100,'Actual100-step thermal execution incomplete')
 for name,expected in dictionaries().items():require(e.text(case/name)==expected,'Actual native dictionary changed '+name)
 checks['exact_Euler_solid_material_boundary_and_solver_dictionaries']=True
 manifest=e.read_json(case/'provenance/prepared-input-sha256.json')
 for name,digest in manifest.items():e.touch(case/name);require(sha(case/name)==digest,'Prepared source changed '+name)
 require(sha(e.touch(case/'input-slab.msh'))==declaration['mesh_source_sha256'],'Shared mesh hash differs')
 require(sha(e.touch(case/'provenance/geometry.npz'))==declaration['geometry_sha256'],'Shared source geometry hash differs')
 require(sha(e.touch(case/'provenance/thermal_slab_reference_v001.py'))==REFERENCE_SHA,'Native reference copy hash differs')
 checks['source_mesh_geometry_reference_and_prepared_hashes_match']=True
 pre=e.read_json(case/'provenance/preflight-validation.json');require(pre.get('passed') is True,'Native preflight not passed')
 poly=native_poly(case)
 for name,digest in pre['native_mesh_sha256'].items():e.touch(poly/name);require(sha(poly/name)==digest,'Native static mesh changed')
 for name,digest in pre['native_initial_C_Vc_sha256'].items():e.touch(case/'0'/name);require(sha(case/'0'/name)==digest,'Actual initial C/Vc changed')
 geometry=native_geometry(poly);patches=native_patches(poly);quality=mesh_quality(e.text(case/'mesh-check.log'))
 checks['full_checkMesh_exact_native_topology_positive_oriented_volumes']=True
 centres0=field(e.text(case/'0/C'),CELLS,3,[0,1,0,0,0,0,0]);centres=field(e.text(case/'200/C'),CELLS,3,[0,1,0,0,0,0,0])
 volumes0=field(e.text(case/'0/Vc'),CELLS,1,[0,3,0,0,0,0,0]);volumes=field(e.text(case/'200/Vc'),CELLS,1,[0,3,0,0,0,0,0])
 require(np.all(volumes>0) and np.max(np.abs(volumes-np.asarray(geometry['volumes'])))<=1e-16 and np.array_equal(volumes,volumes0),'Actual emitted Vc differs from native oriented/static volumes')
 require(np.array_equal(centres,centres0),'Actual initial/final native centres differ on static mesh')
 checks['actual_initial_final_C_Vc_match_static_oriented_geometry']=True
 m,mapping_report,quads=mapping(centres,geometry['points'],volumes,case/'provenance/geometry.npz');verify_cell_vertices(geometry,m,quads)
 checks['actual_cells_grid_vertices_and_source_quads_bijective']=True
 areas=face_areas(geometry)
 residuals=residual_history(e.text(case/'solver.log'));initial_residual=np.array([r['initial_energy_linear_residual'] for r in residuals]);final_residual=np.array([r['final_energy_linear_residual'] for r in residuals])
 checks['every_step_final_linear_energy_residual_le_1e_10']=bool(np.all(final_residual<=1e-10))
 logs=sorted((case/'postProcessing/wallHeatRates').glob('*/wallHeatFlux.dat'));require(len(logs)==1,'Fresh thermal run requires one native heat-rate table')
 table,duplicates=heat_table(e.text(logs[0]))
 times=np.arange(0,201,2);T=[];analytical=[];series_delta=[];Q=[];tableQ=[];cross=[];spread=[];metadata=[];heat_stats=[]
 for n,t in enumerate(times.tolist()):
  memory_floor();temp=temperature(e.text(case/str(t)/'T'),patches);T.append(temp)
  exact=cell_average_temperature(m['x_left_m'],m['x_right_m'],t,80);exact160=cell_average_temperature(m['x_left_m'],m['x_right_m'],t,160)
  analytical.append(exact);series_delta.append(float(np.max(np.abs(exact-exact160))))
  grid=temp[m['native_cell_index_grid']];spread.append(float(np.max(np.ptp(grid,axis=0))))
  if t>0:
   metadata.append(time_metadata(e.text(case/str(t)/'uniform/time'),t,n))
   power,stats=heat_field(e.text(case/str(t)/'wallHeatFlux'),patches,areas);reported=np.asarray([table[(t,name)] for name in PATCHES])
   Q.append(power);tableQ.append(reported[:,2]);heat_stats.append(stats)
   # Table streams may have fewer printed digits than native ASCII fields.
   # Scale by Qref as declared before Run; never multiply native table W by area.
   scale=np.maximum(np.maximum(np.abs(power),np.abs(reported[:,2])),2.)
   cross.append(np.abs(power-reported[:,2])/scale)
 T=np.asarray(T);analytical=np.asarray(analytical);Q=np.asarray(Q);tableQ=np.asarray(tableQ);cross=np.asarray(cross);spread=np.asarray(spread)
 require(T.shape==(101,CELLS) and Q.shape==(100,4),'Actual native thermal history shape differs')
 checks['all101_actual_native_T_states_finite_and_correct_counts']=bool(np.isfinite(T).all())
 checks['initial_actual_interior_temperature_exact300K']=bool(np.all(T[0]==300))
 checks['all_native_T_within_300_400_K_tolerance1e_6']=bool(np.all((T>=300-1e-6)&(T<=400+1e-6)))
 checks['all_times_max_y_spread_le_1e_6_K']=bool(np.max(spread)<=1e-6)
 checks['actual_native_times_indices_and_dt2_match_all100_steps']=True
 checks['native_wall_heat_field_table_Q_crosscheck_le_1e_5_scaled']=bool(np.max(cross)<=1e-5)
 checks['actual_hot_heats_and_cold_cools_with_inward_Q_convention']=bool(np.all(Q[:,0]>0) and np.all(Q[:,1]<=1e-12))
 checks['analytic_80_vs160_series_terms_le_1e_9_K']=max(series_delta)<=1e-9
 energy=np.sum(1e6*volumes[None,:]*(T-300),axis=1);storage=np.diff(energy)/2;netQ=Q.sum(axis=1);balance=storage-netQ
 denominator=np.maximum(np.maximum(np.sum(np.abs(Q),axis=1),np.abs(storage)),2.);relative=np.abs(balance)/denominator
 checks['every_Euler_step_global_energy_relative_error_le_1e_5']=bool(np.max(relative)<=1e-5)
 error=T-analytical;max_error=np.max(np.abs(error),axis=1);rms=np.sqrt(np.average(error**2,axis=1,weights=volumes))
 checks['final_t200_cell_average_analytic_max_error_le_0p5K']=bool(max_error[-1]<=.5)
 checks['final_t200_cell_average_analytic_RMS_error_le_0p25K']=bool(rms[-1]<=.25)
 auth=e.read_json(case/'provenance/root-run-authorization.json')
 require(auth.get('approved') is True and auth.get('approved_by')=='root' and auth.get('action')=='run_thermal_slab_openfoam_v002' and auth.get('allowed_mode')=='Run','Missing exact root thermal Run authorization')
 require(sha(e.touch(case/'provenance/run-openfoam-executed.sh'))==sha(e.touch(case/'provenance/run-openfoam-v002.sh'))==auth.get('wrapper_sha256'),'Executed thermal wrapper hash differs')
 require(sha(case/'provenance/prepare_openfoam_thermal_slab_v002.py')==auth.get('preparation_sha256') and sha(case/'provenance/Launch-OpenFOAM-v002.ps1')==auth.get('launcher_sha256'),'Reviewed thermal preparation/launcher hashes differ')
 checks['reviewed_thermal_execution_sources_match_root_authorization']=True
 rows=[line.split(None,1) for line in e.text(case/'provenance/installed-source-sha256.txt').splitlines() if line.strip()]
 require(len(rows)==len(INSTALLED_SOURCES) and {p[len('/opt/openfoam14/'):] for d,p in rows}==set(INSTALLED_SOURCES),'Installed solid source list differs')
 for digest,path in rows:
  require(path.startswith('/opt/openfoam14/'),'Installed source root differs')
  require(sha(e.touch(case/'provenance/installed-sources'/path[len('/opt/openfoam14/'):]))==digest,'Installed source snapshot hash differs')
 checks['actual_installed_solid_thermo_heat_flux_sources_hashes_match']=True
 input_manifest=e.finish();checks['all_native_and_source_inputs_unchanged_during_read_only_audit']=True
 arrays=dict(times_s=times,heat_times_s=times[1:],residual_times_s=times[1:],T=T,centres=centres,volumes=volumes,areas=volumes/.01,grid_shape=np.array([20,100]),**m,
  patch_names=np.array(PATCHES),wall_names=np.array(PATCHES),energy_J=energy,native_Q_inward_W=Q,native_Q_table_W=tableQ,storage_rate_W=storage,energy_balance_residual_W=balance,
  energy_balance_relative=relative,initial_linear_energy_residual=initial_residual,final_linear_energy_residual=final_residual,
  analytic_cell_average_T=analytical,analytic_error=error,analytic_max_error_K=max_error,analytic_RMS_error_K=rms,max_y_spread_K=spread)
 out.mkdir(parents=True);np.savez_compressed(out/'fields.npz',**arrays)
 all_pass=all(checks.values())
 record=dict(solver='OpenFOAM Foundation14 solid',actual_settings=declaration,case=str(case),actual_native_mesh_quality=quality,actual_native_patches=patches,
  actual_cell_mapping=mapping_report,actual_times_s=times.tolist(),actual_time_metadata=metadata,actual_saved_temperature_states=101,actual_saved_temperature_values=int(T.size),
  native_volume_field='Vc',native_initial_volume_field='0/Vc',native_final_volume_field='200/Vc',physical_time=True,flow_equations_solved=False,
  field_ranges=dict(Tmin_K=float(T.min()),Tmax_K=float(T.max()),final_Tmin_K=float(T[-1].min()),final_Tmax_K=float(T[-1].max()),max_y_spread_K=float(spread.max())),
  residuals=dict(definition='Printed final PCG energy linear-system residual for each of100actual Euler steps; initial residual recorded separately and is not the gate',
   maximum_printed_initial_energy_residual=float(initial_residual.max()),maximum_printed_final_energy_residual=float(final_residual.max()),
   terminal_initial_energy_residual=float(initial_residual[-1]),terminal_final_energy_residual=float(final_residual[-1]),norm_equivalent_to_fluent=False),
  native_heat_rates=dict(patch_order=list(PATCHES),convention='positive INTO solid; wallHeatFlux field W/m2 integrated once over actual nativefaceareas; native table Q alreadyW',
   table=str(logs[0]),same_time_identical_duplicate_rows=duplicates,maximum_field_table_scaled_Q_difference=float(cross.max()),field_table_scaled_Q_target=1e-5,
   maximum_abs_insulated_Q_W=float(np.max(np.abs(Q[:,2:]))),final_Q_inward_W=Q[-1].tolist(),final_energy_J=float(energy[-1])),
  global_energy_balance=dict(definition='Euler storage=sumrhoCvV(Tn-Tprev)/2s; residualstorage-sumQin; denominator max(sumabsQwalls,absstorage,Qref2W)',
   target=1e-5,maximum_relative_error=float(relative.max()),maximum_absolute_residual_W=float(np.abs(balance).max()),final_relative_error=float(relative[-1]),
   all_cell_thermal_flux_conservation_claimed=False),
  analytical_verification=dict(reference_sha256=REFERENCE_SHA,cell_average=True,reference_terms=80,series_check_terms=160,max_series_difference_K=max(series_delta),
   final_t200_max_error_K=float(max_error[-1]),final_t200_volume_weighted_RMS_error_K=float(rms[-1]),
   all_time_max_error_K=float(max_error.max()),all_time_max_error_time_s=int(times[int(np.argmax(max_error))]),all_time_accuracy_gate=False,
   final_max_error_target_K=.5,final_RMS_error_target_K=.25,initial_T_reference_explicit300=True),
  checks=checks,check_results=checks,failed_checks=[n for n,v in checks.items() if not v],all_checks_passed=all_pass,all_declared_checks_passed=all_pass,
  audited_fields_npz='fields.npz',audited_fields_npz_sha256=sha(out/'fields.npz'),physical_validation=False,mesh_independence_established=False,time_step_independence_established=False,
  limitations='Ideal generic constant-property 1D conduction; boundary jump initially incompatible and early transient coarsely resolved by Euler2s; all-time errors retained',
  no_invented_initial_heat_rate_or_residual=True)
 json_new(out/'checks.json',record);json_new(out/'source-manifest.json',input_manifest)
 with (out/'time-history.csv').open('x',encoding='ascii',newline='') as stream:
  w=csv.writer(stream);w.writerow(['time_s','energy_J','Tmin_K','Tmax_K','max_y_spread_K','analytic_max_error_K','analytic_RMS_error_K'])
  for i,t in enumerate(times):w.writerow([int(t),float(energy[i]),float(T[i].min()),float(T[i].max()),float(spread[i]),float(max_error[i]),float(rms[i])])
 with (out/'energy-step-history.csv').open('x',encoding='ascii',newline='') as stream:
  w=csv.writer(stream);w.writerow(['time_s',*['Qin_'+n+'_W' for n in PATCHES],'storage_rate_W','energy_balance_residual_W','energy_balance_relative','initial_linear_energy_residual','final_linear_energy_residual'])
  for i,t in enumerate(times[1:]):w.writerow([int(t),*Q[i].tolist(),float(storage[i]),float(balance[i]),float(relative[i]),float(initial_residual[i]),float(final_residual[i])])
 for path in (source,prep,reference):copied(path,out/('executed-'+path.name))
 print(json.dumps(dict(output=str(out),all_checks_passed=all_pass,failed_checks=record['failed_checks'],fields_sha256=record['audited_fields_npz_sha256']),indent=2))
 return 0 if all_pass else 2
if __name__=='__main__':sys.exit(main())

