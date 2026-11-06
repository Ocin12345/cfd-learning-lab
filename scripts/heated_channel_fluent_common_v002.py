"""Heated-channel pilot native readers/gates; no solver launch or production loop."""
from pathlib import Path
import hashlib,json,re,math
import h5py,numpy as np
from analyse_nozzle_fluent_v005 import chunk_array,common_mapping
WALLS=('top','bottom');BOUNDARIES=('inlet','outlet','top','bottom');EQUATIONS=('continuity','x-velocity','y-velocity')
COUNTS=dict(interior=4488,inlet=24,outlet=24,top=96,bottom=96)
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
TARGETS_SHA='0a163061a9836c096d8b30b16c81410285e29bf3507f10d231619c43b64cab02'
TARGETS=dict(flow_iterations=2000,flow_checkpoints=[500,1000,1500,1800,2000],flow_residual_max=1e-6,flow_mass_relative_max=1e-5,
 flow_dimensionless_divergence_max=1e-5,flow_dimensionless_divergence_RMS_max=1e-6,flow_Dh_over_U_s=.2,flow_velocity1800to2000_L2_max=.001,flow_stability_denominator='final2000field',minimum_orthogonal_quality=.15,
 pilot_thermal_steps=1,pilot_final_time_s=.005,time_step_s=.005,thermal_inner_cap=60,energy_residual_max=1e-8,temperature_min_K=300.,temperature_max_K=350.,temperature_bound_tolerance_K=.01,
 frozen_flow_absolute_tolerance_SI=1e-12,global_Euler_energy_relative_max=1e-5,reconstructed_PDE_energy_relative_max=.005,power_normalization_floor_W_per_m=5000.,
 rho_kg_m3=1.,mu_Pa_s=.01,cp_J_kg_K=1000.,k_W_m_K=10.,alpha_m2_s=.01,reference_depth_m=1.,inlet_diffusion=True,
 production_steps_NOT_AUTHORIZED=600,production_final_time_s_NOT_AUTHORIZED=3.,physical_validation=False,mesh_independence=False,time_step_independence=False)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')
def source_files():
 p=Path(__file__).resolve().parent
 return {n:p/n for n in ('heated_channel_fluent_common_v002.py','analyse_nozzle_fluent_v005.py')}
def geometry_check(g):
 assert g['quads'].shape==(2304,4) and g['nodes'].shape==(2425,2) and np.array_equal(g['grid_shape'],[24,96])
 assert np.isfinite(g['areas']).all() and np.all(g['areas']>0) and abs(g['areas'].sum()-.1)<1e-12
 assert np.array_equal(g['ix'],np.tile(np.arange(96),24)) and np.array_equal(g['iy'],np.repeat(np.arange(24),96))
 return True
def native_time(data):
 from ansys.fluent.core.offline import _lispy
 with h5py.File(data,'r') as f:values={v[0]:v[1] for v in _lispy.parse(f['settings/Data Variables'][0].decode())[1]}
 result={k:float(values[k]) for k in ('flow-time','time-step','physical-time-step')};assert all(map(math.isfinite,result.values())) and result['time-step']==int(result['time-step'])
 return result
def outward_geometry(raw_area,face_mid,owner_centres):
 # Geometry only: derive outward normal from actual native owner mapped to source centre.
 # No native flux, heat power or field value is an input or changed here.
 raw=np.asarray(raw_area,float);mid=np.asarray(face_mid,float);owner=np.asarray(owner_centres,float)
 assert raw.shape==mid.shape==owner.shape and raw.ndim==2 and raw.shape[1]==2
 alignment=np.sum((mid-owner)*raw,axis=1);assert np.isfinite(alignment).all() and np.all(abs(alignment)>1e-24)
 flipped=alignment<0;area=raw.copy();area[flipped]*=-1
 assert np.all(np.sum((mid-owner)*area,axis=1)>0)
 return area,flipped
def snapshot(case,data,g,thermal=False):
 geometry_check(g);mapping=common_mapping(case,g);result=dict(mapping=mapping)
 with h5py.File(case,'r') as c,h5py.File(data,'r') as d:
  root=c['meshes/1'];faces=root['faces'];zt=faces['zoneTopology'];names=b''.join(np.asarray(zt['name'][:]).ravel()).decode().split(';')
  lo=np.asarray(zt['minId'][:]).ravel().astype(int)-1;hi=np.asarray(zt['maxId'][:]).ravel().astype(int);types=np.asarray(zt['zoneType'][:]).ravel().astype(int)
  assert set(names)==set(COUNTS);zones={n:(int(a),int(b)) for n,a,b in zip(names,lo,hi)}
  assert all(b-a==COUNTS[n] for n,(a,b) in zones.items()) and dict(zip(names,map(int,types)))==dict(interior=2,inlet=10,outlet=5,top=3,bottom=3)
  c0=chunk_array(faces['c0'],4728).astype(int);c1=np.zeros(4728,int);seen=np.zeros(4728,bool)
  for ds in faces['c1'].values():
   a=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;b=int(np.asarray(ds.attrs['maxId']).ravel()[0]);v=np.asarray(ds[:]).reshape(-1)
   assert not seen[a:b].any() and np.array_equal(v,v.astype(int));c1[a:b]=v;seen[a:b]=True
  assert np.all((c0>=1)&(c0<=2304)) and np.all((c1>=0)&(c1<=2304))
  nodes=chunk_array(root['nodes/coords'],2425,(2,));face_nodes=np.empty((4728,2),int);seen=np.zeros(4728,bool)
  for ds in faces['nodes'].values():
   a=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;b=int(np.asarray(ds.attrs['maxId']).ravel()[0]);ids=np.asarray(ds['nodes'][:]).reshape(-1,2).astype(int)-1
   assert len(ids)==b-a and not seen[a:b].any() and np.all(np.asarray(ds['nnodes'][:])==2);face_nodes[a:b]=ids;seen[a:b]=True
  assert seen.all();mid=nodes[face_nodes].mean(axis=1);edge=nodes[face_nodes[:,1]]-nodes[face_nodes[:,0]];raw_area=np.column_stack([edge[:,1],-edge[:,0]])
  owner_centres=g['centres'][mapping[c0-1]];area,flipped=outward_geometry(raw_area,mid,owner_centres)
  orientation_proof=dict(geometric_normal_flip_count=int(flipped.sum()),total_native_faces=4728,all_geometric_normal_dot_face_mid_minus_owner_positive=True,minimum_positive_alignment=float(np.sum((mid-owner_centres)*area,axis=1).min()),raw_native_flux_and_heat_fields_unmodified=True)
  cell=d['results/1/phase-1/cells'];face=d['results/1/phase-1/faces'];native={};common={}
  pairs=[('u','SV_U'),('v','SV_V'),('p','SV_P'),('rho','SV_DENSITY'),('mu','SV_MU_LAM')]
  if thermal:pairs += [('T','SV_T'),('h','SV_H')]
  for key,name in pairs:
   native[key]=chunk_array(cell[name],2304);common[key]=np.empty(2304);common[key][mapping]=native[key]
  flux=chunk_array(face['SV_FLUX'],4728)
  result.update(native=native,common=common,flux=flux,c0=c0,c1=c1,zones=zones,face_nodes=face_nodes,face_area=area,raw_face_area=raw_area,face_orientation_proof=orientation_proof,face_mid=mid,cell_field_names=list(cell),face_field_names=list(face))
 return result
def frozen_compare(ref,new):
 assert np.array_equal(ref['face_nodes'],new['face_nodes']) and np.array_equal(ref['c0'],new['c0']) and np.array_equal(ref['c1'],new['c1']) and np.array_equal(ref['face_area'],new['face_area'])
 errors={k:float(abs(ref['common'][k]-new['common'][k]).max()) for k in ('u','v','p')};errors['all_face_mass_flux']=float(abs(ref['flux']-new['flux']).max())
 return dict(maximum_absolute_changes_SI=errors,passed=all(v<=1e-12 for v in errors.values()),definition='Every mapped native U,V,p and all4728 native face mass fluxes; no reconstructed flux')
def terminal_residuals(data,equations):
 result={};indices={}
 with h5py.File(data,'r') as f:
  for name in equations:
   group=f['results/residuals/phase-1/'+name];i=np.asarray(group['iterations'][:]).ravel();a=np.asarray(group['data'][:]);assert len(i) and np.array_equal(i,i.astype(int)) and np.all(np.diff(i)>0) and a.shape[0]==len(i) and a.shape[1]>=2 and np.isfinite(a[:,:2]).all()
   assert a[-1,1]>0;value=float(a[-1,0]/a[-1,1]);assert math.isfinite(value) and value>=0;result[name]=value;indices[name]=i.astype(int).tolist()
 return result,indices
def flow_gates(final,baseline,g,data):
 r,index=terminal_residuals(data,EQUATIONS);assert all(v==list(range(1,2001)) for v in index.values())
 flux=final['flux'];c0=final['c0'];c1=final['c1'];sums=np.zeros(2304);np.add.at(sums,c0-1,flux);inside=c1>0;np.add.at(sums,c1[inside]-1,-flux[inside]);vol=g['areas'][final['mapping']];div=sums/vol*.2
 q={n:float(flux[a:b].sum()) for n,(a,b) in final['zones'].items() if n in BOUNDARIES};net=sum(q.values());inlet=-q['inlet'];assert inlet>0
 b=np.column_stack([baseline['common']['u'],baseline['common']['v']]);f=np.column_stack([final['common']['u'],final['common']['v']]);den=float(np.sum(g['areas'][:,None]*f*f));assert den>0;change=float(np.sqrt(np.sum(g['areas'][:,None]*(f-b)**2)/den))
 a,b=final['zones']['outlet'];outlet_no_backflow=bool(np.all(flux[a:b]>=0));wallflux={n:float(abs(flux[a:b]).max()) for n,(a,b) in final['zones'].items() if n in WALLS}
 metrics=dict(final_scaled_residuals=r,native_face_signed_mass_flow_kg_s=q,actual_inlet_mass_flow_kg_s=inlet,relative_mass_imbalance=abs(net)/inlet,maximum_dimensionless_divergence=float(abs(div).max()),RMS_dimensionless_divergence=float(np.sqrt(np.sum(vol*div*div)/vol.sum())),velocity1800to2000_relativeL2_FINAL_denominator=change,outlet_no_backflow=outlet_no_backflow,wall_maximum_absolute_face_mass_flux_kg_s=wallflux)
 checks=dict(complete_native_flow_residual_history1to2000=True,flow_final_residuals_le1e_6=all(v<=1e-6 for v in r.values()),flow_mass_relative_le1e_5=metrics['relative_mass_imbalance']<=1e-5,flow_local_divergence_max_le1e_5=metrics['maximum_dimensionless_divergence']<=1e-5,flow_local_divergence_RMS_le1e_6=metrics['RMS_dimensionless_divergence']<=1e-6,flow_velocity_change_le1e_3=change<=.001,actual_inlet_inward_outlet_outward=bool(q['inlet']<0 and q['outlet']>0 and outlet_no_backflow),zero_wall_normal_flux=all(v<=1e-12 for v in wallflux.values()),native_rho_mu_constant=bool(np.allclose(final['common']['rho'],1,atol=1e-12,rtol=0) and np.allclose(final['common']['mu'],.01,atol=1e-12,rtol=0)))
 return dict(checks=checks,all_checks_passed=all(checks.values()),metrics=metrics)
def flux_report(path,unit):
 text=Path(path).read_text(encoding='utf-8',errors='replace');assert '['+unit+']' in text
 result={}
 for name in (*BOUNDARIES,'Net'):
  matches=re.findall(r'^\s*'+re.escape(name)+r'\s+('+NUMBER+r')\s*$',text,re.M);assert matches,'Missing native flux report row '+name;result[name]=float(matches[-1]);assert math.isfinite(result[name])
 return result
def boundary_array(data,field,zones):
 values=np.zeros(4728);seen=np.zeros(4728,bool)
 with h5py.File(data,'r') as f:
  for ds in f['results/1/phase-1/faces/'+field].values():
   a=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;b=int(np.asarray(ds.attrs['maxId']).ravel()[0]);v=np.asarray(ds[:]).reshape(-1);assert 0<=a<b<=4728 and len(v)==b-a and not seen[a:b].any() and np.isfinite(v).all();values[a:b]=v;seen[a:b]=True
 boundary=np.zeros(4728,bool)
 for n in BOUNDARIES:a,b=zones[n];boundary[a:b]=True
 assert np.all(seen[boundary]),'Every physical boundary value required'
 return values,{n:float(values[a:b].sum()) for n,(a,b) in zones.items() if n in BOUNDARIES},dict(saved_faces=int(seen.sum()),physical_boundary_faces=int(boundary.sum()),only_boundary_storage=bool(np.array_equal(seen,boundary)))
def thermal_metrics(initial,final,g,initial_data,final_data,total_report,mass_report,initial_total_report=None):
 T0=initial['common']['T'];T=final['common']['T'];h0=initial['common']['h'];h=final['common']['h'];assert np.max(abs(T0-300))<=1e-10 and np.ptp(h0)<=1e-6
 h300=float(h0.mean())
 initial_q_array,initial_q,initial_coverage=boundary_array(initial_data,'SV_HEAT_FLUX',initial['zones'])
 initial_q_report=flux_report(initial_total_report,'W') if initial_total_report is not None else None
 initial_report_match=initial_q_report is not None and all(abs(initial_q[n]-initial_q_report[n])<=1e-6*max(1,abs(initial_q_report[n])) for n in BOUNDARIES)
 initial_open_expected={n:-float(initial['flux'][a:b].sum())*h300 for n,(a,b) in initial['zones'].items() if n in ('inlet','outlet')}
 initial_open_errors={n:initial_q[n]-initial_open_expected[n] for n in initial_open_expected}
 initial_open_reference_match=initial_report_match and all(abs(initial_open_errors[n])<=1e-6*max(1,abs(initial_open_expected[n])) for n in initial_open_expected)
 h_error=float(abs((h-h300)-1000*(T-300)).max());q_array,q_raw,coverage=boundary_array(final_data,'SV_HEAT_FLUX',final['zones']);qr=flux_report(total_report,'W');mr=flux_report(mass_report,'kg/s')
 representation_errors={n:abs(q_raw[n]-qr[n]) for n in BOUNDARIES};report_match=all(v<=1e-6*max(1,abs(qr[n])) for n,v in representation_errors.items()) and abs(sum(q_raw.values())-qr['Net'])<=1e-6*max(1,abs(qr['Net']))
 mout={n:float(final['flux'][a:b].sum()) for n,(a,b) in final['zones'].items() if n in BOUNDARIES};mass_match=all(abs(-mout[n]-mr[n])<=1e-6*max(1,abs(mr[n])) for n in BOUNDARIES)
 q_above300={n:q_raw[n]+h300*mout[n] for n in BOUNDARIES};E0=float(np.sum(g['areas']*1000*(T0-300)));E=float(np.sum(g['areas']*1000*(T-300)));storage=(E-E0)/.005;normalization=max(sum(abs(v) for v in q_above300.values()),abs(storage),5000.);defect=storage-sum(q_above300.values());relative=abs(defect)/normalization
 # Separate orthogonal two-point boundary/upwind PDE reconstruction, never called native heat.
 nativeT=final['native']['T'];nativecentres=g['centres'][final['mapping']];reconstructed={};components={}
 for n in BOUNDARIES:
  a,b=final['zones'][n];owner=final['c0'][a:b]-1;m=final['flux'][a:b];ownerT=nativeT[owner];faceA=np.linalg.norm(final['face_area'][a:b],axis=1);norm=final['face_area'][a:b]/faceA[:,None];distance=np.sum((final['face_mid'][a:b]-nativecentres[owner])*norm,axis=1);assert np.all(distance>0)
  incomingT=300.;upwind=np.where(m>=0,ownerT,incomingT);convection=-m*1000*(upwind-300)
  diffusion=np.zeros(b-a) if n=='outlet' else 10*((350. if n in WALLS else 300.)-ownerT)/distance*faceA
  reconstructed[n]=float(np.sum(convection+diffusion));components[n]=dict(diffusion_in_W=float(diffusion.sum()),convection_in_above300_W=float(convection.sum()))
 pde_defect=storage-sum(reconstructed.values());pde_relative=abs(pde_defect)/normalization
 a,b=final['zones']['outlet'];weights=final['flux'][a:b];owner=final['c0'][a:b]-1;assert weights.sum()>0;bulk=float(np.sum(weights*nativeT[owner])/weights.sum())
 residuals,indices=terminal_residuals(final_data,('energy',));energy_index=indices['energy'];assert 1<=len(energy_index)<=60
 time0=native_time(initial_data);time1=native_time(final_data)
 checks=dict(native_initial_open_boundary_h300_and_inward_enthalpy_reference_proven=initial_open_reference_match,native_initial_T300=bool(np.max(abs(T0-300))<=1e-10),native_T_bounds=bool(T.min()>=299.99 and T.max()<=350.01),temperature_changed_first_step=bool(np.max(abs(T-T0))>1e-10),native_h300_and_cpDeltaT_match=h_error<=1e-6,native_total_heat_report_HDF_integrated_W_all_boundaries=report_match,native_mass_report_HDF_sign_all_boundaries=mass_match,native_hotwall_Qin_positive=all(q_above300[n]>0 for n in WALLS),native_global_Euler_energy_relative_le1e_5=relative<=1e-5,separate_reconstructed_PDE_energy_relative_le0p005=pde_relative<=.005,native_energy_final_scaled_le1e_8=residuals['energy']<=1e-8,actual_one_physical_step_dt0p005=bool(time0['flow-time']==0 and time0['time-step']==0 and abs(time1['flow-time']-.005)<=1e-12 and time1['time-step']==1 and abs(time1['physical-time-step']-.005)<=1e-12))
 return dict(checks=checks,all_checks_passed=all(checks.values()),metrics=dict(native_time_initial=time0,native_time_step1=time1,native_energy_iteration_indices=energy_index,native_energy_final_scaled_residual=residuals['energy'],temperature_min_K=float(T.min()),temperature_max_K=float(T.max()),native_h300_J_kg=h300,native_initial_heat_report_HDF_match=initial_report_match,native_initial_raw_Qin_W=initial_q,native_initial_heat_report_W=initial_q_report,native_initial_open_h300_advection_expected_W=initial_open_expected,native_initial_open_h300_reference_errors_W=initial_open_errors,native_h_cpDeltaT_max_error_J_kg=h_error,SV_HEAT_FLUX_coverage=coverage,native_raw_Qin_W=q_raw,native_total_heat_report_W=qr,native_heat_report_HDF_absolute_errors_W=representation_errors,native_mass_outward_HDF_kg_s=mout,native_mass_inward_report_kg_s=mr,native_reference_corrected_Qin_above300_W_per_m=q_above300,energy_initial_J_per_m=E0,energy_step1_J_per_m=E,storage_W_per_m=storage,native_Euler_defect_W_per_m=defect,native_Euler_normalization_W_per_m=normalization,native_Euler_relative_defect=relative,reconstructed_PDE_Qin_above300_W_per_m=reconstructed,reconstructed_PDE_components_W_per_m=components,reconstructed_PDE_relative_defect=pde_relative,direct_native_minus_reconstructed_boundary_power_W_per_m={n:q_above300[n]-reconstructed[n] for n in BOUNDARIES},actual_outward_flux_weighted_outlet_bulk_temperature_K=bulk,energy_note='Raw native total heat includes convection if present; reference correction uses actualcellh300 and signednative mass flux. No extraenthalpy transport orarea multiplier added. Global conservation only.'))


def thermal_storage_inventory(data,zones):
 """Read actual saved storage without assigning unsupported field units."""
 result={'cell_fields':{},'face_fields':{},'unsolved_initial_heat_fields_are_not_a_transport_reference':True}
 with h5py.File(data,'r') as f:
  cells=f['results/1/phase-1/cells'];faces=f['results/1/phase-1/faces']
  for name in ('SV_T','SV_H','SV_T_M1','SV_ENERGY_M1'):
   if name in cells:
    a=chunk_array(cells[name],2304);result['cell_fields'][name]=dict(count=len(a),minimum=float(a.min()),maximum=float(a.max()),mean=float(a.mean()),peak_to_peak=float(np.ptp(a)),all_finite=bool(np.isfinite(a).all()))
  for name in ('SV_T','SV_H','SV_HEAT_FLUX','SV_HEAT_FLUX_SENSIBLE','SV_RAD_HEAT_FLUX','SV_RAD_HEAT_FLUX_EXTERIOR'):
   if name not in faces:continue
   values=np.zeros(4728);seen=np.zeros(4728,bool)
   for ds in faces[name].values():
    a=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;b=int(np.asarray(ds.attrs['maxId']).ravel()[0]);v=np.asarray(ds[:]).reshape(-1)
    assert 0<=a<b<=4728 and len(v)==b-a and not seen[a:b].any() and np.isfinite(v).all();values[a:b]=v;seen[a:b]=True
   result['face_fields'][name]=dict(saved_faces=int(seen.sum()),minimum=float(values[seen].min()),maximum=float(values[seen].max()),zones={n:dict(saved_faces=int(seen[a:b].sum()),zone_faces=b-a,sum_native_values=float(values[a:b][seen[a:b]].sum()),all_zone_faces_saved=bool(seen[a:b].all())) for n,(a,b) in zones.items() if n in BOUNDARIES},units='native storage; heat-flux W interpretation separately checked against each native report')
 return result
