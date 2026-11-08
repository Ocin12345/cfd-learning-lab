"""Audit every logged step and saved native state of refined sloshing; launch no solver."""
import argparse,csv,hashlib,json,re
from pathlib import Path
import h5py,numpy as np
from analyse_nozzle_fluent_v005 import common_mapping
from sloshing_dt_common_v001 import NX,NY,DX,DT,MATERIALS,RAW_TARGETS,alpha_group,cell_values
NUM=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
P=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def field(path,components):
    text=path.read_text(encoding='ascii');n=NX*NY
    match=re.search(r'internalField\s+uniform\s+([^;]+);',text)
    if match:
        v=np.fromstring(match[1].replace('(',' ').replace(')',' '),sep=' ');assert len(v)==components
        return np.tile(v,(n,1)).reshape(n,components)
    match=re.search(r'internalField\s+nonuniform\s+List<\w+>\s+(\d+)\s*\(',text);assert match and int(match[1])==n,path
    closing=re.search(r'\n\s*\)\s*;',text[match.end():]);assert closing,path
    v=np.fromstring(text[match.end():match.end()+closing.start()].replace('(',' ').replace(')',' '),sep=' ');assert v.size==n*components,path
    return v.reshape(n,components)
def expected(xy,tilt):return np.clip((.05+tilt*(xy[:,0]-.1)-(xy[:,1]-DX/2))/DX,0,1)
def saved_check(step,time,alpha,velocity,expect,mode,pressure):
    assert np.isfinite(alpha).all() and np.isfinite(velocity).all() and np.isfinite(pressure).all()
    area=float(alpha.sum()*DX**2);speed=float(np.linalg.norm(velocity,axis=1).max());difference=float(np.max(abs(alpha-expect)))
    if step==0:assert difference<1e-10 and abs(area-.01)<1e-10
    bounds=bool(alpha.min()>=-1e-6 and alpha.max()<=1+1e-6);volume=abs(area/.01-1)<=.005
    quiet=mode!='flat' or (speed<=1e-4 and difference<=.01)
    return dict(step=step,time_s=time,water_area_m2=area,max_speed_m_s=speed,alpha_min=float(alpha.min()),alpha_max=float(alpha.max()),alpha_bounds_pass=bounds,volume_pass=bool(volume),quiet_pass=bool(quiet),pressure_min=float(pressure.min()),pressure_max=float(pressure.max()))
def fluent(run,mode,steps,interval):
    setup=json.loads((run/'setup.json').read_text());assert setup['case_mode']==mode and setup['cells']==[NX,NY]
    assert setup['total_steps']==steps and setup['frame_interval']==interval and setup['time_step_s']==DT
    tilt=0 if mode=='flat' else .025;assert setup['initial_tilt']==tilt
    assert setup['model']['model']=='vof' and setup['model']['vof_parameters']['vof_formulation']=='explicit'
    assert setup['model']['phase_interaction']['forces']['surface_tension_model'] is False
    assert setup['operating_conditions']['gravity']['components']==[0,-9.81]
    assert setup['vof_patch_options']['patch_reconstructed_interface'] is False and setup['vof_patch_options']['use_volumetric_smoothing'] is False
    for name,(rho,mu) in MATERIALS.items():assert setup['materials'][name]['density']['value']==rho and setup['materials'][name]['viscosity']['value']==mu
    limits=RAW_TARGETS if mode=='flat' else {k:1e-4 for k in RAW_TARGETS}
    reporting=setup['residual_reporting'];assert reporting['absolute_targets']==limits
    assert reporting['options']['normalize'] is False and reporting['options']['residual_values']['scale_residuals']==(mode!='flat')
    for name,target in limits.items():assert reporting['equations'][name]['check_convergence'] is True and reporting['equations'][name]['absolute_criteria']==target
    assert json.loads((run/'status.json').read_text())['stage']=='closed'
    text=(run/'solver-transcript.txt').read_text(errors='replace')
    quality=float(re.search(r'Minimum Orthogonal Quality\s*=\s*('+NUM+')',text)[1]);assert quality>.2
    assert float(re.search(r'minimum volume \(m3\):\s*('+NUM+')',text)[1])>0
    boundaries=list(re.finditer(r'Flow time = ('+NUM+r')s, time step = (\d+)',text));assert len(boundaries)==steps
    pattern=re.compile(r'^\s*(\d+)\s+('+NUM+r')\s+('+NUM+r')\s+('+NUM+r')\s+\d+:\d+:\d+\s+\d+\s*$',re.M)
    histories=[];start=0;previous=0
    for step,b in enumerate(boundaries,1):
        assert int(b[2])==step and abs(float(b[1])-step*DT)<1e-10
        rows=list(pattern.finditer(text,start,b.start()));assert rows
        last=rows[-1];iteration=int(last[1]);values=[float(last[k]) for k in [2,3,4]];assert np.isfinite(values).all()
        assert 0<iteration-previous<=400
        histories.append(dict(step=step,time_s=float(b[1]),residuals=dict(zip(limits,values)),inner_iterations=iteration-previous,passed=all(v<=t for v,t in zip(values,limits.values()))));previous=iteration;start=b.end()
    co=[float(x) for x in re.findall(r'Global Courant Number \[Explicit VOF Criteria\]\s*:\s*('+NUM+')',text)];assert co
    with np.load(run/'geometry.npz') as g:fm=common_mapping(run/'active.cas.h5',g);xy=g['centres']
    assert len(np.unique(fm))==NX*NY
    expect=expected(xy,tilt)
    with (run/'frames.csv').open(newline='') as stream:ledger=list(csv.DictReader(stream))
    assert [int(r['step']) for r in ledger]==list(range(0,steps+1,interval))
    from ansys.fluent.core.offline import _lispy
    records=[];sources=[run/name for name in ['setup.json','geometry.npz','active.cas.h5','tank-final.cas.h5','tank-final.dat.h5','frames.csv','solver-transcript.txt','initial-field-check.json','executed-driver.py']]
    for row in ledger:
        step=int(row['step']);path=run/row['file']
        with h5py.File(path) as f:
            clock={v[0]:v[1] for v in _lispy.parse(f['settings/Data Variables'][0].decode())[1]}
            assert float(clock['time-step'])==step and abs(float(clock['flow-time'])-step*DT)<1e-10 and abs(float(clock['physical-time-step'])-DT)<1e-14
            assert abs(float(row['time_s'])-float(clock['flow-time']))<1e-10
            raw=cell_values(alpha_group(f));a=np.empty_like(raw);a[fm]=raw
            uv=[]
            for name in ['SV_U','SV_V']:
                raw=cell_values(f['results/1/phase-1/cells/'+name]);aligned=np.empty_like(raw);aligned[fm]=raw;uv.append(aligned)
            pressure=cell_values(f['results/1/phase-1/cells/SV_P'])
        records.append(saved_check(step,step*DT,a,np.column_stack(uv),expect,mode,pressure));sources.append(path)
    return histories,records,max(co),dict(minimum_orthogonal_quality=quality,residual_reporting='unscaled' if mode=='flat' else 'scaled',residual_targets=limits),sources

def openfoam(run,mode,steps,interval):
    d=json.loads((run/'declared-settings.json').read_text());assert d['case_mode']==mode and d['cells']==[NX,NY] and d['steps']==steps and d['frame_interval']==interval and d['dt_s']==DT
    tilt=0 if mode=='flat' else .025;assert d['initial_tilt']==tilt and d['momentum_predictor']
    assert (run/'run-status.txt').read_text().strip()=='complete'
    mesh=(run/'mesh-check.log').read_text();assert 'Mesh OK.' in mesh and re.search(r'\bcells:\s+12800\b',mesh)
    text=(run/'solver.log').read_text();assert re.search(r'^End\s*$',text,re.M) and 'FOAM FATAL' not in text
    parts=re.split(r'^Time = ([0-9.eE+-]+)s?\s*$',text,flags=re.M);histories=[];limits={'Ux':1e-4,'Uy':1e-4,'p_rgh':1e-4}
    for index in range(1,len(parts),2):
        t=float(parts[index]);chunk=parts[index+1];assert 'PIMPLE: Iteration 8' in chunk;last=chunk.rsplit('PIMPLE: Iteration 8',1)[1];res={}
        for m in re.finditer(r'Solving for (\w+), Initial residual = ('+NUM+r'), Final residual = ('+NUM+')',last):
            if m[1] in limits:res[m[1]]=max(res.get(m[1],0),float(m[2]))
        assert set(res)==set(limits)
        alpha_rows=re.findall(r'Min\(alpha.water\)\s*=\s*('+NUM+r').*?Max\(alpha.water\)\s*=\s*('+NUM+')',last);assert alpha_rows
        lo,hi=map(float,alpha_rows[-1]);assert np.isfinite([lo,hi,*res.values()]).all()
        histories.append(dict(step=len(histories)+1,time_s=t,residuals=res,alpha_min=lo,alpha_max=hi,passed=all(res[k]<=limits[k] for k in limits) and lo>=-1e-6 and hi<=1+1e-6))
    assert len(histories)==steps and np.allclose([r['time_s'] for r in histories],np.arange(1,steps+1)*DT,atol=1e-10,rtol=0)
    co=[float(x) for x in re.findall(r'^Interface Courant Number mean: '+NUM+r' max: ('+NUM+')',text,re.M)];assert co
    xy=field(run/'0/C',3)[:,:2];ij=np.rint(xy/DX-.5).astype(int)
    assert np.all((ij>=0)&(ij<[NX,NY])) and len(np.unique(ij[:,0]+NX*ij[:,1]))==NX*NY and np.max(abs(xy-(ij+.5)*DX))<1e-10
    expect=expected(xy,tilt);records=[]
    sources=[run/name for name in ['declared-settings.json','solver.log','mesh-check.log','0/C','system/fvSolution','system/fvSchemes','system/controlDict','system/blockMeshDict','constant/g','constant/momentumTransport','constant/phaseProperties','constant/physicalProperties.water','constant/physicalProperties.air']]
    for step in range(0,steps+1,interval):
        folder=run/f'{step*DT:.8g}';a=field(folder/'alpha.water',1).ravel();uv=field(folder/'U',3);pressure=field(folder/'p_rgh',1).ravel()
        if step:
            p=folder/'uniform/time';time=p.read_text();v=re.search(r'^value\s+('+NUM+r')\s*;',time,re.M);n=re.search(r'^index\s+(\d+)\s*;',time,re.M);dt=re.search(r'^deltaT\s+('+NUM+r')\s*;',time,re.M)
            assert v and n and dt and abs(float(v[1])-step*DT)<1e-10 and int(n[1])==step and abs(float(dt[1])-DT)<1e-14;sources.append(p)
        records.append(saved_check(step,step*DT,a,uv,expect,mode,pressure));sources.extend([folder/'alpha.water',folder/'U',folder/'p_rgh'])
    return histories,records,max(co),dict(checkMesh_pass=True,residual_reporting='last-outer initial native',residual_targets=limits),sources


def audit(solver,run,out,steps=80,interval=1):
 from sloshing_courant_v001 import courant_records
 from sloshing_fixed_settings_v001 import compare_fixed
 from sloshing_dt_common_v001 import B
 assert not out.exists()
 history,records,co,extra,sources=(fluent if solver=='fluent' else openfoam)(run,'qualification',steps,interval)
 text=(run/('solver-transcript.txt' if solver=='fluent' else 'solver.log')).read_text(errors='replace')
 coverage=courant_records(text,solver,DT,0,steps)
 if solver=='fluent':
  a=json.loads((run/'fixed-settings-comparison.json').read_text());assert a['passed']
  with h5py.File(B/'fluent/tank-0000.dat.h5') as f,h5py.File(run/'state-0000.dat.h5') as g:
   for name in ['SV_U','SV_V','SV_P']:
    assert np.array_equal(cell_values(f['results/1/phase-1/cells/'+name]),cell_values(g['results/1/phase-1/cells/'+name])),name
   assert np.array_equal(cell_values(alpha_group(f)),cell_values(alpha_group(g)))
  sources += [run/'fixed-settings-comparison.json',run/'loaded-fixed-settings.json']
 else:
  for part in ['processor0','processor1']:
   p=run/part/f'{steps*DT:.8g}'/'uniform/time';raw=p.read_text()
   for key,value in [('value',steps*DT),('index',steps),('deltaT',DT)]:assert abs(float(re.search(r'^\s*'+key+r'\s+('+NUM+');',raw,re.M)[1])-value)<1e-10
   sources.append(p)
  for key in ['U','p_rgh','alpha.water']:
   assert np.array_equal(field(run/'0'/key,3 if key=='U' else 1),field(B/'openfoam/0'/key,3 if key=='U' else 1))
  from sloshing_phi_v001 import phi_wall_flux
  flux=[]
  for k in range(1,steps+1):
   p=run/f'{k*DT:.8g}'/'phi';q=phi_wall_flux(p.read_text(),480,25360)/.01;assert q<=1e-12;flux.append(q);sources.append(p)
  extra['max_absolute_wall_volume_flux_per_depth_m2_s']=max(flux)
  # Exact copied native mesh/partition inputs remain bound by launch manifest.
 fields=all(r['alpha_bounds_pass'] and r['volume_pass'] for r in records);passed=fields and all(r['passed'] for r in history) and co<=.25
 sources += [Path(__file__),P/'sloshing_dt_common_v001.py',P/'sloshing_courant_v001.py',P/'SPEC.md',P/'qualification-inputs-v003.json']
 result=dict(numerical_pass=bool(passed),qualification_only=True,solver=solver,steps=steps,dt_s=DT,states=len(records),end_time_s=steps*DT,records=records,step_records=history,courant_records=coverage,extra=extra,exact_original_initial_fields=True,physical_validation=False,source_sha256={str(p):sha(p) for p in sources})
 out.write_text(json.dumps(result,indent=2)+'\n');assert passed,'Native qualification numerical MISS; unchanged criteria'
 return result
if __name__=='__main__':
 import sys
 a=audit(sys.argv[1],Path(sys.argv[2]),Path(sys.argv[3]));print(json.dumps({k:v for k,v in a.items() if k not in ['records','step_records','courant_records','source_sha256']},indent=2))
