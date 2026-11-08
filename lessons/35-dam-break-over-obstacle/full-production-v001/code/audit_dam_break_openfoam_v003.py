"""Audit native open-boundary dam-break clocks, fields, residuals and mass flux."""
from pathlib import Path
import argparse,hashlib,json,re
import numpy as np
from dam_break_common_v002 import *
from dam_break_admission_v004 import verify
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);ap.add_argument('out',type=Path);ap.add_argument('--mode',choices=['flat','collapse','production'],required=True);args=ap.parse_args();run=args.run;assert not args.out.exists()
    verify('openfoam-'+args.mode)
    declared=json.loads((run/'declared-settings.json').read_text());steps=10000 if args.mode=='production' else 400;interval=100 if args.mode=='production' else 20
    assert declared['mode']==args.mode and declared['cells']==NCELLS and declared['steps']==steps and declared['save_interval']==interval and declared['dt_s']==DT
    assert declared['surface_tension_N_m']==SIGMA and declared['static_atmosphere_gauge_Pa']==0 and declared['backflow_water_fraction']==0
    assert (run/'run-status.txt').read_text().strip()=='complete'
    mesh=(run/'mesh-check.log').read_text();assert 'Mesh OK.' in mesh and re.search(r'\bcells:\s+9072\b',mesh)
    text=(run/'solver.log').read_text();assert re.search(r'^End\s*$',text,re.M) and 'FOAM FATAL' not in text
    chunks=re.split(r'^Time = ('+NUMBER+r')s?\s*$',text,flags=re.M);history=[];targets={'Ux':1e-4,'Uy':1e-4,'p_rgh':1e-4}
    for index in range(1,len(chunks),2):
        t=float(chunks[index]);chunk=chunks[index+1];assert 'PIMPLE: Iteration 8' in chunk;last=chunk.rsplit('PIMPLE: Iteration 8',1)[1];res={}
        for m in re.finditer(r'Solving for (\w+), Initial residual = ('+NUMBER+r'), Final residual = ('+NUMBER+')',last):
            if m[1] in targets:res[m[1]]=max(res.get(m[1],0),float(m[2]))
        assert set(res)==set(targets)
        rows=re.findall(r'Min\(alpha.water\)\s*=\s*('+NUMBER+r').*?Max\(alpha.water\)\s*=\s*('+NUMBER+')',last);assert rows
        lo,hi=map(float,rows[-1]);assert np.isfinite([*res.values(),lo,hi]).all()
        history.append(dict(step=len(history)+1,time_s=t,residuals=res,alpha_min=lo,alpha_max=hi,passed=all(res[k]<=targets[k] for k in targets) and lo>=-1e-6 and hi<=1+1e-6))
    assert len(history)==steps and np.allclose([r['time_s'] for r in history],np.arange(1,steps+1)*DT,atol=1e-12,rtol=0)
    courant=[float(v) for v in re.findall(r'^Interface Courant Number mean: '+NUMBER+r' max: ('+NUMBER+')',text,re.M)];assert len(courant)==steps+1 and np.isfinite(courant).all()
    g=geometry();order=g['openfoam_to_common'];c=native_field(run/'0/C',3);vol=native_field(run/'0/Vc',1).ravel();area=vol/DEPTH
    assert np.max(abs(c[:,:2]-g['centres'][order]))<1e-12 and np.max(abs(area-g['areas'][order]))<1e-15
    initial=initial_fields(g,args.mode);initial_alpha=initial['alpha'][order];M0=initial['mixture_mass_per_m']
    fluxpath=run/'postProcessing/atmosphereFlux/0/surfaceFieldValue.dat';fluxtext=fluxpath.read_text();assert 'orientedSum(rhoPhi)' in fluxtext and 'orientedSum(phi)' in fluxtext and 'patch(atmosphere)' in fluxtext
    assert int(re.search(r'# Faces\s*:\s*(\d+)',fluxtext)[1])==92
    assert abs(float(re.search(r'# Area\s*:\s*('+NUMBER+')',fluxtext)[1])-L*DEPTH)<1e-12
    flux=np.loadtxt(fluxpath,comments='#',ndmin=2);assert flux.shape==(steps+1,3) and np.isfinite(flux).all() and np.allclose(flux[:,0],np.arange(steps+1)*DT,atol=1e-12,rtol=0)
    assert np.max(abs(flux[0,1:]))<1e-12
    # Native patch flux is outward-positive and includes the actual extrusion depth.
    qmass=flux[:,1]/DEPTH;qvolume=flux[:,2]/DEPTH
    qwater=(qmass-MATERIALS['air'][0]*qvolume)/(MATERIALS['water-liquid'][0]-MATERIALS['air'][0])
    integrated_mass=np.r_[0,np.cumsum(qmass[1:])*DT];integrated_water=np.r_[0,np.cumsum(qwater[1:])*DT]
    records=[];sources=[run/name for name in ['solver.log','mesh-check.log','0/C','0/Vc','declared-settings.json']]+[fluxpath,Path(__file__)]
    for step in range(0,steps+1,interval):
        folder=run/f'{step*DT:.12g}';a=native_field(folder/'alpha.water',1).ravel();u=native_field(folder/'U',3);pr=native_field(folder/'p_rgh',1).ravel()
        rho=MATERIALS['air'][0]+a*(MATERIALS['water-liquid'][0]-MATERIALS['air'][0]);mass=float(rho@area);water=float(a@area);speed=float(np.linalg.norm(u,axis=1).max());change=float(abs(a-initial_alpha).max())
        if step:
            clock=folder/'uniform/time';s=clock.read_text();value=re.search(r'^value\s+('+NUMBER+r')\s*;',s,re.M);idx=re.search(r'^index\s+(\d+)\s*;',s,re.M);dt=re.search(r'^deltaT\s+('+NUMBER+r')\s*;',s,re.M)
            assert value and idx and dt and abs(float(value[1])-step*DT)<1e-12 and int(idx[1])==step and abs(float(dt[1])-DT)<1e-14;sources.append(clock)
        else:
            assert change<1e-12 and abs(water-WATER_AREA)<1e-12
            assert np.max(abs(pr-rho*G*c[:,1]-initial['pressure'][order]))<1e-9
            assert speed==0
        closure=(mass-M0+integrated_mass[step])/M0;water_closure=(water-WATER_AREA+integrated_water[step])/WATER_AREA
        bounds=a.min()>=-1e-6 and a.max()<=1+1e-6;quiet=args.mode!='flat' or (speed<=1e-4 and change<=.01)
        records.append(dict(step=step,time_s=step*DT,water_area_m2=water,mixture_mass_per_m=mass,integrated_outward_mass_per_m=float(integrated_mass[step]),relative_mixture_mass_closure=float(closure),relative_water_flux_closure=float(water_closure),max_speed_m_s=speed,maximum_alpha_change=change,alpha_min=float(a.min()),alpha_max=float(a.max()),fields_pass=bool(bounds and quiet and abs(closure)<=.005)))
        sources.extend(folder/name for name in ['alpha.water','U','p_rgh'])
    base=P/'solver-inputs-v004.json';proof=json.loads(base.read_text())['source_sha256'].copy();sources.append(base)
    for source in sources:proof[str(source)]=sha(source)
    passed=all(r['passed'] for r in history) and all(r['fields_pass'] for r in records) and max(courant)<=.25
    result=dict(solver='OpenFOAM14',mode=args.mode,cells=NCELLS,steps=steps,dt_s=DT,saved_states=len(records),final_time_s=steps*DT,numerical_checks_pass=bool(passed),reported_interface_courant_count=len(courant),expected_courant_reports=steps+1,peak_interface_courant=max(courant),native_flux_rows=len(flux),maximum_relative_mixture_mass_closure=max(abs(r['relative_mixture_mass_closure']) for r in records),maximum_relative_water_flux_closure=max(abs(r['relative_water_flux_closure']) for r in records),maximum_saved_speed_m_s=max(r['max_speed_m_s'] for r in records),mass_flux_sign='Native orientedSum rhoPhi outward-positive; divided by extrusion depth',mass_integration='Euler: end-of-step native flux, every actual step',records=records,step_records=history,native_source_sha256=proof,physical_validation=False,production_complete=args.mode=='production' and bool(passed))
    args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['records','step_records','native_source_sha256']},indent=2))
    if not passed:raise SystemExit(2)
if __name__=='__main__':main()