"""Continue a verified 5000-step prefix; never relabel the failed full parent."""
from pathlib import Path
import datetime,hashlib,json,re,shutil,subprocess,sys,traceback
import numpy as np
from dam_break_common_v002 import *
from dam_break_admission_v004 import verify,wait_resources
R=Path('E:/CFD-Learning-Lab');P=R/'local/dam-break-over-obstacle-v001';N=R/'native/runs/dam-break-over-obstacle-v001';OLD=N/'production/openfoam-dt005ms-v001';OUT=N/'production/openfoam-scheduled-urelax07-v001';STATUS=P/'scheduled-openfoam-v001-status.json';OWNED=False
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mark(state,**kw):
 tmp=STATUS.with_suffix('.tmp');tmp.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),**kw),indent=2)+'\n');tmp.replace(STATUS)
def prove():
 for name,h in json.loads((P/'openfoam-scheduled-inputs-v001.json').read_text())['source_sha256'].items():assert sha(Path(name))==h,name
def prepare():
 parent=json.loads((P/'openfoam-production-dt005ms-audit.json').read_text());assert not parent['numerical_checks_pass'];assert len(parent['step_records'])==10000 and all(r['passed'] for r in parent['step_records'][:5000])
 fa=json.loads((P/'fluent-production-dt005ms-audit.json').read_text());assert fa['numerical_checks_pass'] and fa['production_complete']
 window=N/'diagnostics/openfoam-pressure-cycle-urelax07-v002/window-audit.json';w=json.loads(window.read_text());assert w['window_checks_pass']
 for a in [parent,fa,w]:
  for name,h in a['native_source_sha256'].items():assert sha(Path(name))==h,name
 assert not OUT.exists();OUT.mkdir(parents=True);files=[Path(__file__),P/'verify_scheduled_continuation_v001.py',P/'run_scheduled_openfoam_v001.sh',P/'openfoam-production-dt005ms-audit.json',P/'fluent-production-dt005ms-audit.json',window]
 for folder in ['constant','system']:shutil.copytree(OLD/folder,OUT/folder)
 for step in range(0,5001,100):
  name=f'{step*DT:.12g}';shutil.copytree(OLD/name,OUT/name)
  for p in (OUT/name).rglob('*'):
   if p.is_file():assert sha(p)==sha(OLD/name/p.relative_to(OUT/name))
 for proc in ['processor0','processor1']:
  for folder in ['constant','0.25']:shutil.copytree(OLD/proc/folder,OUT/proc/folder)
 c=OUT/'system/controlDict';s=c.read_text();original=s
 s,n=re.subn(r'\bstartTime\s+[^;]+;','startTime 0.25;',s,count=1);assert n==1 and s.split('functions',1)[1]==original.split('functions',1)[1]
 with c.open('w',newline='\n') as f:f.write(s)
 fv=OUT/'system/fvSolution';s=fv.read_text();s,n=re.subn(r'relaxationFactors\s*\{\s*equations\s*\{\s*"\.\*"\s+1;\s*\}\s*\}', 'relaxationFactors { equations { U 0.7; UFinal 0.7; } }',s);assert n==1
 with fv.open('w',newline='\n') as f:f.write(s)
 d=json.loads((OLD/'declared-settings.json').read_text());d.update(momentum_relaxation_schedule=[dict(start_s=0,end_s=.25,value=1),dict(start_s=.25,end_s=.5,value=.7)],prefix_steps=5000,continued_steps=5000,parent_full_audit_retained_as_MISS=True)
 (OUT/'declared-settings.json').write_text(json.dumps(d,indent=2)+'\n')
 files += [p for p in OUT.rglob('*') if p.is_file() and ('uniform' not in p.parts or p.name=='time')]
 (P/'openfoam-scheduled-inputs-v001.json').write_text(json.dumps(dict(source_sha256={str(p):sha(p) for p in files},damping_window_audit=str(window),prefix_parent=str(OLD),prefix_steps=5000,suffix_steps=5000,relaxation_schedule_predeclared=True,thresholds_unchanged=True),indent=2)+'\n')
def audit():
 parent=json.loads((P/'openfoam-production-dt005ms-audit.json').read_text());history=parent['step_records'][:5000].copy()
 text=(OUT/'solver.log').read_text();assert re.search(r'^End\s*$',text,re.M) and 'FOAM FATAL' not in text
 chunks=re.split(r'^Time = ('+NUMBER+r')s?\s*$',text,flags=re.M)
 for i in range(1,len(chunks),2):
  t=float(chunks[i]);last=chunks[i+1].rsplit('PIMPLE: Iteration 8',1);assert len(last)==2;res={}
  for m in re.finditer(r'Solving for (\w+), Initial residual = ('+NUMBER+r'), Final residual = ('+NUMBER+')',last[1]):
   if m[1] in ['Ux','Uy','p_rgh']:res[m[1]]=max(res.get(m[1],0),float(m[2]))
  assert set(res)=={'Ux','Uy','p_rgh'}
  bounds=re.findall(r'Min\(alpha.water\)\s*=\s*('+NUMBER+r').*?Max\(alpha.water\)\s*=\s*('+NUMBER+')',last[1]);assert bounds;lo,hi=map(float,bounds[-1])
  assert np.isfinite([*res.values(),lo,hi]).all()
  history.append(dict(step=len(history)+1,time_s=t,residuals=res,alpha_min=lo,alpha_max=hi,passed=max(res.values())<=1e-4 and lo>=-1e-6 and hi<=1+1e-6))
 assert len(history)==10000 and np.allclose([r['time_s'] for r in history],np.arange(1,10001)*DT,atol=1e-12,rtol=0)
 def getco(t):return [float(v) for v in re.findall(r'^Interface Courant Number mean: '+NUMBER+r' max: ('+NUMBER+')',t,re.M)]
 pc=getco((OLD/'solver.log').read_text())[:5001];sc=getco(text);assert len(pc)==len(sc)==5001
 co_records=[dict(time_s=k*DT,value=v,segment='prefix') for k,v in enumerate(pc)]+[dict(time_s=.25+k*DT,value=v,segment='continuation') for k,v in enumerate(sc)]
 oldflux=OLD/'postProcessing/atmosphereFlux/0/surfaceFieldValue.dat';newflux=OUT/'postProcessing/atmosphereFlux/0.25/surfaceFieldValue.dat'
 pf=np.loadtxt(oldflux,comments='#',ndmin=2)[:5001];sf=np.loadtxt(newflux,comments='#',ndmin=2)
 assert pf.shape==sf.shape==(5001,3) and np.allclose(sf[:,0],.25+np.arange(5001)*DT,atol=1e-12,rtol=0)
 flux=np.concatenate([pf,sf[1:]]);assert flux.shape==(10001,3) and np.isfinite(flux).all() and np.allclose(flux[:,0],np.arange(10001)*DT,atol=1e-12,rtol=0)
 for p in [oldflux,newflux]:
  ft=p.read_text();assert 'patch(atmosphere)' in ft and 'orientedSum(rhoPhi)' in ft and 'orientedSum(phi)' in ft
  assert int(re.search(r'# Faces\s*:\s*(\d+)',ft)[1])==92 and abs(float(re.search(r'# Area\s*:\s*('+NUMBER+')',ft)[1])-L*DEPTH)<1e-12
 g=geometry();order=g['openfoam_to_common'];c=native_field(OUT/'0/C',3);areas=native_field(OUT/'0/Vc',1).ravel()/DEPTH
 assert np.max(abs(c[:,:2]-g['centres'][order]))<1e-12 and np.max(abs(areas-g['areas'][order]))<1e-15
 mass_integral=np.r_[0,np.cumsum(flux[1:,1]/DEPTH)*DT];water_integral=np.r_[0,np.cumsum((flux[1:,1]-flux[1:,2])/(999*DEPTH))*DT]
 expected=initial_fields(g,'production');M0=expected['mixture_mass_per_m'];records=[]
 sources=[OUT/'solver.log',OUT/'mesh-check.log',OUT/'0/C',OUT/'0/Vc',OUT/'declared-settings.json',oldflux,newflux,P/'openfoam-production-dt005ms-audit.json',P/'openfoam-scheduled-inputs-v001.json',Path(__file__)]
 for step in range(0,10001,100):
  folder=OUT/f'{step*DT:.12g}';a=native_field(folder/'alpha.water',1).ravel();u=native_field(folder/'U',3);pr=native_field(folder/'p_rgh',1).ravel()
  mass=float((1+999*a)@areas);water=float(a@areas);speed=float(np.linalg.norm(u,axis=1).max())
  if step:
   clock=folder/'uniform/time';st=clock.read_text();idx=int(re.search(r'^index\s+(\d+)\s*;',st,re.M)[1]);dt=float(re.search(r'^deltaT\s+('+NUMBER+r')\s*;',st,re.M)[1]);value=float(re.search(r'^value\s+('+NUMBER+r')\s*;',st,re.M)[1]);assert idx==step and abs(dt-DT)<1e-14 and abs(value-step*DT)<1e-12;sources.append(clock)
  else:
   assert np.max(abs(a-expected['alpha'][order]))<1e-12 and speed==0 and np.max(abs(pr-(1+999*a)*G*c[:,1]-expected['pressure'][order]))<1e-9
  closure=(mass-M0+mass_integral[step])/M0
  records.append(dict(step=step,time_s=step*DT,water_area_m2=water,mixture_mass_per_m=mass,integrated_outward_mass_per_m=float(mass_integral[step]),relative_mixture_mass_closure=float(closure),relative_water_flux_closure=float((water-WATER_AREA+water_integral[step])/WATER_AREA),max_speed_m_s=speed,maximum_alpha_change=float(abs(a-expected['alpha'][order]).max()),alpha_min=float(a.min()),alpha_max=float(a.max()),fields_pass=bool(a.min()>=-1e-6 and a.max()<=1+1e-6 and abs(closure)<=.005)))
  sources.extend(folder/n for n in ['alpha.water','U','p_rgh'])
 proof=parent['native_source_sha256'].copy()
 for n,h in json.loads((P/'openfoam-scheduled-inputs-v001.json').read_text())['source_sha256'].items():assert sha(Path(n))==h,n;proof[n]=h
 for p in sources:proof[str(p)]=sha(p)
 passed=all(r['passed'] for r in history) and all(r['fields_pass'] for r in records) and max(pc+sc)<=.25
 result=dict(solver='OpenFOAM14',mode='production',cells=NCELLS,steps=10000,dt_s=DT,saved_states=101,final_time_s=.5,numerical_checks_pass=bool(passed),production_complete=bool(passed),reported_interface_courant_count=10002,expected_courant_reports=10002,unique_physical_time_courant_samples=10001,peak_interface_courant=max(pc+sc),native_flux_rows=10002,unique_flux_time_rows=10001,maximum_relative_mixture_mass_closure=max(abs(r['relative_mixture_mass_closure']) for r in records),maximum_relative_water_flux_closure=max(abs(r['relative_water_flux_closure']) for r in records),maximum_saved_speed_m_s=max(r['max_speed_m_s'] for r in records),mass_flux_sign='Native orientedSum outward-positive; divided by extrusion depth',mass_integration='Euler native end-step flux; duplicate restart initial row not integrated',records=records,step_records=history,courant_records=co_records,native_source_sha256=proof,continued_calculation=True,preserved_prefix_steps=5000,newly_computed_steps=5000,relaxation_schedule=[dict(until_s=.25,value=1),dict(after_s=.25,value=.7)],parent_full_miss_preserved=True,physical_validation=False)
 path=P/'openfoam-production-scheduled-audit.json';assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n');return path,result
def main():
 global OWNED
 assert not STATUS.exists();OWNED=True;mark('preparing_verified_continuation');wait_resources('openfoam-production');verify('openfoam-production');prepare()
 resources=wait_resources('openfoam-production');verify('openfoam-production');prove();mark('running_continuation',resources=resources,preserved_prefix_steps=5000,new_steps=5000)
 with (P/'scheduled-openfoam-v001-wrapper.log').open('x',encoding='utf-8') as f:r=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/dam-break-over-obstacle-v001/run_scheduled_openfoam_v001.sh'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert r.returncode==0,(r.returncode,'scheduled-openfoam-v001-wrapper.log');mark('auditing_continued_calculation');path,a=audit();assert a['numerical_checks_pass'],'Continued case misses unchanged numerical criteria'
 f_a=P/'fluent-production-dt005ms-audit.json';fa=json.loads(f_a.read_text());assert fa['production_complete'] and fa['numerical_checks_pass']
 selection=dict(fluent_case=str(N/'production/fluent-coupled-dt005ms-v001'),openfoam_case=str(OUT),fluent_audit=str(f_a),openfoam_audit=str(path),audit_sha256={str(p):sha(p) for p in [f_a,path]},full_pair_numerical_checks_pass=True,dt_s=DT,steps_each=10000,states_each=101,openfoam_continued_from_verified_prefix=True,physical_validation=False)
 out=P/'accepted-full-pair-selection-v001.json';assert not out.exists();out.write_text(json.dumps(selection,indent=2)+'\n');mark('accepted_pair_selected_media_pending',selection=str(out))
if __name__=='__main__':
 try:main()
 except Exception:
  if OWNED:mark('failed',error=traceback.format_exc())
  raise