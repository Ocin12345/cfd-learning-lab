"""Same-checkpoint 100-step control and momentum-relaxed diagnostic."""
from pathlib import Path
import datetime,hashlib,json,re,shutil,subprocess,sys,traceback
import numpy as np
from dam_break_common_v002 import *
from dam_break_admission_v004 import verify,wait_resources
R=Path('E:/CFD-Learning-Lab');P=R/'local/dam-break-over-obstacle-v001';N=R/'native/runs/dam-break-over-obstacle-v001';OLD=N/'production/openfoam-dt005ms-v001';STATUS=P/'openfoam-pressure-cycle-v002-status.json';OWNED=False
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mark(state,**kw):
 tmp=STATUS.with_suffix('.tmp');tmp.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),**kw),indent=2)+'\n');tmp.replace(STATUS)
def prove_sources():
 a=json.loads((P/'openfoam-pressure-cycle-inputs-v002.json').read_text())
 for name,h in a['source_sha256'].items():assert sha(Path(name))==h,name
def audit(case,label):
 text=(case/'solver.log').read_text();assert re.search(r'^End\s*$',text,re.M) and 'FOAM FATAL' not in text
 parts=re.split(r'^Time = ('+NUMBER+r')s?\s*$',text,flags=re.M);rows=[]
 for i in range(1,len(parts),2):
  t=float(parts[i]);chunk=parts[i+1];last=chunk.rsplit('PIMPLE: Iteration 8',1);assert len(last)==2;res={}
  for m in re.finditer(r'Solving for (\w+), Initial residual = ('+NUMBER+r'), Final residual = ('+NUMBER+')',last[1]):
   if m[1] in ['Ux','Uy','p_rgh']:res[m[1]]=max(res.get(m[1],0),float(m[2]))
  assert set(res)=={'Ux','Uy','p_rgh'}
  bounds=re.findall(r'Min\(alpha.water\)\s*=\s*('+NUMBER+r').*?Max\(alpha.water\)\s*=\s*('+NUMBER+')',last[1]);assert bounds;lo,hi=map(float,bounds[-1]);assert np.isfinite([*res.values(),lo,hi]).all()
  rows.append(dict(time_s=t,residuals=res,alpha_min=lo,alpha_max=hi,passed=max(res.values())<=1e-4 and lo>=-1e-6 and hi<=1+1e-6))
 assert len(rows)==100 and np.allclose([r['time_s'] for r in rows],.255+np.arange(1,101)*DT,atol=1e-12,rtol=0)
 co=[float(v) for v in re.findall(r'^Interface Courant Number mean: '+NUMBER+r' max: ('+NUMBER+')',text,re.M)];assert len(co)==101
 fluxpath=case/'postProcessing/atmosphereFlux/0.255/surfaceFieldValue.dat';flux=np.loadtxt(fluxpath,comments='#',ndmin=2);assert flux.shape==(101,3) and np.allclose(flux[:,0],.255+np.arange(101)*DT,atol=1e-12,rtol=0)
 g=geometry();om=g['openfoam_to_common'];areas=g['areas'][om];q=flux[:,1]/DEPTH;integrated=np.r_[0,np.cumsum(q[1:])*DT];field_records=[];M0=None;sources=[case/'solver.log',case/'mesh-check.log',fluxpath,Path(__file__),case/'declared-diagnostic.json']
 for k in range(0,101,20):
  folder=case/f'{.255+k*DT:.12g}';a=native_field(folder/'alpha.water',1).ravel();u=native_field(folder/'U',3);pr=native_field(folder/'p_rgh',1).ravel();mass=float((1+999*a)@areas)
  if M0 is None:M0=mass
  clock=folder/'uniform/time';s=clock.read_text();idx=int(re.search(r'^index\s+(\d+)\s*;',s,re.M)[1]);dt=float(re.search(r'^deltaT\s+('+NUMBER+r')\s*;',s,re.M)[1]);value=float(re.search(r'^value\s+('+NUMBER+r')\s*;',s,re.M)[1]);assert idx==5100+k and abs(dt-DT)<1e-14 and abs(value-(.255+k*DT))<1e-12
  closure=(mass-M0+integrated[k])/M0
  field_records.append(dict(time_s=value,alpha_min=float(a.min()),alpha_max=float(a.max()),relative_mass_closure=float(closure),fields_pass=bool(a.min()>=-1e-6 and a.max()<=1+1e-6 and abs(closure)<=.005)))
  sources += [folder/n for n in ['alpha.water','U','p_rgh','uniform/time']]
 result=dict(label=label,diagnostic_only=True,steps=100,start_time_s=.255,end_time_s=.26,dt_s=DT,outer_correctors=8,momentum_equation_relaxation=1 if label=='baseline' else .7,missed_steps=[r for r in rows if not r['passed']],peak_interface_courant=max(co),maximum_relative_mass_closure=max(abs(r['relative_mass_closure']) for r in field_records),window_checks_pass=bool(all(r['passed'] for r in rows) and max(co)<=.25 and all(r['fields_pass'] for r in field_records)),step_records=rows,field_records=field_records,native_source_sha256={str(p):sha(p) for p in sources},full_production_complete=False)
 out=case/'window-audit.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');return result
def main():
 global OWNED
 assert not STATUS.exists();OWNED=True;mark('preparing')
 for name,expected in [('fluent-production-dt005ms-audit.json',True),('openfoam-production-dt005ms-audit.json',False)]:
  a=json.loads((P/name).read_text());assert a['numerical_checks_pass']==expected
  for source,h in a['native_source_sha256'].items():assert sha(Path(source))==h,source
 wait_resources('openfoam-production');verify('openfoam-production')
 sources=[Path(__file__),P/'verify_pressure_cycle_recipe_v002.py',P/'run_openfoam_pressure_cycle_v002.sh',P/'fluent-production-dt005ms-audit.json',P/'openfoam-production-dt005ms-audit.json'];cases=[]
 for label in ['baseline','urelax07']:
  case=N/'diagnostics'/('openfoam-pressure-cycle-'+label+'-v002');assert not case.exists();case.mkdir(parents=True);cases.append((label,case))
  for folder in ['0','constant','system','0.255']:shutil.copytree(OLD/folder,case/folder)
  for proc in ['processor0','processor1']:
   for folder in ['constant','0.255']:shutil.copytree(OLD/proc/folder,case/proc/folder)
  c=case/'system/controlDict';s=c.read_text()
  for key,val in [('startTime','0.255'),('endTime','0.26'),('writeInterval','20')]:s,n=re.subn(r'\b'+key+r'\s+[^;]+;',key+' '+val+';',s,count=1);assert n==1,key
  assert s.split('functions',1)[1]==c.read_text().split('functions',1)[1],'Native per-step flux monitor changed'
  with c.open('w',newline='\n') as f:f.write(s)
  fv=case/'system/fvSolution';s=fv.read_text()
  if label=='urelax07':
   s,n=re.subn(r'relaxationFactors\s*\{\s*equations\s*\{\s*"\.\*"\s+1;\s*\}\s*\}', 'relaxationFactors { equations { U 0.7; UFinal 0.7; } }',s);assert n==1
   with fv.open('w',newline='\n') as f:f.write(s)
  declare=dict(label=label,source_case=str(OLD),source_time_s=.255,source_global_step=5100,dt_s=DT,steps=100,end_time_s=.26,changed_candidate_control='U/UFinal momentum equation relaxation only, 1 to0.7',outer_correctors=8,source_fields='all stored partitioned fields at0.255 and unchanged partition meshes',physical_validation=False)
  (case/'declared-diagnostic.json').write_text(json.dumps(declare,indent=2)+'\n')
  for p in case.rglob('*'):
   if p.is_file() and ('uniform' not in p.parts or p.name=='time'):sources.append(p)
 recipe=P/'openfoam-pressure-cycle-inputs-v002.json';assert not recipe.exists();recipe.write_text(json.dumps(dict(source_sha256={str(p):sha(p) for p in sources},full_simulation=False),indent=2)+'\n')
 for label,case in cases:
  resources=wait_resources('openfoam-production');verify('openfoam-production');prove_sources();mark('running',variant=label,resources=resources)
  posix=case.as_posix().replace('E:/','/mnt/e/')
  with (case/'wrapper.log').open('x',encoding='utf-8') as f:r=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/dam-break-over-obstacle-v001/run_openfoam_pressure_cycle_v002.sh',posix],stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,creationflags=subprocess.CREATE_NO_WINDOW)
  assert r.returncode==0,(r.returncode,str(case/'wrapper.log'));a=audit(case,label);mark('audited',variant=label,passed=a['window_checks_pass'],missed_steps=len(a['missed_steps']),peak_courant=a['peak_interface_courant'])
 mark('both_variants_computed_and_audited',full_production_complete=False)
if __name__=='__main__':
 try:main()
 except Exception:
  if OWNED:mark('failed',error=traceback.format_exc())
  raise