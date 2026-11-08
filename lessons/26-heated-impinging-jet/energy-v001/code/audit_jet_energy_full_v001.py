"""Complete actual600-step merged native energy audit. Never launches solver."""
import re,json
import numpy as np,h5py
from jet_energy_common_v001 import *
from jet_energy_ledger_v001 import ledger,fluent_heat300,foam_patch_heat300,fluent_steps,NUM
from jet_energy_restart_v001 import merge_steps,same_restart
from audit_jet_energy_qualification_v002 import clock,report,sources
from analyse_heated_jet_v002 import values,partial,counted,field,patch_values,common_mapping,TDIM,PDIM,N as NC,NF
from jet_common_v001 import ZONES,RHO,DX,DEPTH,NX,NY
CP=1006.;K=.0262

def locate(code,step):
 run=N/(code+('-qualification-v001' if step<=20 else '-production-v001'))
 return run/(f'heat-{step:04d}.dat.h5' if code=='fluent' else f'{step*.005:g}/T')
def accepted_prefix(code):
 p=N/(code+'-qualification-v001')/'qualification-audit-v002.json';a=json.loads(p.read_text());assert a['numerical_pass']
 for f,h in a['source_sha256'].items():assert sha(f)==h,f
 return a,p

def audit_fluent():
 o=N/'fluent-production-v001';q=N/'fluent-qualification-v001';a,ap=accepted_prefix('fluent');assert json.loads((o/'status.json').read_text())['stage']=='closed'
 assert sha(o/'source.cas.h5')==sha(q/'qualification-final.cas.h5') and sha(o/'source-0020.dat.h5')==sha(q/'qualification-final.dat.h5')
 clock(o/'source-0020.dat.h5',20);clock(o/'heat-0020.dat.h5',20)
 defs=[('SV_T','cells',NC),('SV_H','cells',NC),('SV_U','cells',NC),('SV_V','cells',NC),('SV_P','cells',NC),('SV_FLUX','faces',NF)]
 same_restart({k:values(q/'heat-0020.dat.h5',k,b,n) for k,b,n in defs},{k:values(o/'heat-0020.dat.h5',k,b,n) for k,b,n in defs})
 with np.load(PARENT/'geometry.npz') as g:
  assert np.array_equal(common_mapping(o/'active.cas.h5',g),common_mapping(q/'active.cas.h5',g))
 with h5py.File(o/'source.cas.h5') as h:
  z=h['meshes/1/faces/zoneTopology'];zones={n:slice(int(lo)-1,int(hi)) for n,lo,hi in zip(z['name'][0].decode().split(';'),z['minId'][...],z['maxId'][...]) if n in ZONES}
 assert set(zones)==set(ZONES)
 h0=values(q/'source-0000.dat.h5','SV_H');assert np.ptp(h0)<1e-10;href=float(h0[0]);baseline={k:values(q/'source-0000.dat.h5',k,b,n) for k,b,n in defs if k not in ['SV_T','SV_H']}
 E=[];Q={n:[] for n in ZONES};mass=[];bounds=[];frozen=0.;replay={};paths=[ap,o/'source.cas.h5',o/'source-0020.dat.h5',o/'active.cas.h5',o/'heat-0020.dat.h5',o/'solver-transcript.txt',o/'loaded-fixed-settings.json',o/'fixed-settings-comparison.json',o/'run-settings.json'];states=[]
 for s in range(601):
  p=locate('fluent',s);paths.append(p);clock(p,s);T=values(p,'SV_T');h=values(p,'SV_H');assert np.isfinite(T).all() and np.isfinite(h).all();assert np.max(abs(h-href-CP*(T-300)))<=1e-8
  bounds.append([float(T.min()),float(T.max())]);E.append(float(RHO*DX**2*np.sum(h-href)));states.append(dict(step=s,time_s=s*.005,path=str(p)))
  for k,b,n in defs:
   if k in baseline:
    v=values(p,k,b,n);assert np.isfinite(v).all();frozen=max(frozen,float(abs(v-baseline[k]).max()))
  if s%20==0 or s==599:
   original=FBASE/f'heat-{s:04d}.dat.h5';replay[str(s)]=float(abs(T-values(original,'SV_T')).max());paths.append(original)
  if not s:continue
  flux=values(p,'SV_FLUX','faces',NF)
  with h5py.File(p) as hf:heat=partial(hf['results/1/phase-1/faces/SV_HEAT_FLUX'],NF)
  raw={n:float(heat[zones[n]].sum()) for n in ZONES};mo={n:float(flux[zones[n]].sum()) for n in ZONES};assert np.isfinite([*raw.values(),*mo.values()]).all()
  run=p.parent;hr=run/f'heat-report-{s:04d}.txt';mr=run/f'mass-report-{s:04d}.txt';report(hr,raw);report(mr,{n:-v for n,v in mo.items()});paths += [hr,mr]
  for n,v in fluent_heat300(raw,mo,href).items():Q[n].append(v)
  mass.append(abs(sum(mo.values()))/abs(mo['inlet']))
 continued=fluent_steps((o/'solver-transcript.txt').read_text(),21,600,.005,30);rows=merge_steps(a['step_residuals'],continued)
 return finish('fluent',o,E,Q,mass,frozen,bounds,rows,replay,1e-8,paths,states)

def time_meta(p,s):
 text=p.read_text()
 for key,expect in [('value',s*.005),('index',s),('deltaT',.005)]:
  a=re.findall(r'^\s*'+key+r'\s+('+NUM+r')\s*;',text,re.M);assert len(a)==1 and abs(float(a[0])-expect)<=1e-10,(p,key,a,expect)

def audit_openfoam():
 o=N/'openfoam-production-v001';q=N/'openfoam-qualification-v001';a,ap=accepted_prefix('openfoam');assert (o/'run-status.txt').read_text().strip()=='complete'
 assert 'Mesh OK.' in (o/'mesh-check.log').read_text();mesh=o/'constant/polyMesh';npch,b=counted(mesh/'boundary');patch={n:dict(count=int(re.search(r'nFaces\s+(\d+)',v)[1]),start=int(re.search(r'startFace\s+(\d+)',v)[1])) for n,v in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',b)};assert set(patch)==set(ZONES)|{'frontAndBack'}
 num,b=counted(mesh/'owner');owner=np.fromstring(b,sep=' ',dtype=int);assert len(owner)==num;ni=min(v['start'] for v in patch.values())
 defs=[('U',3,[0,1,-1,0,0,0,0],NC),('p',1,[0,2,-2,0,0,0,0],NC),('phi',1,PDIM,ni)]
 same_restart({k:field(q/'0.1'/k,n,w,d) for k,w,d,n in defs}|{'T':field(q/'0.1/T',NC,1,TDIM)},{k:field(o/'0.1'/k,n,w,d) for k,w,d,n in defs}|{'T':field(o/'0.1/T',NC,1,TDIM)})
 time_meta(o/'0.1/uniform/time',20)
 ph={n:patch_values(o/'1250/phi',n,patch[n]['count']) for n in ZONES};mo={n:float(ph[n].sum()*RHO/DEPTH) for n in ZONES};massratio=abs(sum(mo.values()))/abs(mo['inlet']);baseline={k:field(o/'1250'/k,n,w,d) for k,w,d,n in defs}
 paths=[ap,o/'solver.log',o/'mesh-check.log']+[p for p in (o/'constant').rglob('*') if p.is_file()];E=[];Q={n:[] for n in ZONES};bounds=[];frozen=0.;replay={};states=[]
 for s in range(601):
  p=locate('openfoam',s);run=p.parent.parent;time=p.parent.name;T=field(p,NC,1,TDIM);paths.append(p);assert np.isfinite(T).all();bounds.append([float(T.min()),float(T.max())]);E.append(float(RHO*CP*DX**2*np.sum(T-300)));states.append(dict(step=s,time_s=s*.005,path=str(p)))
  if s%20==0 or s==599:
   original=OBASE/time/'T';replay[str(s)]=float(abs(T-field(original,NC,1,TDIM)).max());paths.append(original)
  for rank in [0,1]:
   rp=run/f'processor{rank}'/time/'T';assert rp.exists();paths.append(rp);txt=rp.read_text();loc=re.search(r'location\s+"([^"]+)"',txt);assert loc and abs(float(loc[1])-s*.005)<=1e-10
   if s:tm=rp.parent/'uniform/time';time_meta(tm,s);paths.append(tm)
  if not s:continue
  if s:tm=p.parent/'uniform/time';time_meta(tm,s);paths.append(tm)
  for n in ZONES:
   b=patch[n];owners=owner[b['start']:b['start']+b['count']];Q[n].append(foam_patch_heat300(n,T[owners],ph[n],dx=DX,depth=DEPTH,rho=RHO,cp=CP,k=K))
  for k,w,d,n in defs:
   fp=p.parent/k;paths.append(fp);v=field(fp,n,w,d);assert np.isfinite(v).all();frozen=max(frozen,float(abs(v-baseline[k]).max()))
  for n in ZONES:frozen=max(frozen,float(abs(patch_values(p.parent/'phi',n,patch[n]['count'])-ph[n]).max()))
 text=(o/'solver.log').read_text();times=[float(v) for v in re.findall(r'^Time = ('+NUM+r')s?\s*$',text,re.M)];res=[float(v) for v in re.findall(r'Solving for T, Initial residual = '+NUM+r', Final residual = ('+NUM+r')',text)];assert len(times)==len(res)==580
 continued=[dict(step=i,time_s=t,residual=v) for i,t,v in zip(range(21,601),times,res)];rows=merge_steps(a['step_residuals'],continued)
 return finish('openfoam',o,E,Q,[massratio]*600,frozen,bounds,rows,replay,1e-12,paths,states)

def finish(code,o,E,Q,mass,frozen,bounds,rows,replay,reslimit,paths,states):
 l=ledger(list(range(1,601)),[i*.005 for i in range(1,601)],E,Q,.005)
 checks=dict(allstep_and_prefix_energy=l['passed'],frozen_flow=frozen<=1e-12,mass=max(mass)<=1e-4,temperature=min(v[0] for v in bounds)>=299.99 and max(v[1] for v in bounds)<=320.01,residual=max(v['residual'] for v in rows)<=reslimit,all32_original_state_replays=len(replay)==32 and max(replay.values())<=1e-8)
 record=dict(code=code,numerical_pass=all(checks.values()),checks=checks,actual_steps=600,prefix_steps=20,continued_steps=580,unique_states=601,dt_s=.005,end_time_s=3.,continuation_not_independent_full_rerun=True,energy_ledger=l,patch_input_W_per_m=Q,initial_energy_J_per_m=E[0],max_mass_relative=max(mass),max_frozen_change=frozen,temperature_range_K=[min(v[0] for v in bounds),max(v[1] for v in bounds)],step_residuals=rows,thermal_residual_limit=reslimit,max_terminal_residual=max(v['residual'] for v in rows),original_state_T_max_abs_replay_difference_K=replay,native_states=states,source_sha256=sources(sorted(set(paths))),physical_validation=False,time_independence=False,heat_scope=('Native Fluent integrated face heat verified against all600 actual heat/mass reports; 300K datum shift' if code=='fluent' else 'Reconstruction from actual every-step T/frozen phi under unchanged upwind/orthogonal scheme; not native heat report'))
 dump(o/'production-audit-v001.json',record);return record
if __name__=='__main__':
 import sys
 a=audit_fluent() if sys.argv[1]=='fluent' else audit_openfoam();print(json.dumps({k:a[k] for k in ['code','numerical_pass','checks','actual_steps','unique_states','max_terminal_residual']}));assert a['numerical_pass']
