import csv,re,json
import numpy as np,h5py
from jet_energy_common_v001 import *
from jet_energy_ledger_v001 import ledger,fluent_heat300,foam_patch_heat300,fluent_steps,NUM
from analyse_heated_jet_v002 import values,partial,counted,field,patch_values,common_mapping,TDIM,PDIM,N as NC,NF
from jet_common_v001 import ZONES,RHO,DX,DEPTH

CP=1006.;K=.0262

def clock(p,step):
 from ansys.fluent.core.offline import _lispy
 with h5py.File(p) as h:tree=_lispy.parse(h['settings/Data Variables'][0].decode())
 current={v[0]:v[1] for v in tree[1]}
 for name,expected in [('flow-time',step*.005),('time-step',step),('physical-time-step',.005)]:
  value=float(current[name]);assert np.isfinite(value) and abs(value-expected)<=1e-10,(name,value,expected)
def report(path,expected):
 t=path.read_text()
 for n,v in expected.items():
  a=re.findall(r'^\s*'+n+r'\s+('+NUM+r')\s*$',t,re.M)
  assert len(a)==1 and np.isclose(float(a[0]),v,rtol=1e-7,atol=1e-11),(path,n,a,v)
def sources(paths):return {str(p):sha(p) for p in paths}
def audit_fluent():
 o=N/'fluent-qualification-v001';assert json.loads((o/'status.json').read_text())['stage']=='closed'
 paths=[o/'active.cas.h5',o/'source.cas.h5',o/'source-0000.dat.h5',o/'solver-transcript.txt',o/'loaded-fixed-settings.json',o/'fixed-settings-comparison.json']
 assert sha(o/'source.cas.h5')==sha(FBASE/'heated-jet.cas.h5') and sha(o/'source-0000.dat.h5')==sha(FBASE/'heat-0000.dat.h5')
 with np.load(PARENT/'geometry.npz') as geo:assert np.array_equal(common_mapping(o/'active.cas.h5',geo),common_mapping(o/'source.cas.h5',geo))
 with h5py.File(o/'source.cas.h5') as h:
  z=h['meshes/1/faces/zoneTopology'];zones={n:slice(int(a)-1,int(b)) for n,a,b in zip(z['name'][0].decode().split(';'),z['minId'][...],z['maxId'][...]) if n in ZONES}
 assert set(zones)==set(ZONES)
 h300=values(o/'source-0000.dat.h5','SV_H');assert np.ptp(h300)<1e-10;href=float(h300[0]);E=[];Q={n:[] for n in ZONES};mass=[];frozen=0.;bounds=[]
 for step in range(21):
  p=o/f'heat-{step:04d}.dat.h5';paths.append(p);clock(p,step);T=values(p,'SV_T');h=values(p,'SV_H');assert np.isfinite(T).all() and np.isfinite(h).all();assert np.max(abs(h-href-CP*(T-300)))<=1e-8;bounds.append([float(T.min()),float(T.max())]);E.append(float(RHO*DX**2*np.sum(h-href)))
  for k,kind,n in [('SV_U','cells',NC),('SV_V','cells',NC),('SV_P','cells',NC),('SV_FLUX','faces',NF)]:frozen=max(frozen,float(np.max(abs(values(p,k,kind,n)-values(o/'source-0000.dat.h5',k,kind,n)))))
  if not step:continue
  flux=values(p,'SV_FLUX','faces',NF)
  with h5py.File(p) as h:heat=partial(h['results/1/phase-1/faces/SV_HEAT_FLUX'],NF)
  raw={n:float(heat[zones[n]].sum()) for n in ZONES};mo={n:float(flux[zones[n]].sum()) for n in ZONES};assert all(np.isfinite(list(raw.values())))
  a=o/f'heat-report-{step:04d}.txt';b=o/f'mass-report-{step:04d}.txt';report(a,raw);report(b,{n:-x for n,x in mo.items()});paths += [a,b]
  qi=fluent_heat300(raw,mo,href)
  for n in Q:Q[n].append(qi[n])
  mass.append(abs(sum(mo.values()))/abs(mo['inlet']))
 rows=fluent_steps((o/'solver-transcript.txt').read_text(),1,20,.005,30);replay=float(abs(values(o/'heat-0020.dat.h5','SV_T')-values(FBASE/'heat-0020.dat.h5','SV_T')).max())
 return finish('fluent',o,E,Q,mass,frozen,bounds,rows,replay,1e-8,paths)

def audit_openfoam():
 o=N/'openfoam-qualification-v001';assert (o/'run-status.txt').read_text().strip()=='complete'
 assert 'Mesh OK.' in (o/'mesh-check.log').read_text();mesh=o/'constant/polyMesh';npch,body=counted(mesh/'boundary');patch={n:dict(count=int(re.search(r'nFaces\s+(\d+)',b)[1]),start=int(re.search(r'startFace\s+(\d+)',b)[1])) for n,b in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',body)}
 assert set(patch)==set(ZONES)|{'frontAndBack'}
 n,b=counted(mesh/'owner');owner=np.fromstring(b,sep=' ',dtype=int);assert len(owner)==n;ni=min(p['start'] for p in patch.values())
 basephi=o/'1250/phi';ph={n:patch_values(basephi,n,patch[n]['count']) for n in ZONES};massout={n:float(ph[n].sum()*RHO/DEPTH) for n in ZONES};massratio=abs(sum(massout.values()))/abs(massout['inlet'])
 paths=[p for p in (o/'constant').rglob('*') if p.is_file()]+[o/'solver.log',o/'mesh-check.log',basephi];E=[];Q={n:[] for n in ZONES};bounds=[];frozen=0.
 for step in range(21):
  name=f'{step*.005:g}';p=o/name/'T';T=field(p,NC,1,TDIM);paths.append(p);assert np.isfinite(T).all();bounds.append([float(T.min()),float(T.max())]);E.append(float(RHO*CP*DX**2*np.sum(T-300)))
  for rank in [0,1]:
   rp=o/f'processor{rank}'/name/'T';assert rp.exists();paths.append(rp)
   # Reconstructed values plus both actual partition timestamps, no interpolated states.
   txt=rp.read_text();a=re.search(r'location\s+"([^"]+)"',txt);assert a and abs(float(a[1])-step*.005)<1e-10
  if not step:continue
  for n in ZONES:
   pch=patch[n];own=owner[pch['start']:pch['start']+pch['count']];Q[n].append(foam_patch_heat300(n,T[own],ph[n],dx=DX,depth=DEPTH,rho=RHO,cp=CP,k=K))
  for k,w,dim,n in [('U',3,[0,1,-1,0,0,0,0],NC),('p',1,[0,2,-2,0,0,0,0],NC),('phi',1,PDIM,ni)]:
   p=o/name/k;paths.append(p);frozen=max(frozen,float(abs(field(p,n,w,dim)-field(o/'1250'/k,n,w,dim)).max()))
  for n in ZONES:frozen=max(frozen,float(abs(patch_values(o/name/'phi',n,patch[n]['count'])-ph[n]).max()))
 text=(o/'solver.log').read_text();times=[float(x) for x in re.findall(r'^Time = ('+NUM+r')s?\s*$',text,re.M)];res=[float(x) for x in re.findall(r'Solving for T, Initial residual = '+NUM+r', Final residual = ('+NUM+r')',text)];assert len(times)==20 and len(res)==20 and np.allclose(times,np.arange(1,21)*.005,rtol=0,atol=1e-10)
 rows=[dict(step=s,time_s=t,residual=v) for s,t,v in zip(range(1,21),times,res)];replay=float(abs(field(o/'.1/T',NC,1,TDIM)-field(OBASE/'.1/T',NC,1,TDIM)).max()) if (o/'.1/T').exists() else float(abs(field(o/'0.1/T',NC,1,TDIM)-field(OBASE/'0.1/T',NC,1,TDIM)).max())
 return finish('openfoam',o,E,Q,[massratio]*20,frozen,bounds,rows,replay,1e-12,paths)

def finish(code,o,E,Q,mass,frozen,bounds,rows,replay,reslimit,paths):
 l=ledger(list(range(1,21)),[i*.005 for i in range(1,21)],E,Q,.005);checks=dict(energy_every_step_and_prefix=l['passed'],frozen_flow=frozen<=1e-12,mass=max(mass)<=1e-4,temperature=min(x[0] for x in bounds)>=299.99 and max(x[1] for x in bounds)<=320.01,residual=max(r['residual'] for r in rows)<=reslimit,output_only_replay=replay<=1e-8)
 a=dict(code=code,actual_steps=20,actual_states=21,dt_s=.005,end_time_s=.1,checks=checks,numerical_pass=all(checks.values()),initial_energy_J_per_m=E[0],energy_ledger=l,patch_input_W_per_m=Q,max_mass_relative=max(mass),max_frozen_change=frozen,max_replay_T_difference_K=replay,thermal_residual_limit=reslimit,max_terminal_residual=max(r['residual'] for r in rows),step_residuals=rows,source_sha256=sources(paths),heat_evidence=('Native Fluent face-integrated heat, verified against every step actual flux report and shifted to300K' if code=='fluent' else 'Reconstructed from native owner T and frozen phi using unchanged upwind/orthogonal scheme; not native heat report'))
 dump(o/'qualification-audit-v002.json',a);return a
if __name__=='__main__':
 import sys
 a=audit_fluent() if sys.argv[1]=='fluent' else audit_openfoam();print(json.dumps({k:a[k] for k in ['code','numerical_pass','checks','max_replay_T_difference_K','max_terminal_residual']}));assert a['numerical_pass']
