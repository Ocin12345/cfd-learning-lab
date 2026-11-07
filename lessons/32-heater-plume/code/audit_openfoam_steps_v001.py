"""Audit completed native OpenFOAM coupling loops without equating its normalization to Fluent's."""
from pathlib import Path
import re,json,sys
p=Path(sys.argv[1]);outer=int(sys.argv[2]);expected=int(sys.argv[3]);out=Path(sys.argv[4])
raw=(p/'solver.log').read_text();assert re.search(r'^End\s*$',raw,re.M)
blocks=re.split(r'^Time = ([0-9.eE+-]+)s?\s*$',raw,flags=re.M);lim={'Ux':1e-6,'Uy':1e-6,'h':1e-8,'p_rgh':1e-6};rows=[]
for i in range(1,len(blocks),2):
 t=float(blocks[i]);chunk=blocks[i+1];assert f'PIMPLE: Iteration {outer}' in chunk
 last=chunk.split(f'PIMPLE: Iteration {outer}')[-1];vals={}
 for m in re.finditer(r'Solving for (\w+), Initial residual = ([-+\d.eE]+), Final residual = ([-+\d.eE]+)',last):
  name=m[1]
  if name in lim:vals[name]=max(vals.get(name,0),float(m[2]))
 assert set(vals)==set(lim)
 rows.append(dict(time_s=t,residuals=vals,passed=all(v<=lim[k] for k,v in vals.items())))
a=dict(expected_steps=expected,steps_found=len(rows),targets=lim,outer_correctors=outer,all_passed=len(rows)==expected and all(r['passed'] for r in rows),maxima={k:max(r['residuals'][k] for r in rows) for k in lim},misses=[r for r in rows if not r['passed']],steps=rows)
out.write_text(json.dumps(a,indent=2)+'\n');print(json.dumps({k:v for k,v in a.items() if k not in ['steps','misses']},indent=2));print('misses',len(a['misses']));assert a['all_passed']
