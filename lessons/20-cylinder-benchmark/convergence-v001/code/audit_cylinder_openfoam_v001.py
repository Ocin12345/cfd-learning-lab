from pathlib import Path
import csv,hashlib,json,math,re
from analyse_vortex_openfoam_v002 import residual_steps,native_boundary,boundary_mass,native_field,force_data,fatal_solver_error,FLOAT
from plot_vortex_shedding_v002 import openfoam_mesh,foam_labels

def verify_residual_window(rows,start,steps,enforce=True):
 if len(rows)!=steps:raise ValueError('Missing/extra physical residual steps')
 for k,r in enumerate(rows,1):
  if abs(r['time_s']-(start+k*.005))>1e-8 or not r['variables_complete'] or not r['all_printed_residuals_finite']:raise ValueError('Native residual time/coverage/finite mismatch')
  passed=all(r['variables'][n][metric]<=1e-5 for n in ['Ux','Uy','p'] for metric in ['terminal_initial','terminal_final'])
  if enforce and not passed:raise ValueError('Native terminal residual misses unchanged1e-5')
 return rows

def audit_openfoam(case,start=4.,steps=100):
 case=Path(case);text=(case/'solver.log').read_text();assert not fatal_solver_error(text) and re.search(r'^End\s*$',text,re.M)
 rows=verify_residual_window(residual_steps(text,'solver.log'),start,steps,enforce=False);misses=[r for r in rows if any(r['variables'][n][metric]>1e-5 for n in ['Ux','Uy','p'] for metric in ['terminal_initial','terminal_final'])]
 sections=re.split(r'^Time =\s*('+FLOAT+r')\s*s?\s*$',text,flags=re.M)
 for j in range(1,len(sections),2):assert list(map(int,re.findall(r'PIMPLE: Iteration\s+(\d+)',sections[j+1])))==list(range(1,9))
 control=(case/'system/controlDict').read_text();assert re.search(r'\bdeltaT\s+0.005\s*;',control) and re.search(r'\bwriteInterval\s+1\s*;',control)
 mesh=openfoam_mesh(case);assert mesh['cells']==21598 and abs(mesh['span_m']-1)<1e-10
 meshdir=next(p for p in [case/'constant/polyMesh',case/'0/polyMesh'] if p.exists());patches=native_boundary((meshdir/'boundary').read_text());ninternal=len(foam_labels(meshdir/'neighbour'))
 dirs={round(float(p.name),8):p for p in case.iterdir() if p.is_dir() and re.fullmatch(r'[0-9.]+',p.name) and float(p.name)>=start-1e-8};times=[round(start+k*.005,8) for k in range(steps+1)];assert set(dirs)==set(times)
 fields=[];inputfiles=[]
 for t in times:
  folder=dirs[t]
  for rank in [0,1]:
   clock=case/('processor'+str(rank))/folder.name/'uniform/time';raw=clock.read_text();m=re.search(r'\bvalue\s+('+FLOAT+r')\s*;',raw);assert m and abs(float(m[1])-t)<1e-8;inputfiles.append(clock)
  stats={n:native_field((folder/n).read_text(),21598,width,dims) for n,width,dims in [('U',3,[0,1,-1,0,0,0,0]),('p',1,[0,2,-2,0,0,0,0])]}
  raw=(folder/'phi').read_text();native_field(raw,ninternal,1,[0,3,-1,0,0,0,0]);mass=boundary_mass(raw,patches);assert mass['inlet_inflow_outlet_outflow'];fields.append(dict(time_s=t,fields=stats,mass=mass));inputfiles +=[folder/n for n in ['U','p','phi']]
 massmax=max(r['mass']['relative_all_physical_boundary_imbalance'] for r in fields)
 fp=list((case/'postProcessing/forceCoefficients').rglob('forceCoeffs.dat'));assert fp;force_rows=[]
 for p in fp:force_rows+=force_data(p.read_text());inputfiles.append(p)
 force_rows=sorted(force_rows);fmap={round(r[0],8):r for r in force_rows};assert len(fmap)==len(force_rows) and set(times[1:])<=set(fmap)
 with (case/'force-history.csv').open('x',newline='') as f:
  w=csv.writer(f);w.writerow(['time_s','Cd','Cl']);w.writerows(force_rows)
 co=[];state=start
 for line in text.splitlines():
  m=re.fullmatch(r'Time =\s*('+FLOAT+r')s?\s*',line)
  if m:state=float(m[1])
  m=re.search(r'Courant Number mean:\s*('+FLOAT+r')\s+max:\s*('+FLOAT+r')',line)
  if m:
   mean,maximum=map(float,m.groups());assert math.isfinite(mean) and math.isfinite(maximum);co.append(dict(native_pre_advance_state_time_s=state,mean=mean,maximum=maximum))
 expectedco=[start]+[start+k*.005 for k in range(steps)];assert len(co)==len(expectedco) and all(abs(r['native_pre_advance_state_time_s']-t)<1e-8 for r,t in zip(co,expectedco))
 inputfiles +=[p for d in ['system','constant','0/polyMesh'] for p in (case/d).rglob('*') if p.is_file()];inputfiles +=[case/'solver.log',case/'force-history.csv',case/'mesh-check.log']
 sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputfiles}
 result=dict(start_time_s=start,end_time_s=start+steps*.005,steps=steps,saved_states=len(times),residual_definition='Last printed initial/final p/Ux/Uy linear residual per physical step after fixed8outer corrections; same as original cylinder record. Not equivalent to Fluent scaled norms.',residual_miss_steps=len(misses),residual_misses=misses,step_records=rows,fields_and_mass=fields,maximum_relative_mass_imbalance=massmax,courant_records=co,courant_time_convention='Native pre-advance state time; initial duplicate report retained. No imposed CFL target for this implicit cylinder model.',numerical_window_pass=not misses and massmax<=1e-4,inherited_startup_not_repaired=True,physical_validation=False,benchmark_pass=False,native_source_sha256=sources)
 (case/'window-audit.json').write_text(json.dumps(result,indent=2)+'\n');return result