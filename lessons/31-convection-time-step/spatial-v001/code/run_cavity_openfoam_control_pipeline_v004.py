from pathlib import Path
from datetime import datetime
import hashlib,json,os,subprocess,shutil,time,traceback,re
import numpy as np
from cavity_spatial_common_v001 import *
from cavity_prelaunch_v001 import prelaunch,verify
from wing_fluent_v001 import available_memory
from cavity_native_fields_v001 import foam,metrics
from analyse_nozzle_openfoam_v002 import field,FIELDS,patch_values
from plot_openfoam_airfoil_fields_v002 import coordinate_bijection
S=Path(__file__).parent;N=ROOT/'native/runs/showcase-convection-spatial-v001';O=N/'openfoam-control-v004';F=N/'fluent-qualification-v004'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mark(stage,**v):
 p=S/'openfoam-control-v004-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def parent():
 a=json.loads((F/'qualification-audit-v002.json').read_text());assert a['numerical_pass']
 assert json.loads((F/'control-checks.json').read_text())['control']['max_speed_m_s']<=1e-6
 for n,h in a['source_sha256'].items():assert sha(Path(n))==h,n
 return dict(path=str(F/'qualification-audit-v002.json'),sha256=sha(F/'qualification-audit-v002.json'))
def gate():
 while available_memory()<4*2**30 or shutil.disk_usage(ROOT).free<10*2**30:mark('waiting_resources');time.sleep(30)
 p=parent();a=prelaunch('openfoam-control-inputs-v004.json');return p,a
def run(stage):
 with (S/('openfoam-control-'+stage+'-v004.log')).open('xb') as f:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/convection-spatial-v001/run_cavity_openfoam_control_v004.sh',stage],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,(stage,p.returncode)
def main():
 assert not (S/'openfoam-control-v004-status.json').exists();p,a=gate();mark('reusing_verified_native_mesh_before_control',resources=a,parent=p)
 with np.load(ROOT/'native/meshes/showcase-convection-spatial-v001/fine-v001/geometry.npz') as g:
  from plot_nozzle_native_v004 import foam_list
  from cavity_extrusion_geometry_v001 import verify_extrusion
  C=field(O/'0/C',9216,*FIELDS['C']);mesh=O/'0/polyMesh';mesh=mesh if mesh.exists() else O/'constant/polyMesh'
  _,body=foam_list(mesh/'points');points=np.fromstring(body.replace('(',' ').replace(')',' '),sep=' ').reshape(-1,3);mapping=verify_extrusion(C,points,g['centres'])
 assert 'Mesh OK.' in (O/'mesh-check.log').read_text() and re.search(r'\bcells:\s+9216\b',(O/'mesh-check.log').read_text())
 meshfiles=[p for d in O.rglob('polyMesh') for p in d.rglob('*') if p.is_file()];native={str(p):sha(p) for p in meshfiles};native[str(O/'0/C')]=sha(O/'0/C')
 (O/'native-mesh-audit.json').write_text(json.dumps(dict(bijection=mapping,extrusion_depth_m=.005,source_sha256=native),indent=2)+'\n')
 p,a=gate()
 for n,h in native.items():assert sha(Path(n))==h,n
 assert available_memory()>=4*2**30 and shutil.disk_usage(ROOT).free>=10*2**30
 mark('running_600_steady_isothermal_iterations',resources=a,parent=p);run('solve');verify('openfoam-control-inputs-v004.json');parent()
 for n,h in native.items():assert sha(Path(n))==h,n
 q=foam(O,1200,96);m=metrics(q,96);drift=abs(m['mass_per_depth_kg_m']/(RHO0*L*L)-1)
 checks=dict(quiet=m['Umax_m_s']<=1e-6,temperature=max(abs(m['Tmin_K']-300),abs(m['Tmax_K']-300))<=1e-4,mass=drift<=1e-5,eos=m['eos_relative_error']<=1e-6)
 walls={name:float(abs(patch_values(O/'1200/phi',name,96)).sum())/DEPTH for name in WALLS};checks['impermeable_walls']=sum(walls.values())<=1e-12
 log=(O/'solver.log').read_text();steps=[int(float(t)) for t in re.findall(r'^Time = ([0-9.eE+-]+)s?\s*$',log,re.M)];assert steps==list(range(601,1201))
 record=dict(control_pass=all(checks.values()),steady_control_prefix_iterations=600,additional_control_iterations=600,total_control_iterations=1200,not_physical_heating_time=True,metrics=m,mass_relative_drift=drift,check_results=checks,native_wall_mass_absolute_per_depth_kg_ms=walls,source_sha256={str(p):sha(p) for p in O.rglob('*') if p.is_file()})
 (O/'control-audit-v001.json').write_text(json.dumps(record,indent=2)+'\n');mark('isothermal_control_pass' if record['control_pass'] else 'isothermal_control_numerical_miss',audit=str(O/'control-audit-v001.json'),heated_qualification_not_launched=True,full20s_not_admitted=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise