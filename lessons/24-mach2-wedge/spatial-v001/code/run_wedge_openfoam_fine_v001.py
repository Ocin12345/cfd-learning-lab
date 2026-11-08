from pathlib import Path
from datetime import datetime
import hashlib,json,os,shutil,subprocess,time,traceback
import numpy as np
from wing_fluent_v001 import available_memory
from wedge_prelaunch_v001 import prelaunch,verify
from wedge_foam_mapping_v002 import native_mapping
from analyse_nozzle_openfoam_v002 import field,FIELDS
from audit_wedge_openfoam_fine_v001 import audit
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-mach2-wedge-spatial-v001';O=N/'openfoam-fine-v001';G=R/'native/meshes/showcase-mach2-wedge-spatial-v001/fine-v001/geometry.npz'
def mark(stage,**values):
 p=S/'openfoam-fine-v001-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**values),indent=2)+'\n');q.replace(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted_fluent():
 p=N/'fluent-fine-v001/endpoint-audit.json';a=json.loads(p.read_text());assert a['numerical_endpoint_pass'] and not a['failed_checks']
 for n,h in a['source_sha256'].items():assert sha(Path(n))==h,n
 return dict(audit=str(p),sha256=sha(p),verified_native_inputs=len(a['source_sha256']))
def run(mode):
 with (S/('openfoam-'+mode+'-v001.log')).open('xb') as f:r=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/mach2-wedge-spatial-v001/run_wedge_openfoam_fine_v001.sh',mode],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 if r.returncode:raise RuntimeError('Native stage failed: '+mode)
def gate():
 while available_memory()<4*2**30 or shutil.disk_usage(R).free<10*2**30:mark('waiting_resources');time.sleep(30)
 verify('openfoam-fine-inputs-v001.json');a=accepted_fluent();r=prelaunch('openfoam-fine-inputs-v001.json');return a,r
def main():
 assert not (S/'openfoam-fine-v001-status.json').exists();a,r=gate();mark('building_native_mesh',resources=r);run('mesh')
 with np.load(G) as g:
  C=field(O/'0/C',89600,*FIELDS['C']);order,metrics=native_mapping(C,g['centres'],.01)
  error=float(np.max(np.linalg.norm(C[order,:2]-g['centres'],axis=1)));assert error<=1e-10
 np.save(O/'native-to-common-cells.npy',order)
 native={str(p):sha(p) for folder in [O/'constant/polyMesh',O/'processor0/constant/polyMesh',O/'processor1/constant/polyMesh'] for p in folder.glob('*') if p.is_file()}
 native[str(O/'0/C')]=sha(O/'0/C')
 (O/'native-mesh-audit.json').write_text(json.dumps(dict(cell_mapping=metrics,maximum_actual_centroid_error_m=error,native_mesh_sha256=native,geometry_sha256=sha(G)),indent=2)+'\n')
 a,r=gate()
 for n,h in native.items():assert sha(Path(n))==h,n
 assert available_memory()>=4*2**30 and shutil.disk_usage(R).free>=10*2**30
 (S/'openfoam-fine-launch-admission-v001.json').write_text(json.dumps(dict(resources=r,accepted_fluent=a,native_mapping=metrics),indent=2)+'\n')
 mark('running_transient_0_to_8ms',resources=r);run('solve');verify('openfoam-fine-inputs-v001.json')
 for n,h in native.items():assert sha(Path(n))==h,n
 mark('auditing_native_fields');a=audit(O)
 mark('openfoam_endpoint_pass' if a['numerical_endpoint_pass'] else 'openfoam_endpoint_numerical_miss',audit=str(O/'endpoint-audit.json'),failed_checks=a['failed_checks'],late_CFL_diagnostic_not_yet_admitted=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise
