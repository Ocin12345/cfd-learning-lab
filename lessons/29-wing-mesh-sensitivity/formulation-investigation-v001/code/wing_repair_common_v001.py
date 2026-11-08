from pathlib import Path
import json,hashlib,datetime,os
R=Path('E:/CFD-Learning-Lab');W=Path(__file__).resolve().parent;N=R/'native/runs/showcase-wing-convergence-v001';F=R/'native/runs/finite-wing-refined-v001/fluent';OFO=F.parent/'openfoam';M=R/'native/meshes/finite-wing-refined-v002';PY=R/'environments/fluent/Scripts/python.exe'
EQ=('continuity','x-velocity','y-velocity','z-velocity')
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')
def mark(p,stage,**v):
 p=Path(p);q=p.with_suffix('.tmp');dump(q,dict(stage=stage,observed_at=datetime.datetime.now().astimezone().isoformat(),pid=os.getpid(),**v));q.replace(p)
