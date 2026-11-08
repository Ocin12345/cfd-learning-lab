from pathlib import Path
import json,hashlib,os,datetime
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).resolve().parent;N=R/'native/runs/showcase-heated-jet-energy-v001'
FBASE=R/'native/runs/heated-impinging-jet-fluent-v002';OBASE=R/'native/runs/heated-impinging-jet-openfoam-v002';PARENT=R/'native/runs/impinging-jet-fluent-v001'
PY=R/'environments/fluent/Scripts/python.exe'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')
def mark(p,stage,**v):
 p=Path(p);q=p.with_suffix('.tmp');dump(q,dict(stage=stage,observed_at=datetime.datetime.now().astimezone().isoformat(),pid=os.getpid(),**v));q.replace(p)
