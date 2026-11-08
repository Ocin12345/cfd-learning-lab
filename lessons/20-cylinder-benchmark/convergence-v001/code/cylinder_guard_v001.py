from pathlib import Path
import hashlib,json

def verify_inputs(manifest,out=None):
 if out is not None:
  out=Path(out).resolve()
  if out.drive.upper()!='E:' or out.exists():raise ValueError('Fresh E: output required')
 sources=json.loads(Path(manifest).read_text())
 for name,expected in sources.items():
  if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=expected:raise ValueError('Frozen source changed: '+name)
 return len(sources)