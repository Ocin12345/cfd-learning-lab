from pathlib import Path
import datetime,hashlib,json,os,re,shutil
S=Path(__file__).resolve().parent

def native_path(name):
 if os.name=='nt':return Path(name)
 s=name.replace('\\','/')
 if re.match(r'^[Ee]:/',s):return Path('/mnt/e/'+s[3:])
 return Path(s)

def verify():
 a=json.loads((S/'openfoam-qualification-inputs-v001.json').read_text())
 for name,h in a.items():
  p=native_path(name)
  if hashlib.sha256(p.read_bytes()).hexdigest()!=h:raise ValueError('OpenFOAM frozen input changed: '+str(p))
 return len(a)
def main():
 count=verify();root=native_path('E:/CFD-Learning-Lab');case=root/'native/runs/showcase-cylinder-convergence-v001/openfoam-qualification-v001'
 assert not (case/'solver.log').exists() and not (case/'run-status.txt').exists()
 assert re.search(r'numberOfSubdomains\s+2\s*;', (case/'system/decomposeParDict').read_text())
 assert shutil.disk_usage(root).free>=10*1024**3
 if os.name!='nt':
  m=re.search(r'MemAvailable:\s+(\d+)\s+kB',Path('/proc/meminfo').read_text());assert m and int(m[1])*1024>=4*1024**3
 print(json.dumps(dict(checked_at=datetime.datetime.now().astimezone().isoformat(),verified_inputs=count,E_free_bytes=shutil.disk_usage(root).free,passed=True)))
if __name__=='__main__':main()