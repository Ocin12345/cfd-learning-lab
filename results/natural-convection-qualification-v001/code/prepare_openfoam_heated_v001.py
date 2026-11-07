from pathlib import Path
import shutil,json,re,hashlib
from cavity_buoyancy_common_v002 import *
base=ROOT/'native/runs/natural-convection-qualification-openfoam-v002';src=base/'control';dst=base/'heated'
assert all(json.loads((src/'control-checks.json').read_text())['checks'].values());assert not dst.exists();dst.mkdir()
for f in ['system','constant']:shutil.copytree(src/f,dst/f)
poly=src/'constant/polyMesh'
if not poly.exists():shutil.copytree(src/'0/polyMesh',dst/'constant/polyMesh')
(dst/'0').mkdir()
for name in ['U','p','p_rgh','T']:
 s=(src/'300'/name).read_text();s=re.sub(r'location\s+"300";', 'location "0";',s)
 if name=='T':
  s=s[:s.index('boundaryField')]+'''boundaryField {
left { type fixedValue; value uniform 300.4; }
right { type fixedValue; value uniform 299.6; }
top { type zeroGradient; }
bottom { type zeroGradient; }
frontAndBackPlanes { type empty; }
}
'''
 (dst/'0'/name).write_text(s)
p=dst/'system/controlDict';s=p.read_text().replace('endTime 300','endTime 1').replace('deltaT 1;','deltaT 0.02;').replace('writeInterval 100;','writeInterval 1;');p.write_text(s)
p=dst/'system/fvSchemes';s=p.read_text().replace('steadyState','Euler').replace('bounded Gauss','Gauss');p.write_text(s)
p=dst/'system/fvSolution';s=p.read_text().replace('nOuterCorrectors 1','nOuterCorrectors 4').replace('relTol 0.01','relTol 0');s=s[:s.index('relaxationFactors')];p.write_text(s)
write(dst/'declared-settings.json',specification());write(dst/'source-control-hashes.json',{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (src/'300').iterdir() if p.is_file()})
write(dst/'prepared-input-hashes.json',{p.relative_to(dst).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in dst.rglob('*') if p.is_file()})
