from pathlib import Path
import hashlib,json,shutil,re
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-mach2-wedge-spatial-v001';P=N/'openfoam-fine-v001';O=N/'openfoam-late-cfl-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare():
 a=json.loads((P/'endpoint-audit.json').read_text());assert a['numerical_endpoint_pass'] and not a['failed_checks']
 for n,h in a['source_sha256'].items():assert sha(Path(n))==h,n
 assert not O.exists();O.mkdir();copied={}
 for folder in ('0','0.007','constant','system','processor0/0','processor0/0.007','processor0/constant','processor1/0','processor1/0.007','processor1/constant'):
  src=P/folder;assert src.is_dir(),src
  for p in src.rglob('*'):
   if not p.is_file():continue
   relative=p.relative_to(P);q=O/relative;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(p)==sha(q);copied[str(relative)]=sha(p)
 control=O/'system/controlDict';old=control.read_text();new,n=re.subn(r'\bstartTime\s+0;', 'startTime 0.007;',old);assert n==1
 new,n=re.subn(r'\bmaxCo\s+0\.3;', 'maxCo 0.15;',new);assert n==1;control.write_text(new)
 for folder in ('0.007','processor0/0.007','processor1/0.007'):
  t=(O/folder/'uniform/time').read_text();assert re.search(r'name\s+"0.007";',t);assert re.search(r'index\s+8044;',t)
 assert not (O/'0.008').exists()
 record=dict(parent=str(P),parent_audit_sha256=sha(P/'endpoint-audit.json'),copied_native_file_count=len(copied),source_and_original_copy_sha256=copied,changed_control_only=dict(before_sha256=hashlib.sha256(old.encode()).hexdigest(),after_sha256=sha(control),edits=['startTime 0 -> 0.007','maxCo 0.3 -> 0.15']),scope='Exact saved-checkpoint late-time continuation, not independent full run',complete_inputs={str(p):sha(p) for p in O.rglob('*') if p.is_file()})
 (O/'checkpoint-provenance.json').write_text(json.dumps(record,indent=2)+'\n');print('Exact-copy files',len(copied))
if __name__=='__main__':prepare()