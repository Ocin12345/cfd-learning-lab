from pathlib import Path
import hashlib,json,re,shutil
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');R=Path('E:/CFD-Learning-Lab');old=R/'native/runs/2026-11-01-vortex-shedding-v001/openfoam-extended';new=R/'native/runs/showcase-cylinder-convergence-v001/openfoam-qualification-v001';assert not new.exists();new.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
proof={}
for folder in ['constant','system','0','4','processor0/0','processor0/4','processor1/0','processor1/4']:
 src=old/folder;assert src.is_dir(),src;dst=new/folder;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(src,dst)
 for p in src.rglob('*'):
  if p.is_file():q=new/p.relative_to(old);assert sha(p)==sha(q);proof[str(q)]=dict(source=str(p),sha256=sha(p))
c=new/'system/controlDict';before=c.read_text();after=before
for key,value in [('startFrom','startTime'),('startTime','4'),('endTime','4.5'),('writeInterval','1')]:
 # Change first/top-level occurrence only, preserving functions block.
 after,count=re.subn(r'(?m)^(\s*'+key+r'\s+)[^;]+;',lambda m:m[1]+value+';',after,count=1);assert count==1
assert before[before.index('functions'):]==after[after.index('functions'):]
c.write_text(after,encoding='utf-8')
(S/'openfoam-qualification-copy-proof-v001.json').write_text(json.dumps(dict(initial_copy_files_verified=proof,only_modified_input='system/controlDict',modified_source_sha256=sha(old/'system/controlDict'),modified_copy_sha256=sha(c),functions_block_unchanged=True,originals_preserved=True),indent=2)+'\n')
a=(S/'audit_cylinder_openfoam_v001.py').read_text();a=a.replace("case=Path(case);text=", "case=Path(case);assert 'Mesh OK.' in (case/'mesh-check.log').read_text();text=")
(S/'audit_cylinder_openfoam_v002.py').write_text(a,encoding='utf-8')
with (S/'progress.md').open('a') as f:f.write('\nOpenFOAMpreparationonly: exactrootand2partition4srestart/historyand0meshes copied/byteverified, originals preserved. Onlytop-levelcontrolstart/end/write modified, functionblockidentical. Prepared100step4..4.5qualification, NOTLAUNCHEDbeforeFluent8sacceptance. Nativeallstepcoverage/residualmissingregressionsREDthenGREEN2/2; sameoriginalterminalinitial/finalresidualdefinition1e-5, enforcedallnewsteps. Courantnativepre-advancetimesandinitialduplicatereportexplicit.\n')
print('Prepared OpenFOAM qualification case: '+str(new))