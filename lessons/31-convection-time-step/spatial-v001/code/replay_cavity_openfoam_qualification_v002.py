from pathlib import Path
import hashlib,json,re
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-convection-spatial-v001';O=N/'openfoam-qualification-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 original=json.loads((S/'openfoam-heated-qualification-inputs-v001.json').read_text());target=O/'constant/polyMesh/cellProc';source=N/'openfoam-control-v004/0/polyMesh/cellProc'
 assert sha(source)==original[str(target)]
 before=source.read_bytes();expected,count=re.subn(rb'(location\s+")0/polyMesh(";)',rb'\g<1>constant/polyMesh\2',before);assert count==1 and expected==target.read_bytes()
 for p,h in original.items():
  if Path(p)!=target:assert sha(Path(p))==h,p
 # Exact metadata correction only; all decomposition labels and all other973 source files unchanged.
 corrected=dict(original);corrected[str(target)]=sha(target);p=S/'openfoam-heated-postdecomposition-inputs-v002.json';assert not p.exists();p.write_text(json.dumps(corrected,indent=2)+'\n')
 (S/'decomposition-header-correction-v002.json').write_text(json.dumps(dict(original_manifest='openfoam-heated-qualification-inputs-v001.json',original_failure_retained=True,allowed_exact_rewrite='cellProc location header0/polyMesh -> constant/polyMesh; entire remaining bytes identical',before_sha256=sha(source),after_sha256=sha(target),all_other_inputs_unchanged=len(original)-1,no_solver_rerun=True,no_numerical_limit_relaxed=True),indent=2)+'\n')
 from audit_cavity_openfoam_qualification_v001 import audit
 a=audit(O);assert a['numerical_pass'];print('Complete native100step/101state auditPASS',a['maximum_native_courant'],a['Nu_hot_at_1s'])
if __name__=='__main__':main()