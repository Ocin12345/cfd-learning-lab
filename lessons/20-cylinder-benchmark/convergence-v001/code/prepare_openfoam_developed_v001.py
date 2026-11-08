from pathlib import Path
import hashlib,json,re,shutil
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-cylinder-convergence-v001';old=N/'openfoam-qualification-v001';new=N/'openfoam-developed-v001';assert not new.exists();new.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a=json.loads((old/'window-audit.json').read_text());assert a['numerical_window_pass'] and a['steps']==100 and a['saved_states']==101
for p,h in a['native_source_sha256'].items():assert sha(Path(p))==h
F=R/'exports/figures/2026-11-08-cylinder-paired-qualification-v001';(S/'paired-qualification-visual-review-v001.json').write_text(json.dumps(dict(reviewed_by='Codex',actual_figure=str(F/'paired-cylinder-check.png'),sha256=sha(F/'paired-cylinder-check.png'),visually_inspected=True,shown_in_chat=True,shared_speed_scale_and_physical_time=True,native_force_phase_difference_retained=True,fluent_cycle_drift_readable=True,physical_validation=False),indent=2)+'\n')
proof={}
for folder in ['constant','system','0','4.5','processor0/0','processor0/4.5','processor1/0','processor1/4.5']:
 src=old/folder;assert src.is_dir(),src;dst=new/folder;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(src,dst)
 for p in src.rglob('*'):
  if p.is_file():q=new/p.relative_to(old);assert sha(p)==sha(q);proof[str(q)]=dict(source=str(p),sha256=sha(p))
c=new/'system/controlDict';before=c.read_text();after=before
for key,value in [('startTime','4.5'),('endTime','8')]:after,count=re.subn(r'(?m)^(\s*'+key+r'\s+)[^;]+;',lambda m:m[1]+value+';',after,count=1);assert count==1
assert before[before.index('functions'):]==after[after.index('functions'):];c.write_text(after)
(S/'openfoam-developed-copy-proof-v001.json').write_text(json.dumps(dict(initial_copy_files_verified=proof,only_modified_input='system/controlDict startTime4.5/endTime8',source_control_sha256=sha(old/'system/controlDict'),configured_control_sha256=sha(c),functions_unchanged=True,physical_numerics_unchanged=True),indent=2)+'\n')
s=(S/'verify_openfoam_qualification_v001.py').read_text().replace('openfoam-qualification-inputs-v001','openfoam-developed-inputs-v001').replace('openfoam-qualification-v001','openfoam-developed-v001');(S/'verify_openfoam_developed_v001.py').write_text(s)
s=(S/'run_cylinder_openfoam_qualification_v001.sh').read_text().replace('openfoam-qualification-v001','openfoam-developed-v001').replace('cfd-cylinder-qualification-v001','cfd-cylinder-developed-v001').replace('verify_openfoam_qualification_v001','verify_openfoam_developed_v001').replace("4.000001:4.5","4.500001:8").replace('qualification','development interval')
with (S/'run_cylinder_openfoam_developed_v001.sh').open('w',newline='\n') as f:f.write(s)
s=(S/'wait_and_run_openfoam_qualification_v001.py').read_text().replace('verify_openfoam_qualification_v001','verify_openfoam_developed_v001').replace('scheduled-openfoam-qualification-v001-status','cylinder-developed-openfoam-v001-status').replace('openfoam-qualification-v001','openfoam-developed-v001').replace('openfoam-qualification-launch-admission','openfoam-developed-launch-admission').replace('openfoam-qualification-wrapper','openfoam-developed-wrapper').replace('run_cylinder_openfoam_qualification_v001.sh','run_cylinder_openfoam_developed_v001.sh').replace("audit_openfoam(N/'openfoam-developed-v001',4.,100)","audit_openfoam(N/'openfoam-developed-v001',4.5,700)").replace("'openfoam_qualification_pass'","'openfoam_development_interval_pass'").replace("'openfoam_qualification_numerical_miss'","'openfoam_development_interval_numerical_miss'").replace('steps=100,states=101','steps=700,states=701').replace("'running_openfoam_qualification'","'running_openfoam_development_interval'")
# Require the selected accepted half-second checkpoint as well as the existing full Fluent gate.
s=s.replace("verify();mark('waiting_for_accepted_fluent_interval')","verify();qa=json.loads((N/'openfoam-qualification-v001/window-audit.json').read_text());assert qa['numerical_window_pass'] and qa['steps']==100;mark('waiting_for_accepted_fluent_interval')")
compile(s,'developedOcontroller','exec');(S/'run_cylinder_openfoam_developed_v001.py').write_text(s)
files=json.loads((S/'openfoam-qualification-inputs-v001.json').read_text())
for n in ['verify_openfoam_developed_v001.py','run_cylinder_openfoam_developed_v001.sh','run_cylinder_openfoam_developed_v001.py','openfoam-developed-copy-proof-v001.json','paired-qualification-visual-review-v001.json']:
 p=S/n;files[str(p)]=sha(p)
for p,h in a['native_source_sha256'].items():files[p]=h
for p in new.rglob('*'):
 if p.is_file():files[str(p)]=sha(p)
for n in ['window-audit.json','force-history.csv']:
 p=old/n;files[str(p)]=sha(p)
for p in (S/'analysis-fluent-v001').iterdir():
 if p.is_file():files[str(p)]=sha(p)
(S/'openfoam-developed-inputs-v001.json').write_text(json.dumps(files,indent=2)+'\n')
assert b'\r' not in (S/'run_cylinder_openfoam_developed_v001.sh').read_bytes()
with (S/'progress.md').open('a') as f:f.write('\nFullFluent4..8interval: all800newstepsPASS/801uniquefields; duplicateboundaryU/V/Pidentical/forcevaluesmatch.11completecycles; last6criterionPASS(periodrange0.002265%,liftamplituderange0.033566%,meanCdrange0.002914%). Inherited4sinitialresidualnotrepaired. OpenFOAM100step4..4.5qualificationPASS101states, actualpaired4.5sfield/nativeforce/Fcyclefigureinspected/shown. ExactO4.5root+bothnativepartitioncheckpoints/history/meshes copiedandverified; onlystart/endcontrols changed for700steps to8. No thresholdrelaxed; finalfullpairrecord/reviewpending.\n')
print(json.dumps(dict(new_case=str(new),frozen_inputs=len(files),qualification_mass_max=a['maximum_relative_mass_imbalance'])))