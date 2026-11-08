"""Read-only native Fluent audit for the 5-degree locally refined airfoil."""
from pathlib import Path
import argparse,hashlib,json,math
import analyse_openfoam_airfoil_v002 as av

def check_targets(t):
 required=dict(residual_target=1e-7,relative_mass_imbalance_max=1e-6,drag_relative_range_max=1e-4,lift_absolute_range_max=1e-5,minimum_orthogonal_quality=.2,last_window_samples=20,additional_iterations=100)
 for k,v in required.items():
  if t.get(k)!=v:raise ValueError('Changed acceptance target: '+k)

def applicable_checks(raw,angle):
 if angle!=5:raise ValueError('This audit admits only declared 5 degrees')
 return {k:v for k,v in raw.items() if k!='fluent.zero_angle_lift_diagnostic'}

def audit(run,mesh,output):
 R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent
 if output.exists():raise FileExistsError(output)
 assert json.loads((run/'status.json').read_text())['stage']=='completed'
 criteria=json.loads((run/'declared-targets.json').read_text());check_targets(criteria)
 admission=json.loads((S/'native-quality-and-visual-admission-v001.json').read_text())
 for name,h in admission['source_sha256'].items():assert av.sha256(Path(name))==h, name
 import airfoil_prelaunch_v003 as gate
 gate.verify_inputs()
 base=R/'native/runs/2026-11-06-airfoil-incidence-paired-v001'
 template=base/'fluent-alpha5/loaded-settings.json'
 targets_path=base/'openfoam-alpha5-v002/declared-comparison.json'
 targets=json.loads(targets_path.read_text());assert targets['angle_degrees']==5 and targets['reynolds_number']==100
 xy,edges,cells,_=av.read_source_mesh(mesh)
 assert cells in (119580,159300)
 assert json.loads((run/'geometry-setup.json').read_text())==json.loads((mesh.parent/'setup.json').read_text())
 settings=json.loads((run/'loaded-settings.json').read_text());assert settings['methods']==json.loads(template.read_text())['methods']
 # Freeze all native evidence before independent replay; no mutation of any run.
 paths=[x for x in run.iterdir() if x.is_file() and x.suffix in ('.json','.txt','.csv','.xy','.h5','.py')]
 paths += [mesh,mesh.parent/'setup.json',mesh.parent/'geometry-proof.json',template,targets_path,Path(__file__),Path(av.__file__),Path(av.__file__).with_name('plot_openfoam_airfoil_fields_v002.py'),S/'native-quality-and-visual-admission-v001.json']
 hashes={str(p):av.sha256(p) for p in paths}
 record,history,residuals=av.check_fluent(run,targets,xy,edges,cells,retain_failed_checks=True)
 raw=record['check_results'];checks=applicable_checks(raw,5)
 # Independently replay the unit-reference normalization from reported forces.
 for row in history:
  for d in ('drag','lift'):assert math.isclose(row[d+'_coefficient'],2*row[d+'_force_N'],rel_tol=0,abs_tol=1e-7)
 assert all(av.sha256(Path(p))==h for p,h in hashes.items()),'Evidence changed during read-only audit'
 record.update(cells=cells,angle_degrees=5,Re=100,check_results=checks,failed_checks=[k for k,v in checks.items() if not v],numerical_checks_passed=all(checks.values()),numerical_convergence_passed=all(checks.values()),symmetry_diagnostic_passed=None,zero_incidence_symmetry_applicable=False,unapplied_zero_angle_diagnostic=raw.get('fluent.zero_angle_lift_diagnostic'),unit_force_normalization_replayed=True,source_run=str(run),source_mesh=str(mesh),source_sha256=hashes,physical_validation=False,mesh_independence_established=False,numerical_cause_identified=False)
 output.mkdir(parents=True);(output/'checks.json').write_text(json.dumps(record,indent=2)+'\n')
 print(json.dumps({k:record[k] for k in ['cells','Cd','Cl','final_iteration','final_scaled_residuals','relative_mass_imbalance','Cd_relative_range','Cl_absolute_range','minimum_orthogonal_quality','failed_checks','numerical_checks_passed']},indent=2))
 return record

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--mesh',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 r=audit(a.run,a.mesh,a.out)
 if not r['numerical_checks_passed']:raise SystemExit(2)
