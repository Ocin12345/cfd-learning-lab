from pathlib import Path
import csv,hashlib,json
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent
F=R/'exports/figures/2026-11-08-airfoil-TE-paired-v002/force-comparison.csv'
with F.open(newline='') as f:rows=[{k:(v if k=='solver' else float(v)) for k,v in r.items()} for r in csv.DictReader(f)]
mapped={(r['solver'],int(r['level'])):r for r in rows}
gaps=[]
for level in (0,1,2):
 f=mapped['ANSYS Fluent',level];o=mapped['OpenFOAM',level]
 gap=100*(f['Cd']-o['Cd'])/abs(o['Cd'])
 gaps.append(dict(level=level,Cd_difference_percent_OpenFOAM_denominator=gap,cross_code_Cd_learning_target_percent=.5,cross_code_learning_target_pass=abs(gap)<=.5,individual_endpoint_numerical_checks_passed=True))
responses=[]
for code in ('ANSYS Fluent','OpenFOAM'):
 base=mapped[code,0]
 for level in (1,2):
  r=mapped[code,level];responses.append(dict(solver=code,level=level,Cd_change_from_baseline_percent=100*(r['Cd']-base['Cd'])/abs(base['Cd']),Cl_change_from_baseline_percent=100*(r['Cl']-base['Cl'])/abs(base['Cl'])))
p=S/'airfoil-response-metrics-v001.json';assert not p.exists()
record=dict(cross_code=gaps,local_responses=responses,mesh_independence_established=False,physical_validation=False,numerical_cause_identified=False,source_csv=str(F),source_sha256=hashlib.sha256(F.read_bytes()).hexdigest())
p.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
