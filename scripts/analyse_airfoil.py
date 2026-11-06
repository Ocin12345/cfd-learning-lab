"""Check a completed Fluent airfoil run and export small, reviewable records."""
from pathlib import Path
import argparse,csv,hashlib,json,math,re,shutil
p=argparse.ArgumentParser()
p.add_argument('--run',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();run=a.run
targets=json.loads((run/'declared-targets.json').read_text())
summary=json.loads((run/'summary.json').read_text())
geometry=json.loads((run/'geometry-setup.json').read_text())
settings=json.loads((run/'loaded-settings.json').read_text())
assert settings['viscous']['model']=='laminar'
assert settings['material']['density']['value']==1 and settings['material']['viscosity']['value']==.01
assert settings['material']['density']['option']==settings['material']['viscosity']['option']=='value'
assert settings['general']['solver']['time']=='steady'
assert settings['general']['solver']['two_dim_space']=='planar'
assert not settings['general']['operating_conditions']['gravity']['enable']
assert settings['inlet']['momentum']['velocity_specification_method']=='Components'
assert [v['value'] for v in settings['inlet']['momentum']['velocity_components']]==[1,0]
assert settings['outlet']['momentum']['gauge_pressure']['value']==0
for key in ['density','velocity','length','area','depth']:
    assert settings['references'][key]==1,('Coefficient reference mismatch',key)
actual_re=settings['material']['density']['value']*settings['inlet']['momentum']['velocity_components'][0]['value']*settings['references']['length']/settings['material']['viscosity']['value']
assert actual_re==100
assert geometry['reynolds_number']==100 and geometry['angle_degrees']==0
with (run/'force-history.csv').open(newline='') as f:
    rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
assert [r['additional_iteration'] for r in rows]==list(range(targets['additional_iterations']+1))
assert all(math.isfinite(v) for row in rows for v in row.values())
for row in rows:
    for direction,key in [('drag','drag_coefficient'),('lift','lift_coefficient')]:
        native=(run/f"{direction}-{int(row['additional_iteration']):03d}.txt").read_text()
        net=re.findall(r'^\s*Net\s+((?:[-+0-9.eE]+\s+){5}[-+0-9.eE]+)\s*$',native,re.M)
        assert net and float(net[-1].split()[-1])==row[key]
assert rows[-1]==summary['final']
last=rows[-targets['last_window_samples']:]
cd=[r['drag_coefficient'] for r in last];cl=[r['lift_coefficient'] for r in last]
cd_range=(max(cd)-min(cd))/abs(math.fsum(cd)/len(cd))
assert cd_range<=targets['drag_relative_range_max']
assert max(cl)-min(cl)<=targets['lift_absolute_range_max']
assert max(map(abs,cl))<=targets['zero_angle_abs_lift_diagnostic_target']
text=(run/'solver-transcript.txt').read_text(errors='replace')
quality=float(re.search(r'Minimum Orthogonal Quality\s*=\s*([-+0-9.eE]+)',text)[1])
assert quality>=targets['minimum_orthogonal_quality']
residual_rows=re.findall(r'^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)',text,re.M)
final=residual_rows[-1];residuals=list(map(float,final[1:]))
assert max(residuals)<=targets['residual_target'] and 'solution is converged' in text
mass=(run/'mass-balance.txt').read_text()
net=float(re.search(r'^\s*Net\s+([-+0-9.eE]+)',mass,re.M)[1])
inlet=float(re.search(r'^\s*inlet\s+([-+0-9.eE]+)',mass,re.M)[1])
assert abs(net/inlet)<=targets['relative_mass_imbalance_max']
for suffix in ['cas','dat']:assert (run/f'airfoil-Re100-final.{suffix}.h5').is_file()
checks=dict(cells=geometry['cells'],Re=100,angle_degrees=0,
            final_iteration=int(final[0]),final_scaled_residuals=residuals,
            minimum_orthogonal_quality=quality,relative_mass_imbalance=abs(net/inlet),
            final_Cd=rows[-1]['drag_coefficient'],final_Cl=rows[-1]['lift_coefficient'],
            last20_Cd_relative_range=cd_range,last20_Cl_absolute_range=max(cl)-min(cl),
            numerical_checks_passed=True,mesh_independence_established=False,physical_validation=False,
            source_snapshot_sha256=hashlib.sha256((run/'executed-script.py').read_bytes()).hexdigest())
if a.out.exists():raise FileExistsError('Use a fresh output directory')
a.out.mkdir(parents=True)
for name in ['summary.json','declared-targets.json','force-history.csv','residuals.xy','mass-balance.txt','geometry-setup.json']:
    shutil.copy2(run/name,a.out/name)
(a.out/'numerical-settings.json').write_text(json.dumps(settings,indent=2)+'\n',encoding='utf-8')
(a.out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8')
print(json.dumps(checks,indent=2))
