"""Aggregate checked airfoil solutions and measured mesh spacing; no CFD runs."""
from pathlib import Path
import argparse,csv,hashlib,json,math
p=argparse.ArgumentParser()
p.add_argument('--repo',type=Path,required=True)
p.add_argument('--metrics-root',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
inputs=[('coarse','airfoil-sensitivity-coarse-v002','metrics-coarse-v002.json'),
 ('medium','lesson03-alpha0-local-v002','metrics-medium.json'),
 ('fine','airfoil-sensitivity-fine055-v002','metrics-fine055.json')]
records=[];hashes={};mesh_source_hashes=set();reference_methods=None
for label,directory,metric_file in inputs:
 root=a.repo/'results'/directory
 checks=json.loads((root/'checks.json').read_text())
 settings=json.loads((root/'numerical-settings.json').read_text())
 geometry=json.loads((root/'geometry-setup.json').read_text())
 metrics=json.loads((a.metrics_root/metric_file).read_text())
 assert checks['numerical_checks_passed'] and metrics['checks_passed']
 assert checks['Re']==100 and checks['angle_degrees']==0
 assert checks['cells']==metrics['cells']==geometry['cells']
 methods=json.dumps(settings['methods'],sort_keys=True)
 if reference_methods is None:reference_methods=methods
 assert methods==reference_methods,'Numerical methods differ between meshes'
 for item in metrics['inputs'].values():
  assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()==item['sha256'],'Metrics input changed'
 assert geometry['preserved_source_polygon_vertices'] and geometry['local_edge_sizing']
 mesh_source_hashes.add(geometry['source_mesh_sha256'])
 refs=settings['references'];rho=refs['density'];u=refs['velocity'];area=refs['area']
 assert rho==u==area==refs['length']==refs['depth']==1
 with (root/'force-history.csv').open(newline='') as f:force=list(csv.DictReader(f))[-1]
 pressure=float(force['pressure_drag_N'])/(.5*rho*u*u*area)
 viscous=float(force['viscous_drag_N'])/(.5*rho*u*u*area)
 assert math.isclose(pressure+viscous,checks['final_Cd'],abs_tol=1e-7)
 records.append(dict(mesh=label,cells=checks['cells'],nodes=metrics['nodes'],
  size_scale=geometry['size_scale'],Cd=checks['final_Cd'],Cl=checks['final_Cl'],
  pressure_Cd=pressure,viscous_Cd=viscous,minimum_orthogonal_quality=checks['minimum_orthogonal_quality'],
  relative_mass_imbalance=checks['relative_mass_imbalance'],
  final_scaled_residuals=checks['final_scaled_residuals'],wall_faces=metrics['wall_face_count'],
  median_wall_edge_m=metrics['wall_regions']['all']['edge_length_m']['median'],
  median_wall_normal_height_m=metrics['wall_regions']['all']['adjacent_triangle_normal_height_m']['median'],
  median_leading_edge_height_m=metrics['wall_regions']['leading_edge_x_le_0.05']['adjacent_triangle_normal_height_m']['median'],
  median_trailing_edge_height_m=metrics['wall_regions']['trailing_edge_x_ge_0.95']['adjacent_triangle_normal_height_m']['median'],
  median_wake_sqrt_area_m=metrics['wake']['triangle_sqrt_area_m']['median']))
 for file in [root/'checks.json',root/'force-history.csv',root/'geometry-setup.json',root/'numerical-settings.json',a.metrics_root/metric_file]:
  hashes[f'{label}/{file.name}']=hashlib.sha256(file.read_bytes()).hexdigest()
assert len(mesh_source_hashes)==1
assert all(records[i]['cells']<records[i+1]['cells'] for i in range(2))
assert all(records[i]['wall_faces']<records[i+1]['wall_faces'] for i in range(2))
for metric in ['median_wall_edge_m','median_wall_normal_height_m','median_leading_edge_height_m','median_trailing_edge_height_m','median_wake_sqrt_area_m']:
 assert all(records[i][metric]>records[i+1][metric] for i in range(2)),metric
fine=records[-1]
for r in records:
 r['Cd_difference_from_fine_percent']=100*(r['Cd']-fine['Cd'])/abs(fine['Cd'])
change=abs(records[1]['Cd_difference_from_fine_percent'])
result=dict(conditions='Re100, zero incidence, unchanged original polygon, outer radius20c',
 records=records,medium_fine_Cd_difference_percent=change,
 declared_Cd_change_target_percent=.5,sensitivity_target_pass=change<.5,
 uniform_refinement=False,GCI_established=False,mesh_independence_established=False,
 domain_independence_established=False,physical_validation=False,
 source_polygon_mesh_sha256=next(iter(mesh_source_hashes)),source_hashes=hashes,
 note='Measured differences for this nonuniform mesh family; inspect pressure and viscous cancellation, lift bias and regional spacing. No formal error bound.')
if a.out.exists():raise FileExistsError('Choose a fresh output directory')
a.out.mkdir(parents=True)
(a.out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
with (a.out/'comparison.csv').open('w',newline='',encoding='utf-8') as f:
 writer=csv.DictWriter(f,fieldnames=[k for k in records[0] if k!='final_scaled_residuals'])
 writer.writeheader();writer.writerows({k:v for k,v in r.items() if k!='final_scaled_residuals'} for r in records)
print(json.dumps(result,indent=2))
