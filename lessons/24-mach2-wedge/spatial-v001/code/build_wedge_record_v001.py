from pathlib import Path
from datetime import datetime
import hashlib,json,csv,shutil
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/mach2-wedge-spatial-v001');R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-mach2-wedge-spatial-v001'
O=R/'repo/cfd-learning-lab/lessons/24-mach2-wedge/spatial-v001'
assert not O.exists();O.mkdir();copied={}
def copy(p,relative):
 q=O/relative;q.parent.mkdir(parents=True,exist_ok=True);assert not q.exists();shutil.copyfile(p,q);assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256(q.read_bytes()).digest();copied[relative]=dict(source=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
for p in S.glob('*'):
 if p.suffix in ('.py','.sh') or p.name in ['SPEC.md','PLAN.md','LATE-CFL-ADMISSION-v001.md','primary-sources-v001.json','baseline-native-hashes-v001.json'] or p.name.endswith('-inputs-v001.json') or p.name.endswith('-red-v001.txt') or p.name.endswith('-green-v001.txt') or p.name in ['foam-mapping-green-v002.txt','foam-mapping-green-v001.txt']:
  copy(p,('code/' if p.suffix in ('.py','.sh') else 'methods/')+p.name)
for directory in ['2026-11-08-wedge-refinement-mesh-v001','2026-11-08-wedge-spatial-paired-v001','2026-11-08-wedge-spatial-numerics-v001','2026-11-08-wedge-late-cfl-v002']:
 for p in (R/'exports/figures'/directory).glob('*'):
  if p.is_file():copy(p,('figures/' if p.suffix in ('.png','.pdf') else 'data/')+p.name)
selected=json.loads((S/'wedge-selected-measurements-v001.json').read_text());late=json.loads((N/'openfoam-late-cfl-v001/endpoint-audit.json').read_text())
copy(S/'wedge-selected-measurements-v001.json','data/selected-measurements.json')
for name in ['fluent-fine-v001','openfoam-fine-v001','openfoam-late-cfl-v001']:copy(N/name/'endpoint-audit.json','checks/'+name+'-audit.json')
for name in ['fluent-quality-v001/native-mesh-quality.json','openfoam-fine-v001/native-mesh-audit.json','openfoam-late-cfl-v001/checkpoint-provenance.json']:
 copy(N/name,'checks/'+Path(name).name)
for sub in ['0/p','0/T','0/U','constant/physicalProperties','constant/momentumTransport','system/blockMeshDict','system/controlDict','system/fvSchemes','system/fvSolution','system/decomposeParDict']:
 if (N/'openfoam-fine-v001'/sub).is_file():copy(N/'openfoam-fine-v001'/sub,'inputs/openfoam-fine/'+sub)
rows=[]
for r in selected:
 for w in r['width']['columns']:rows.append(dict(code=r['code'],grid=r['resolution'],column=w['column'],x_m=w['x_mean_m'],width_m=w['width_m'],projected_column_interval_m=w['projected_column_centroid_interval_m'],width_in_projected_column_intervals=w['width_in_projected_column_intervals']))
for w in late['width']['columns']:rows.append(dict(code='OpenFOAM',grid='Fine late maxCo0.15',column=w['column'],x_m=w['x_mean_m'],width_m=w['width_m'],projected_column_interval_m=w['projected_column_centroid_interval_m'],width_in_projected_column_intervals=w['width_in_projected_column_intervals']))
with (O/'data/shock-width-columns.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
summary=dict(study=24,new_models=0,conditions=dict(Mach=2,wedge_angle_deg=10,gamma=1.4,R_J_kgK=287.05,p_Pa=101325,T_K=300,domain_m=[1,.82],corner_x_m=.2,OpenFOAM_depth_m=.01,assumptions='2D perfect-gas Euler/slip, no viscosity or turbulence'),cases=[dict(code=r['code'],grid=r['resolution'],cells=r['cells'],metrics=r['metrics'],native=r['native'],width={k:v for k,v in r['width'].items() if k!='columns'}) for r in selected],late_diagnostic=dict(time_steps=late['time_steps'],saved_times_s=late['saved_times_s'],native_metrics=late['measurements'],fine_main_to_late_endpoint_relative_L2=late['fine_main_to_late_endpoint_relative_L2'],median_width_m=late['width']['median_width_m'],width_change_fraction=late['width']['median_width_m']/selected[-1]['width']['median_width_m']-1,local_signed_pressure_difference_percent_range=[-5.459251555884516,7.376442262028118],independent_full_rerun=False,scope='Only7..8ms from exact main fine native checkpoint'),spatial_and_adaptive_time_effects_confounded_for_OpenFOAM=True,physical_validation=False,new_upload=False,owner_execution_or_reflection_claimed=False)
(O/'data/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
text='''# A sharper captured shock on a finer grid

This is a refinement of existing **study24**, with actual ANSYS Fluent and OpenFOAM fields. It adds no numbered model. The mesh doubles from200×112 (22,400 cells) to400×224 (89,600 cells) while retaining the same domain, straight10° wedge and boundary conditions.

![Native density-gradient comparison](figures/wedge-coarse-fine-shock-detail.png)

Both fine endpoints passed the original numerical limits. The median connected10–90% numerical pressure-transition width changed from **12.163 to6.094mm in Fluent** and **13.630 to6.804mm in OpenFOAM**. The smaller thickness in metres describes shock capture by these schemes; it is not a molecular shock thickness or a measured order of convergence.

| Code / grid | Cells | Fitted angle (°) | Downstream/upstream pressure | Median numerical width (mm) |
|---|---:|---:|---:|---:|
| Fluent coarse |22,400|39.316659|1.70672068|12.163|
| Fluent fine |89,600|39.314617|1.70662373|6.094|
| OpenFOAM coarse |22,400|39.319189|1.70698758|13.630|
| OpenFOAM fine |89,600|39.320702|1.70662587|6.804|

The ideal weak-shock values are39.313932° and1.70657860 for this perfect-gas setup. The OpenFOAM angle error increases slightly on this refinement; improvement is not monotonic for every measured quantity. The ideal relations are an analytical consistency check, not experimental validation. [NASA describes the inviscid oblique-shock assumptions and jump relations](https://www.grc.nasa.gov/www/BGH/oblique.html).

[Full native pressure/Mach comparison](figures/wedge-coarse-fine-native-fields.png) · [Native shock profiles and widths](figures/wedge-shock-width-and-profiles.png) · [Numerical evidence](figures/wedge-numerical-checks.png) · [Mesh close-up](figures/wedge-refinement-mesh.png)

## What passed

Fine Fluent reached the accepted saved steady endpoint at1520 iterations, with maximum native scaled residual1.803e-9, relative boundary mass imbalance1.870e-12 and pressure-field relative L2 change2.837e-5 across its final200iterations. Earlier checkpoints miss one or more endpoint limits and remain preserved. Its native minimum orthogonal quality is0.984818.

Fine OpenFOAM completed9188 actual time steps to8ms and33 saved states including the prescribed uniform initial field. Native maximum Courant number is0.300025, terminal relative mass imbalance4.825e-5 and7–8ms pressure-field change0.00192968 against the unchanged0.002 limit. This narrowly passing finite endpoint check does not establish complete steady-state convergence. Diagonal-solve zero residuals are not nonlinear convergence evidence. Positivity, ideal-gas consistency, mapped native geometry, supersonic outlet, analytical-angle/pressure and saved-state checks also pass.

## What the additional time-step check establishes

OpenFOAM uses adaptive steps at a fixed Courant target. Refining its grid changes both spacing and actual time steps, so its main coarse/fine result combines spatial and temporal effects. A separate2284-step branch reuses the exact saved7ms fine-grid checkpoint and advances to8ms at maxCo0.15, saving five states. It passes the same endpoint limits and retains all inherited startup Courant reports with their actual pre-advance time. See the [native branch evidence](figures/wedge-late-cfl-native-evidence.png).

Between the fine main and half-CFL endpoints, pressure relative L2 difference is0.00168749 and Mach difference0.00068205. Median width changes from6.80417 to6.74317mm (about0.90%). Local signed pressure differences reach−5.46% to+7.38% of freestream pressure near the captured shock. The displayed difference scale is±8% with no clipped cells; the preserved first export had±5% and7 clipped cells. Small domain-weighted differences must not be confused with small changes everywhere. This late-time branch is neither an independent full rerun nor proof of time independence.

## Methods and reproduction

Read [the frozen specification](methods/SPEC.md), [plan](methods/PLAN.md), [late-CFL admission](methods/LATE-CFL-ADMISSION-v001.md), [numerical summary](data/summary.json), [every column width](data/shock-width-columns.csv) and [native bracket samples](data/selected-measurements.json). The pressure plateau/angle masks match the original study. Widths use connected native10/50/90% brackets, linear interpolation in signed-normal coordinates, every prescribed x=.36–.84m column and no smoothing/extrapolation. Cell units mean projected vertical-column centroid increments, not an isotropic cell size. Numerical schlieren uses a tested physical-coordinate Jacobian of cell-centred density, not a simulated optical experiment.

[Executed sources and tests](code/) and [native audit records](checks/) are included. The original coarse study remains in the parent folder. Executable controllers contain recorded E: paths and immutable launch guards; **do not rerun them over completed directories**. To reproduce on this machine, use the existing E: environment and D: Fluent261, review the frozen spec, allocate new versioned paths and manifests, then qualify geometry/native quality before the finite sequential flow runs. Hidden two-rank sessions, one solver at a time, source/parent checks and RAM≥4GiB/E:free≥10GiB admission apply. See the [environment guide](../../../docs/REPRODUCING.md). Native checkpoints, meshes, full logs and installed dependencies stay outside Git; this packet is an execution/methods archive, not a self-contained restart bundle.

All original coarse/native files and failed helper versions are preserved. The failed first centre-mapping test used a tolerance unsupported by the shared API; the new version tightens it to1e-10. No solver acceptance limit was relaxed. The original geometry nodes are bit-identical members of the fine source mesh; native mapping/connectivity and centroid proofs bind the field plots.

## Limits and contributions

Two grids do not establish observed order, GCI, global mesh/domain/time independence or a unique discrepancy cause. Euler/slip excludes boundary layers, viscosity and turbulence. Agreement with an ideal weak shock does not validate an aircraft, real wind-tunnel flow or this portfolio physically.

The owner requested paired CFD models and visuals and authorized improvement of existing showcases. Codex selected and documented this bounded refinement, wrote/tested the scripts, ran and checked the solvers and prepared the figures and record. ANSYS Fluent and OpenFOAM computed the native fields. Owner execution and reflection are not claimed. This is a locally prepared record; no new GitHub upload is claimed.
'''
(O/'README.md').write_text(text,encoding='utf-8')
(S/'wedge-record-copy-verification-v001.json').write_text(json.dumps(copied,indent=2)+'\n')
# Preserve each touched original as a separate exact local backup before appending.
base=json.loads((S/'repository-before-wedge-improvement-v001.json').read_text());repo=R/'repo/cfd-learning-lab';backup=S/'record-root-backups-v001';backup.mkdir()
for relative in ['README.md','docs/PROJECT_HISTORY.md','docs/REPRODUCING.md','lessons/24-mach2-wedge/README.md']:
 p=repo/relative;assert hashlib.sha256(p.read_bytes()).hexdigest()==base[relative];q=backup/relative;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
 target='spatial-v001/README.md' if relative.startswith('lessons/') else ('../lessons/24-mach2-wedge/spatial-v001/README.md' if relative.startswith('docs/') else 'lessons/24-mach2-wedge/spatial-v001/README.md')
 with p.open('a',encoding='utf-8') as f:f.write('\n\n## Mach2 wedge refinement —8November2026\n\nThe [local two-grid comparison]('+target+') adds paired89,600-cell endpoints, native numerical-shock-width measurements and a bounded saved-state late-CFL diagnostic. It preserves the original study and its historical limits; no global mesh/time independence, experimental validation or new upload is claimed.\n')
print('Prepared',O,'with',len(copied),'verified copied artifacts; four backed-up prior pages appended')