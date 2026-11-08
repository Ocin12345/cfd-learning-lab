from pathlib import Path
import hashlib,importlib.metadata,json,re,shutil,sys
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');REPO=R/'repo/cfd-learning-lab';N=R/'native/runs/showcase-cylinder-convergence-v001';F=R/'exports/figures/2026-11-08-cylinder-developed-paired-v001';V=R/'exports/videos/2026-11-08-cylinder-developed-v001';P=R/'exports/figures/2026-11-08-cylinder-numerical-record-v002';stage=S/'record-staging-v001';dest=REPO/'lessons/20-cylinder-benchmark/convergence-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert not stage.exists() and not dest.exists();base=json.loads((S/'repository-before-cylinder-improvement-v001.json').read_text())
for n,h in base.items():assert sha(REPO/n)==h,'Prior file changed: '+n
assert json.loads((S/'final-cylinder-visual-review-v001.json').read_text())['full_MP4_GIF_frames_all_decoded_and_compared']==161
stage.mkdir();copies={}
def cp(src,name):
 d=stage/name;d.parent.mkdir(parents=True,exist_ok=True);assert not d.exists();shutil.copyfile(src,d);assert sha(src)==sha(d);copies[name]=dict(source=str(src),sha256=sha(d))
for folder in ['analysis-fluent-v001','analysis-openfoam-v001','analysis-pressure-v001']:
 for p in (S/folder).iterdir():
  if p.is_file():cp(p,'data/'+p.name)
for n in ['summary.json','plot-receipt.json']:cp(P/n,'data/'+n)
for n in ['numerical-checks.png','numerical-checks.pdf','forces-and-pressure.png','forces-and-pressure.pdf']:cp(P/n,n)
for k in [0,40,80,120,160]:cp(F/f'frame-{k:03d}.png',f'state-{k:03d}.png')
for k in [0,80,160]:cp(F/f'state-{k:03d}.pdf',f'state-{k:03d}.pdf')
cp(F/'native-field-audit.json','data/native-field-audit.json');cp(V/'cylinder-developed.mp4','cylinder-developed.mp4');cp(V/'compact-gif-v001/cylinder-compact.gif','cylinder-compact.gif')
for src,name in [(V/'video-export-record.json','video-export-record.json'),(V/'gif-preview-v001/gif-export-record.json','external-full-gif-export-record.json'),(V/'compact-gif-v001/gif-export-record.json','compact-gif-export-record.json')]:cp(src,'data/'+name)
for p in sorted(S.glob('*.py')):cp(p,'code/'+p.name)
for p in sorted(S.glob('*.sh')):cp(p,'code/'+p.name)
for n in ['SPEC.md','PLAN.md','COUPLED-ADDENDUM-v001.md','DEVELOPED-INTERVAL-v001.md','progress.md']:cp(S/n,'history/'+n.replace('.md','.txt'))
for p in sorted(S.glob('*.json')):
 if any(k in p.name for k in ['inputs','copy-proof','regression','protocol','visual-review','hypothesis-result','preview-export-rejection','image-review']):cp(p,'provenance/'+p.name)
for p in sorted(S.glob('*green-v001.log')):cp(p,'provenance/'+p.name)
for folder,label in [('fluent-pressure-baseline-v001','pressure-only-baseline'),('fluent-pressure-tight-v001','pressure-only-tight'),('fluent-coupled-pilot-v001','fluent-coupled-pilot'),('fluent-coupled-window-v001','fluent-half-second'),('fluent-developed-v001','fluent-suffix'),('openfoam-qualification-v001','openfoam-qualification'),('openfoam-developed-v001','openfoam-suffix')]:cp(N/folder/'window-audit.json','data/'+label+'-audit.json')
for folder in ['system','constant']:
 for p in (N/'openfoam-developed-v001'/folder).iterdir():
  if p.is_file():cp(p,'openfoam-inputs/'+folder+'/'+p.name)
for n in ['U','p']:cp(N/'openfoam-developed-v001/0'/n,'openfoam-inputs/0/'+n)
cp(REPO/'lessons/20-cylinder-benchmark/README.md','history/earlier-benchmark-README.txt');cp(REPO/'lessons/07-vortex-shedding/README.md','history/earlier-wake-README.txt')
env=dict(python=sys.version,packages={n:importlib.metadata.version(n) for n in ['ansys-fluent-core','numpy','h5py','matplotlib','Pillow']},fluent='ANSYS Fluent2026R1 Student, external D: installation',openfoam='Foundation14, existing CFD-Foundation14 WSL distribution',native_outputs_and_environment_excluded=True);(stage/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
summary=json.loads((P/'summary.json').read_text());f=summary['fluent'];o=summary['openfoam']
readme='''# Cylinder wake: convergence and developed-cycle evidence

![Compact animation from actual native states](cylinder-compact.gif)

[Full-size paired MP4](cylinder-developed.mp4) · [Numerical checks](numerical-checks.png) · [Forces and pressure](forces-and-pressure.png) · [Methods](METHODS.md) · [Data](data/)

Both selected continuations cover **4–8 s, 800 new steps and 801 native states per solver**. Every new step meets the unchanged declared residual and boundary-mass limits. The last six complete force cycles satisfy the predeclared 0.5% range/mean stability criterion for period, sampled lift peak-to-peak and time-mean drag. This improves the numerical evidence in the existing cylinder showcase; it does not add a new numbered model.

Alternating vortices produce oscillating lift and a drag signal at roughly twice the lift frequency. Earlier pictures looked periodic while the lift amplitude was still growing. This update checks the iteration convergence and then measures how successive cycles change.

| Last-six-cycle quantity | ANSYS Fluent | OpenFOAM |
| --- | ---: | ---: |
'''
for label,k,fmt in [('Mean period (s)','mean_period_s','.6f'),('Strouhal number','St','.6f'),('Time-mean drag coefficient','time_mean_Cd','.6f'),('Mean sampled lift peak-to-peak','mean_sampled_Cl_peak_to_peak','.6f'),('Mean estimated front–rear wall pressure (Pa)','last_six_mean_wall_pressure_difference_Pa','.6f'),('Maximum relative boundary mass imbalance over 4–8 s','maximum_relative_boundary_mass_imbalance','.3e')]:readme+=f'| {label} | {f[k]:{fmt}} | {o[k]:{fmt}} |\n'
readme+='''
The lift-amplitude relative ranges over the last six cycles are **0.0336% in Fluent and 0.0279% in OpenFOAM**, below the declared 0.5% diagnostic criterion. This is observed window stability, not a proof of exact periodicity, discretization independence or physical accuracy. Native forces are plotted at their own times without phase alignment or fitted extrema.

![Actual numerical histories](numerical-checks.png)

## What changed

The original PISO Fluent calculation and both new pressure-tolerance-only diagnostics retain their numerical misses. Tightening pressure AMG termination alone did not remove the continuity plateau. A same-checkpoint Coupled algorithm bundle (flow Courant20, explicit pressure/momentum relaxation0.5) passes ten new steps; a guarded half-second window passes, followed by the finite continuation to8s. Physics, mesh, dt, spatial/time schemes and residual definitions remain fixed. This bundle does not uniquely isolate one cause.

OpenFOAM retains its eight PIMPLE outer corrections, two pressure corrections and original physical/numerical settings. A100-step qualification passes before its exact4.5s partitioned checkpoint is copied for the remaining700steps. All restart interfaces and duplicate states/force samples are checked. These are saved-state continuations, not fresh independent calculations from zero.

## Interpretation and limits

The inherited0–4s startup contains retained numerical misses. The4s initial fields are finite and checked, but their earlier residual failure is not retroactively repaired. Mesh and time-step independence are open. The128-segment cylinder and practical pressure outlets are not asserted equivalent to the published benchmark's exact geometry and traction condition; the modern reference also describes a longer development procedure. **No matched benchmark pass or experimental validation is claimed.**

Pressure differences are documented wall-vertex estimates at(0.15,0.2) and(0.25,0.2)m. They use adjacent native wall-face pressures interpolated in polar angle, with matching stencils and weights recorded. Fluent uses native face SV_P; OpenFOAM's zeroGradient evaluated face value equals its owner-cell pressure. Density1 converts its kinematic pressure numerically toPa. These are discretized estimates, not exact continuum point measurements.

The MP4 contains161actual saved states0.025s apart at20fps (8.05s playback). The compactGIF selects every second native frame:81states0.05s apart at20fps (4.05s playback). Every encoded frame was decoded and compared with its source. There is no generated intermediate motion or phase shift. The GIF reset is not simulated flow. Full161frameGIF remains external because of its41.3MB size; full-size PNG/PDF and the4.03MB MP4 remain available here.

Speed scales are fixed0–3m/s. Signed vorticity is an area-weighted reconstruction from native velocity, not a solver-exported exact wall derivative. Its colors clip at±60/s with extension markers; raw extrema and saturated fractions are retained. Red/blue show rotation sign, not flow direction. The force/pressure phase differences are retained and do not rank the solvers' accuracy.

[Reproduction and storage](REPRODUCING.md) · [Methods and primary sources](METHODS.md) · [Preserved history](history/) · [Source provenance](provenance/)

The owner directed the progression and requested paired visuals. Codex prepared, operated, diagnosed, checked, plotted and documented the work; Fluent and OpenFOAM computed their respective fields. Independent owner execution, reflection, research novelty or an admissions claim has not been invented. This is a local GitHub-ready record; no new upload is claimed. Fresh technical review is pending.
'''
(stage/'README.md').write_text(readme,encoding='utf-8')
(stage/'METHODS.md').write_text('''# Conditions, numerical checks and primary sources

The2.2 by0.41m channel contains a nominal0.1m cylinder centred at(0.2,0.2)m, represented by128linear segments. The shared mesh contains21,598triangles. Fluent is planar/unit depth; OpenFOAM is a1m empty-direction extrusion of the same topology. Positive cell areas, no bridge across the body, native mappings and mesh checks are verified. Fluent minimum orthogonal quality is about0.3333; the original OpenFOAM maximum nonorthogonality is about30degrees. These are mesh-quality checks, not mesh independence.

Density is1kg/m3 and dynamic/kinematic viscosity numerically0.001 in SI units. The model is incompressible, isothermal and laminar, nominalRe=100. The inlet is ux=6*y*(0.41-y)/0.41^2, uy=0: analytic mean1m/s, maximum1.5m/s. The41sampled inlet faces give0.4101219512m3/s,0.02974% above the analytic0.41m3/s. No inlet rescaling is applied. Cylinder and channel walls are stationary/no slip; practical zero gauge-pressure outlets are retained. Force references are rho1,U1,D0.1,area0.1m2 and unit depth.

Fixed dt=0.005s, second-order implicit time and second-order spatial schemes are preserved. Fluent uses Coupled pressure/velocity, numerical flow Courant20, explicit pressure/momentum relaxation0.5,120-inner cap, pressure second order and second-order-upwind momentum. Native scaled continuity/Ux/Uy must each be<=1e-5. Residual monitor definitions and normalization are retained from the parent. Each guarded suffix step is saved, then native HDF clock/residual/finite field/material and native mass checks are applied before another advance.

OpenFOAM Foundation14 incompressibleFluid uses backward time, linearUpwind momentum, eight fixed outer/two pressure/one nonorthogonal corrections. The unchanged residual metric is the last printed initial and final linear-system residual for each Ux/Uy/p within a physical step, after the fixed outer corrections. Every new step must meet1e-5. These residual normalizations differ from Fluent and cannot rank accuracy. Mesh/native field clocks for both ranks, finite U/p/phi, force and report coverage are checked. Raw Courant histories retain native pre-advance state times and both initial reports at each restart:802raw reports over800unique pre-advance times. No new arbitrary CFL limit is imposed on this implicit cylinder calculation.

Boundary imbalance is abs(sum of native signed boundary flux)/abs(inlet flux), limit1e-4 at every saved state. Fluent's report uses its native sign convention; OpenFOAM phi is outward-positive and includes inlet/outlet/walls/cylinder. Density and depth are1, so their mass/volume flux values are numerically consistent. All800newsteps and801field states per code are source-bound. Duplicate initial states at segment boundaries are compared before deduplication; initial4s residual convergence is not claimed.

Cycle crossings are successive rising lift zeros, linearly estimated between native force samples. Extrema use sampled values only. Cycle Cd means use trapezoidal integration with linear force interpolation only at cycle boundaries. The last-six criterion (range/absolute mean<=0.5% separately for period, lift peak-to-peak and time-meanCd) was recorded before the8s run completed. It is a declared diagnostic stability criterion, not uncertainty or benchmark acceptance.

Pressure estimates use two adjacent evaluated native wall-face pressures bracketing each stated vertex. Matching face-centre coordinates, angular offsets and weights are in data/wall-pressure-audit.json. The [Foundation14 zeroGradient implementation](https://github.com/OpenFOAM/OpenFOAM-14/blob/master/src/finiteVolume/fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.C) evaluates the wall value from the internal patch value. The selected case's native boundary is verified zeroGradient. Time-mean pressure differences use the same actual cycle boundaries and disclosed temporal endpoint interpolation.

The [primary Featflow periodic Re100 specification](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark2_re100.html) defines the geometry, inlet, traction outlet, pressure points and developed-cycle measurements. Our practical outlets are not proven equivalent; polygon geometry, startup, development duration and discretization differ. The earlier literature-frequency comparison remains a retained startup record, not a newly passed benchmark. The [ANSYS pressure–velocity coupling theory](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/flu_th/flu_th_sec_uns_solve_pv.html) supplies algorithm context; the numerical correction here is established by actual controlled runs, not by that documentation alone.
''',encoding='utf-8')
(stage/'REPRODUCING.md').write_text('''# Reproduction and storage

This is a methods/execution archive with small inputs and derived evidence. It is not a complete native restart bundle or a fresh-machine installation package. Large meshes, all HDF5/processor fields, full native logs, environments, licence/connection information and the41.3MB fullGIF remain outside Git on E:/CFD-Learning-Lab. Original projects and failed attempts are preserved.

Use the source-bound manifests in data/ and provenance/ to identify actual native cases. The code copies retain original E: paths and reject reused output directories. They require existing repository/scripts helpers and the preserved native checkpoints/input receipts; copying code alone does not satisfy those prerequisites. Do not run completed controllers against existing outputs.

On the original workstation, dot-source E:/CFD-Learning-Lab/local/Use-Project.ps1, use environments/fluent/Scripts/python.exe -B, and put local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001 plus repo/cfd-learning-lab/scripts on PYTHONPATH. Fluent remains on D:; Foundation14 is the existing CFD-Foundation14 WSL distribution. Native restart reproduction requires the exact hash-bound original4s and subsequent partitioned checkpoints, including backward time history.

Executed Fluent correction drivers/guards and all-step auditors, OpenFOAM preparation/gates/LFwrappers/auditor, combination/cycle/pressure/plot/export sources and regression tests are archived in code/. Earlier preparation drafts are retained and distinguished by the history ledger; reuse the selected versions identified in provenance. The small OpenFOAM dictionaries start at4.5s and need the external checkpoint; the0/ fields only document the original initialization. No mesh or installed dependency is bundled here.

Recheck hashes, active native processes, RAM>=4GiB and E:free>=10GiB after imports and immediately before any fresh launch. Use one hidden session at a time, two ranks, save before closing, and fresh versioned directories. The actual resource and source gates remain in the archived scripts/receipts. OMPI uses a project-specific/dev/shm directory for its unsupported E:FIFO; other project temp/cache files use E:.

To inspect without a solver, use the PNG/PDF/MP4/GIF/CSV/JSON files. Rendering again requires the selected external native fields and a fresh E:exports directory. Native velocity supplies speed and the documented spatial vorticity reconstruction. Exporters check source PNG hashes and decode every selected frame; full/compact exports have separate frame/time records.
''',encoding='utf-8')
(stage/'history/README.md').write_text('''# Preserved execution history

Original cylinder startup0–2s and earlier2–4s continuations retain their misses and source files in the earlier study records and external native archive. New pressure-tolerance-only diagnostics from the same4s state each miss continuity in all10steps; the tighter pressure solve alone does not cure the plateau. Their audits and source-bound hypothesis record are retained.

The accepted Coupled bundle first passes10steps, then100guarded steps, then690more to8s. Original spatial/time schemes, physics and residual definitions are checked against native readbacks. OpenFOAM first passes100steps4..4.5, then700remainingsteps from the byte-verified partitioned checkpoint. No limit is relaxed. These continuations inherit the earlier startup; no failed parent is relabelled.

Unlaunched preparation drafts are retained. Preflight found an inherited force-CSV schema mismatch and a native clock helper schema mismatch, each reproduced/corrected before a native solve. One controller stopped before case creation because WSL interpreted a pipe pattern via its shell and a localized UTF16 warning broke UTF8 decoding. Direct-exec WSL/bytes handling corrected the gate without advancing any solver. Native coverage/runtime/cycle/wall-stencil regressions and the ledger retain the evidence.

Plot versions with clipped values or overlapping legends/annotations are preserved outside Git; selected version002 plots were inspected. All original solver cases, source-bound scripts, archives and workspace files remain on E:. Nothing here claims independent owner execution or a new publication.
''',encoding='utf-8')
(stage/'copy-verification.json').write_text(json.dumps(dict(copies=copies,all_copied_files_byte_identical=True,native_solver_files_included=False,technical_review_pending=True),indent=2)+'\n')
# Standalone packet links must resolve before modifying the checkout.
for p in stage.rglob('*.md'):
 for target in re.findall(r'!?\[[^\]]*\]\(([^)\n]+)\)',p.read_text()):
  if not target.startswith(('https://','http://','#')):assert (p.parent/target).exists(),(p,target)
changed_names=['README.md','docs/PROJECT_HISTORY.md','docs/REPRODUCING.md','lessons/07-vortex-shedding/README.md','lessons/20-cylinder-benchmark/README.md','manifest-sha256.json'];backup=S/'record-root-backups-v001';assert not backup.exists()
for n in changed_names:d=backup/n;d.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(REPO/n,d)
shutil.copytree(stage,dest)
for p in stage.rglob('*'):
 if p.is_file():assert sha(p)==sha(dest/p.relative_to(stage))
root=(REPO/'README.md').read_text();marker='## Latest visual model: dam break over an obstacle';assert marker in root;root=root.replace(marker,'## Latest showcase improvement: cylinder wake\n\nThe [versioned cylinder update](lessons/20-cylinder-benchmark/convergence-v001/README.md) adds800numerically accepted new steps per code, actual paired video, developed-cycle evidence and transparent wall-pressure estimates. Earlier startup misses and benchmark/validation limits remain.\n\n![Actual paired cylinder continuation](lessons/20-cylinder-benchmark/convergence-v001/cylinder-compact.gif)\n\n'+marker,1);(REPO/'README.md').write_text(root,encoding='utf-8')
for lesson in ['07-vortex-shedding','20-cylinder-benchmark']:
 p=REPO/'lessons'/lesson/'README.md';old=p.read_text();head,rest=old.split('\n',1);link='../20-cylinder-benchmark/convergence-v001/README.md' if lesson.startswith('07') else 'convergence-v001/README.md';head='# Cylinder wake: convergence and developed-cycle evidence' if lesson.startswith('20') else head;p.write_text(head+'\n\n**8November2026 update:** [A new versioned continuation record]('+link+') covers4–8s with800accepted new steps/801nativefields per code, actual paired media, cycle stability and wall-pressure estimates. The inherited startup and unmatched benchmark conditions remain explicit.\n\n## Earlier startup record, retained\n'+rest,encoding='utf-8')
with (REPO/'docs/PROJECT_HISTORY.md').open('a',encoding='utf-8') as f:f.write('\n## Cylinder showcase convergence repair — 8November2026\n\nStudies07/20 gain a versioned4–8s paired continuation:800newsteps/801nativefields per solver, unchanged numerical limits, predeclared last-six-cycle stability, native forces, explicit wall-face pressure interpolation and source-checked media. Original misses/startup are preserved; no matched benchmark pass, physical validation or new upload is claimed. The model count remains35records.\n')
with (REPO/'docs/REPRODUCING.md').open('a',encoding='utf-8') as f:f.write('\n## Cylinder convergence and developed cycles\n\nSee lessons/20-cylinder-benchmark/convergence-v001 for the selected saved-state continuations, source/audit bindings, methods and actual media. Exact native4s/subsequent partitioned checkpoints, environments and meshes remain external on E:. This is an execution archive, not a clean-machine restart bundle. Never rerun frozen controllers over completed output paths.\n')
manifest={p.relative_to(REPO).as_posix():sha(p) for p in sorted(REPO.rglob('*')) if p.is_file() and '.git' not in p.relative_to(REPO).parts and p.name!='manifest-sha256.json'};(REPO/'manifest-sha256.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
changed=[n for n,h in base.items() if sha(REPO/n)!=h];assert set(changed)<=set(changed_names),changed
result=dict(record=str(dest),record_files=sum(p.is_file() for p in dest.rglob('*')),baseline_files=len(base),changed_existing_files=changed,manifest_entries=len(manifest),technical_review_pending=True,github_published=False);(S/'record-preparation-v001.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))