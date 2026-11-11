# CFD model execution workflow

The owner clarified that actual CFD model runs and results take priority over
lessons, and selected continued airfoil simulations. Keep explanations brief,
show actual solver-derived figures after each completed calculation and
continue within the authorized model sequence until told to stop. Formal
prediction/reflection checkpoints are no longer required. Do not invent owner
participation, independent work, research novelty or validation.

For every model, preserve geometry provenance, assumptions, boundary conditions,
coefficient definitions, numerical methods, source code, native solver evidence,
checks, figures and limitations. Numerical convergence, mesh/domain sensitivity,
cross-code agreement and physical validation are distinct claims.

## Model sequence

Continue airfoil mesh sensitivity and matched Fluent/OpenFOAM calculations,
then airfoil domain sensitivity and a small incidence comparison at declared
conditions. Later models include transient cylinder shedding with time-step
checks, a genuinely matching published benchmark and a simple 3D wing.
Concorde and the ornithopter retain their own geometry, mesh and physics
prerequisites; smaller-model results do not validate those aircraft. This is
a model roadmap, not an instruction to launch every case concurrently.

## Execution and storage

Use `E:\CFD-Learning-Lab` as the root. Repository content stays in
`repo/cfd-learning-lab`; imported sources and large native runs stay outside
Git under `native`, Python under `environments`, supported caches/temp files
under `temp`, and large/vector exports under `exports`. Use installed ANSYS
on D:. Copy and verify required sources while preserving originals.

Use one hidden Fluent session at a time, initially with two solver processes.
Start OpenFOAM with two MPI processes. Check disk, RAM, licence availability
and existing processes before runs. Admit overlapping cases only when reserve
RAM is adequate. Save native results before closing. Use fresh versioned paths,
retain failed/rejected trials and never weaken criteria to manufacture a pass.
Common-Licensing feature queries may not enumerate the Student entitlement;
a successful live Student launch is recorded separately.

Match geometry, meshes, flow conditions and force references across codes where
possible; document discretization, initialization and residual-definition
differences. Iteration counts are not a performance comparison with different
initializations. Follow `FIGURE_STANDARD.md`, use real numerical data and retain
captions, source hashes, vector axes and rendered export checks.

## Contributions and publication

The owner selects models and priorities and contributes any work actually
performed. Gmsh computes meshes; Fluent and OpenFOAM compute their numerical solutions. Keep credentials, connection/licence files, installed dependencies and private
application records outside Git. Inspect the complete local contents before
publication.

Older lesson folders remain as history of the project. They do not override the
owner's newer preference for simulation-led work without lesson checkpoints.
