"""Prepare the fresh 8,800-cell backward-step case for OpenFOAM Foundation 14.

Writes inputs/provenance only; conversion, checks, and solver are separate.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil


def prepare(mesh, case):
    mesh, case = mesh.resolve(), case.resolve()
    repo = Path(__file__).resolve().parents[1]
    lab = repo.parent.parent
    expected = lab / "native/runs/2026-11-01-backward-step-v001/openfoam"
    if case != expected.resolve() or case.exists():
        raise ValueError("Use the fresh native backward-step-v001/openfoam directory")
    content = mesh.read_bytes()
    cells = re.search(rb"^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)", content, re.M)
    if cells is None or int(cells[1], 16) != 8800:
        raise ValueError("Expected the shared 8,800-cell Fluent ASCII step mesh")
    if shutil.disk_usage(case.anchor).free < 10 * 1024**3:
        raise RuntimeError("Require 10 GiB free output space")
    local = lab / "local/backward-step-v001"
    for source in (local / "DESIGN.md", local / "run-openfoam.sh"):
        if not source.is_file():
            raise FileNotFoundError(source)
    case.mkdir(parents=True)

    def write(relative, body, kind="dictionary"):
        path = case / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        header = ("FoamFile\n{\n    version 2.0;\n    format ascii;\n"
                  f"    class {kind};\n    object {path.name};\n}}\n\n")
        with path.open("x", encoding="ascii", newline="\n") as stream:
            stream.write(header + body.strip() + "\n")

    write("system/controlDict", """
solver incompressibleFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 2000;
deltaT 1;
writeControl timeStep;
writeInterval 100;
purgeWrite 0;
writeFormat ascii;
writePrecision 14;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable false;
""")
    write("system/fvSchemes", """
ddtSchemes { default steadyState; }
gradSchemes { default leastSquares; }
divSchemes
{
    default none;
    div(phi,U) bounded Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    write("system/fvSolution", """
solvers
{
    p
    {
        solver GAMG;
        tolerance 1e-10;
        relTol 0.01;
        smoother GaussSeidel;
    }
    U
    {
        solver smoothSolver;
        smoother symGaussSeidel;
        tolerance 1e-10;
        relTol 0.01;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent yes;
}
relaxationFactors
{
    equations { U 0.9; }
}
""")
    write("system/decomposeParDict", """
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
""")
    write("constant/physicalProperties", "viscosityModel constant;\nnu 0.001;")
    write("constant/momentumTransport", "simulationType laminar;")
    write("0/U", """
dimensions [0 1 -1 0 0 0 0];
internalField uniform (1 0 0);
boundaryField
{
    inlet { type fixedValue; value uniform (1 0 0); }
    outlet { type zeroGradient; }
    walls { type noSlip; }
    "(frontAndBack.*|defaultFaces)" { type empty; }
}
""", "volVectorField")
    write("0/p", """
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    walls { type zeroGradient; }
    "(frontAndBack.*|defaultFaces)" { type empty; }
}
""", "volScalarField")
    shutil.copy2(mesh, case / "input-step.msh")
    provenance = case / "provenance"
    provenance.mkdir()
    for source in (Path(__file__), local / "DESIGN.md", local / "run-openfoam.sh"):
        shutil.copy2(source, provenance / source.name)
    setup = mesh.parent / "setup.json"
    if setup.is_file():
        shutil.copy2(setup, provenance / "mesh-setup.json")
    targets = dict(
        model="2D backward-facing step; steady incompressible laminar",
        x_bounds_m=[-.4, 2], y_bounds_m=[0, .2], step_height_m=.1,
        upstream_floor_m=.1, expansion_ratio=2, extrusion_depth_m=1,
        inlet_velocity_m_s=1, density_kg_m3=1,
        dynamic_viscosity_Pa_s=.001, kinematic_viscosity_m2_s=.001,
        reynolds_step_height=100, upstream_reynolds_hydraulic_diameter=200,
        inlet_note="Uniform inlet develops over the explicit upstream length",
        outlet_gauge_pressure_Pa=0, grid_dx_m=.01, grid_dy_m=.005,
        intended_cells=8800, mesh_source_path=str(mesh),
        mesh_source_sha256=hashlib.sha256(content).hexdigest(),
        solver="OpenFOAM Foundation 14 incompressibleFluid",
        solve="steady SIMPLE consistent; leastSquares; bounded linearUpwind; corrected diffusion",
        mpi_processes=2, initial_iteration_cap=2000, write_interval=100,
        U_equation_relaxation=.9, linear_solver_absolute_tolerance=1e-10,
        residual_target=1e-6, relative_mass_imbalance_target=1e-5,
        final100_checkpoint_iterations=[1900, 2000],
        final100_relative_velocity_change_target=.001,
        final100_relative_developed_pressure_drop_change_target=.001,
        final100_note="Endpoint field difference at 1900/2000 is not an every-iteration stability test",
        residual_early_stop=False, physical_validation=False,
        residual_norms_equivalent_to_fluent=False,
        pressure_field_note="Kinematic p [m2/s2]; Pa gauge = density * p")
    (case / "declared-comparison.json").write_text(
        json.dumps(targets, indent=2) + "\n", encoding="ascii")
    hashes = {str(p.relative_to(case)).replace("\\", "/"):
              hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(case.rglob("*")) if p.is_file()}
    (provenance / "prepared-input-sha256.json").write_text(
        json.dumps(hashes, indent=2) + "\n", encoding="ascii")
    (case / "run-status.json").write_text(
        json.dumps(dict(stage="prepared", solver_started=False)) + "\n",
        encoding="ascii")
    print(json.dumps(dict(case=str(case), intended_cells=8800,
                          mesh_source_sha256=targets["mesh_source_sha256"],
                          solver_started=False), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--case", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.mesh, args.case)
