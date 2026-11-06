"""Prepare a matched zero-incidence airfoil case for OpenFOAM Foundation 14.

Writes dictionaries and copies the supplied mesh; never runs conversion, mesh
checking, MPI, WSL or a CFD solver. The first cap is 2000 steady iterations.
Two nonorthogonal pressure correctors and U equation relaxation 0.7 are declared
before solving to accommodate the repaired triangular mesh, unlike the nearly
orthogonal cylinder. These choices are not evidence of accuracy or convergence.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import shutil


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(mesh, case):
    mesh, case = mesh.resolve(), case.resolve()
    repo = Path(__file__).resolve().parents[1]
    require(not case.exists() and not case.is_relative_to(repo),
            'Use a fresh external native case directory')
    setup_path = mesh.parent / 'setup.json'
    setup = json.loads(setup_path.read_text())
    rho, mu = setup['density_kg_m3'], setup['viscosity_Pa_s']
    velocity, chord = setup['inlet_velocity_m_s'], setup['chord_m']
    require(all(math.isfinite(x) and x > 0 for x in [rho, mu, velocity, chord]),
            'Positive finite physical inputs required')
    nu = mu / rho
    require(setup['angle_degrees'] == 0 and math.isclose(velocity * chord / nu, 100),
            'This preparation supports the matched zero-incidence Re100 case')
    require(rho == velocity == chord == 1 and nu == .01,
            'Expected unit density/speed/chord and nu=0.01')
    text = mesh.read_text()
    match = re.search(r'^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)', text, re.M)
    require(match is not None, 'Unsupported Fluent ASCII cell header')
    cells = int(match[1], 16)
    require(cells == setup['cells'] and setup['mesh_type'] == 'triangles',
            'Source mesh cell count/type differs from setup')
    targets = dict(reynolds_number=velocity * chord / nu, angle_degrees=0,
                   density_kg_m3=rho, dynamic_viscosity_Pa_s=mu,
                   kinematic_viscosity_m2_s=nu, velocity_m_s=velocity,
                   chord_m=chord, span_m=1, reference_area_m2=chord,
                   pressure_reference_Pa_gauge=0, wall_patch='airfoil',
                   mesh_source_path=str(mesh), mesh_source_sha256=sha256(mesh),
                   mesh_setup_sha256=sha256(setup_path), intended_cells=cells,
                   intended_nodes=setup['nodes'], outer_radius_m=setup['outer_radius_m'],
                   mesh_type=setup['mesh_type'],
                   solver='OpenFOAM Foundation 14 incompressibleFluid',
                   solve='steady SIMPLE consistent; leastSquares gradients; bounded linearUpwind momentum; corrected diffusion',
                   mpi_processes=2, initial_iteration_cap=2000,
                   n_nonorthogonal_correctors=2, U_equation_relaxation=.7,
                   numerical_choices_note='Two nonorthogonal pressure correctors and U relaxation 0.7 declared for triangular mesh; no accuracy claim',
                   final_force_window=100, cd_relative_range_target=1e-4,
                   cl_absolute_range_target=1e-5, zero_angle_abs_lift_diagnostic_target=1e-3,
                   relative_boundary_mass_imbalance_target=1e-6,
                   initial_p_U_residual_learning_target=1e-7,
                   pressure_residual_definition='Maximum printed initial p residual across all pressure correction solves within the iteration',
                   cross_code_cd_difference_learning_target_percent=.5,
                   residual_norms_equivalent_to_fluent=False, physical_validation=False,
                   initialization='Uniform U=(1,0,0), kinematic p=0; Fluent reference uses its own saved converged field')
    case.mkdir(parents=True)

    def write(relative, body, kind='dictionary'):
        path = case / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        header = ('FoamFile\n{\n    version 2.0;\n    format ascii;\n'
                  f'    class {kind};\n    object {path.name};\n}}\n\n')
        with path.open('x', encoding='ascii', newline='\n') as stream:
            stream.write(header + body.strip() + '\n')

    write('system/controlDict', """
solver incompressibleFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 2000;
deltaT 1;
writeControl timeStep;
writeInterval 2000;
purgeWrite 0;
writeFormat ascii;
writePrecision 14;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable true;
functions
{
    forceCoefficients
    {
        type forceCoeffs;
        libs ("libforces.so");
        patches (airfoil);
        rho rhoInf;
        rhoInf 1;
        p p;
        U U;
        CofR (0.25 0 0);
        liftDir (0 1 0);
        dragDir (1 0 0);
        pitchAxis (0 0 1);
        magUInf 1;
        lRef 1;
        Aref 1;
        log false;
        writeControl timeStep;
        writeInterval 1;
    }
}
""")
    write('system/fvSchemes', """
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
    write('system/fvSolution', """
solvers
{
    p { solver GAMG; tolerance 1e-10; relTol 0.01; smoother GaussSeidel; }
    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01; }
}
SIMPLE { nNonOrthogonalCorrectors 2; consistent yes; }
relaxationFactors { equations { U 0.7; } }
""")
    write('system/decomposeParDict', """
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
""")
    write('constant/physicalProperties', 'viscosityModel constant;\nnu 0.01;')
    write('constant/momentumTransport', 'simulationType laminar;')
    write('0/U', """
dimensions [0 1 -1 0 0 0 0];
internalField uniform (1 0 0);
boundaryField
{
    inlet { type fixedValue; value uniform (1 0 0); }
    outlet { type zeroGradient; }
    airfoil { type fixedValue; value uniform (0 0 0); }
    "(frontAndBack.*|defaultFaces)" { type empty; }
}
""", 'volVectorField')
    write('0/p', """
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    airfoil { type zeroGradient; }
    "(frontAndBack.*|defaultFaces)" { type empty; }
}
""", 'volScalarField')
    shutil.copy2(mesh, case / 'input-airfoil.msh')
    shutil.copy2(setup_path, case / 'source-setup.json')
    shutil.copy2(Path(__file__), case / 'executed-prepare-script.py')
    with (case / 'declared-comparison.json').open('x', encoding='utf-8') as stream:
        json.dump(targets, stream, indent=2, allow_nan=False)
        stream.write('\n')
    with (case / 'preparation-provenance.json').open('x', encoding='utf-8') as stream:
        json.dump(dict(script_sha256=sha256(Path(__file__)),
                       mesh_source_sha256=targets['mesh_source_sha256'],
                       mesh_setup_sha256=targets['mesh_setup_sha256'],
                       conversion_started=False, mesh_check_started=False, solver_started=False), stream, indent=2)
        stream.write('\n')
    with (case / 'run-status.json').open('x', encoding='utf-8') as stream:
        json.dump(dict(stage='prepared_only', solver_started=False), stream, indent=2)
        stream.write('\n')
    return targets


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--case', type=Path, required=True)
    args = parser.parse_args()
    targets = prepare(args.mesh, args.case)
    print(json.dumps(dict(case=str(args.case.resolve()), cells=targets['intended_cells'],
                          mesh_source_sha256=targets['mesh_source_sha256'], solver_started=False), indent=2))


if __name__ == '__main__':
    main()
