"""Prepare a matched Re=20 cylinder case for OpenFOAM Foundation 14.

Mesh conversion and checkMesh run separately; this never starts a solver.
"""
from pathlib import Path
import argparse,hashlib,json,shutil

p=argparse.ArgumentParser()
p.add_argument('--mesh',required=True,type=Path)
p.add_argument('--case',required=True,type=Path)
a=p.parse_args()
repo=Path(__file__).resolve().parents[1]
case=a.case.resolve()
if case.exists() or case.is_relative_to(repo):raise ValueError('Use a fresh external case directory')
case.mkdir(parents=True)

def write(relative,body,kind='dictionary',dimensions=None):
    path=case/relative
    path.parent.mkdir(parents=True,exist_ok=True)
    header=f'FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {kind};\n    object {path.name};\n}}\n\n'
    path.write_bytes((header+body.strip()+'\n').encode())

write('system/controlDict','''
solver incompressibleFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 1200;
deltaT 1;
writeControl timeStep;
writeInterval 1200;
purgeWrite 0;
writeFormat ascii;
writePrecision 12;
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
        patches (cylinder);
        rho rhoInf;
        rhoInf 1;
        p p;
        U U;
        CofR (0 0 0);
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
''')
write('system/fvSchemes','''
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
''')
write('system/fvSolution','''
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
''')
write('system/decomposeParDict','''
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
''')
write('constant/physicalProperties','viscosityModel constant;\nnu 0.05;')
write('constant/momentumTransport','simulationType laminar;')
write('0/U','''
dimensions [0 1 -1 0 0 0 0];
internalField uniform (1 0 0);
boundaryField
{
    inlet { type fixedValue; value uniform (1 0 0); }
    outlet { type zeroGradient; }
    cylinder { type fixedValue; value uniform (0 0 0); }
    "(frontAndBack.*|defaultFaces)" { type empty; }
}
''','volVectorField')
write('0/p','''
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    cylinder { type zeroGradient; }
    "(frontAndBack.*|defaultFaces)" { type empty; }
}
''','volScalarField')
shutil.copy2(a.mesh,case/'input-cylinder.msh')
targets=dict(reynolds_number=20,density_kg_m3=1,kinematic_viscosity_m2_s=.05,
             diameter_m=1,velocity_m_s=1,span_m=1,reference_area_m2=1,
             mesh_source_sha256=hashlib.sha256(a.mesh.read_bytes()).hexdigest(),
             intended_cells=98304,outer_radius_m=20,
             solver='OpenFOAM Foundation 14 incompressibleFluid',
             solve='steady SIMPLE, least-squares gradients, linear-upwind momentum',
             target_iterations=1200,final_force_window=100,
             cd_relative_range_target=1e-4,relative_boundary_mass_imbalance_target=1e-6,
             initial_p_U_residual_learning_target=1e-7,
             cross_code_cd_difference_learning_target_percent=.5,
             residual_norms_equivalent_to_fluent=False,physical_validation=False)
(case/'declared-comparison.json').write_text(json.dumps(targets,indent=2))
print(json.dumps(dict(case=str(case),source_mesh_sha256=targets['mesh_source_sha256'],solver_started=False),indent=2))
