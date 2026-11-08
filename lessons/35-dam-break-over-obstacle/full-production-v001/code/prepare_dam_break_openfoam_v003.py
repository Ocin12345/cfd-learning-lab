"""Prepare fresh Foundation14 dam-break inputs from shared analytic geometry."""
from pathlib import Path
import argparse,json
import numpy as np
from dam_break_common_v002 import *
R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/dam-break-over-obstacle-v001'
def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='ascii',newline='\n') as f:f.write(text)
def header(name,kind='dictionary',location='system'):
    return f'FoamFile\n{{ format ascii; class {kind}; location "{location}"; object {name}; }}\n'
def field(out,name,values,kind,dimensions,boundaries):
    scalar=values.ndim==1;typ='scalar' if scalar else 'vector'
    rows=[f'{v:.17g}' if scalar else '('+' '.join(f'{x:.17g}' for x in v)+')' for v in values]
    text=header(name,kind,'0')+f'dimensions {dimensions};\ninternalField nonuniform List<{typ}>\n{len(rows)}\n(\n'+'\n'.join(rows)+'\n);\nboundaryField\n{\n'+boundaries+'\n}\n'
    write(out/'0'/name,text)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['flat','collapse','production']);ap.add_argument('out',type=Path);args=ap.parse_args();out=args.out;assert out.resolve().is_relative_to(N.resolve()) and not out.exists();out.mkdir(parents=True)
    g=geometry();init=initial_fields(g,args.mode);order=g['openfoam_to_common'];steps=10000 if args.mode=='production' else 400;interval=100 if args.mode=='production' else 20
    xy=[(0,0),(OB_LEFT,0),(OB_RIGHT,0),(L,0),(0,OB_HEIGHT),(OB_LEFT,OB_HEIGHT),(OB_RIGHT,OB_HEIGHT),(L,OB_HEIGHT),(0,H),(OB_LEFT,H),(OB_RIGHT,H),(L,H)]
    vertices=[(x,y,z) for z in [0,DEPTH] for x,y in xy]
    v='\n'.join('('+' '.join(f'{a:.17g}' for a in row)+')' for row in vertices)
    blocks='''
    hex (0 1 5 4 12 13 17 16) (46 16 1) simpleGrading (1 1 1)
    hex (2 3 7 6 14 15 19 18) (38 16 1) simpleGrading (1 1 1)
    hex (4 5 9 8 16 17 21 20) (46 84 1) simpleGrading (1 1 1)
    hex (5 6 10 9 17 18 22 21) (8 84 1) simpleGrading (1 1 1)
    hex (6 7 11 10 18 19 23 22) (38 84 1) simpleGrading (1 1 1)
'''
    mesh=header('blockMeshDict')+'scale 1;\nvertices\n(\n'+v+'\n);\nblocks\n(\n'+blocks+');\ndefaultPatch { type empty; }\nboundary\n(\nwalls { type wall; faces ((0 12 16 4)(4 16 20 8)(7 19 15 3)(11 23 19 7)(0 1 13 12)(1 5 17 13)(5 6 18 17)(2 14 18 6)(2 3 15 14)); }\natmosphere { type patch; faces ((8 20 21 9)(9 21 22 10)(10 22 23 11)); }\n);\n'
    write(out/'system/blockMeshDict',mesh)
    common='defaultFaces { type empty; }\n'
    field(out,'alpha.water',init['alpha'][order],'volScalarField','[0 0 0 0 0 0 0]',common+'walls { type zeroGradient; }\natmosphere { type inletOutlet; inletValue uniform 0; value uniform 0; }')
    field(out,'U',init['velocity'][order],'volVectorField','[0 1 -1 0 0 0 0]',common+'walls { type fixedValue; value uniform (0 0 0); }\natmosphere { type pressureInletOutletVelocity; value uniform (0 0 0); }')
    field(out,'p_rgh',init['p_rgh'][order],'volScalarField','[1 -1 -2 0 0 0 0]',common+f'walls {{ type fixedFluxPressure; value uniform 0; }}\natmosphere {{ type prghPressure; p uniform 0; value uniform {G*H:.17g}; }}')
    write(out/'constant/g',header('g','uniformDimensionedVectorField','constant')+'dimensions [0 1 -2 0 0 0 0]; value (0 -9.81 0);\n')
    write(out/'constant/phaseProperties',header('phaseProperties',location='constant')+'phases (water air); sigma 0.07;\n')
    for name,(rho,mu) in [('water',MATERIALS['water-liquid']),('air',MATERIALS['air'])]:write(out/('constant/physicalProperties.'+name),header('physicalProperties.'+name,location='constant')+f'viscosityModel constant; nu {mu/rho:.17g}; rho {rho:.17g};\n')
    write(out/'constant/momentumTransport',header('momentumTransport',location='constant')+'simulationType laminar;\n')
    functions='''
functions
{
    atmosphereFlux
    {
        type surfaceFieldValue;
        libs ("libfieldFunctionObjects.so");
        patch atmosphere;
        operation orientedSum;
        fields (rhoPhi phi);
        writeFields false;
        log true;
        executeControl timeStep;
        executeInterval 1;
        writeControl timeStep;
        writeInterval 1;
    }
}
'''
    write(out/'system/controlDict',header('controlDict')+f'application foamRun; solver incompressibleVoF; startFrom startTime; startTime 0; stopAt endTime; endTime {steps*DT:.12g}; deltaT {DT:.12g}; writeControl timeStep; writeInterval {interval}; purgeWrite 0; writeFormat ascii; writePrecision 15; writeCompression off; timeFormat general; timePrecision 12; runTimeModifiable false; adjustTimeStep no; maxCo 0.25; maxAlphaCo 0.25;\n'+functions)
    schemes='''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes
{
    div(rhoPhi,U) Gauss linearUpwind grad(U);
    div(phi,alpha) Gauss interfaceCompression vanLeer 1;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default no; p_rgh; alpha.water; }
'''
    write(out/'system/fvSchemes',header('fvSchemes')+schemes)
    solution='''solvers
{
    "alpha.water.*" { nCorrectors 2; nSubCycles 1; MULESCorr yes; MULES { nIter 10; tolerance 0.01; } solver smoothSolver; smoother symGaussSeidel; tolerance 1e-11; relTol 0; }
    "pcorr.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
    p_rgh { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0.05; }
    p_rghFinal { $p_rgh; relTol 0; }
    "(U|k|epsilon|omega).*" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-6; relTol 0; minIter 1; }
}
PIMPLE { momentumPredictor yes; nOuterCorrectors 8; nCorrectors 3; nNonOrthogonalCorrectors 0; }
relaxationFactors { equations { ".*" 1; } }
'''
    write(out/'system/fvSolution',header('fvSolution')+solution)
    write(out/'system/decomposeParDict',header('decomposeParDict')+'numberOfSubdomains 2; method simple; simpleCoeffs { n (2 1 1); delta 0.001; }\n')
    declared=dict(mode=args.mode,cells=NCELLS,domain_m=[L,H],obstacle_m=[OB_LEFT,OB_RIGHT,OB_HEIGHT],depth_m=DEPTH,dt_s=DT,steps=steps,save_interval=interval,water_area_m2=init['water_area_m2'],mixture_mass_per_m=init['mixture_mass_per_m'],static_atmosphere_gauge_Pa=0,backflow_water_fraction=0,surface_tension_N_m=SIGMA,wall_adhesion=False,materials=MATERIALS,initial_pressure='shared hydrostatic physical gauge; converted to p_rgh',all_native_flux_accounting_pending=True,physical_validation=False)
    write(out/'declared-settings.json',json.dumps(declared,indent=2)+'\n');print(json.dumps(dict(case=str(out),mode=args.mode,cells=NCELLS,steps=steps,water_area_m2=init['water_area_m2'])))
if __name__=='__main__':main()