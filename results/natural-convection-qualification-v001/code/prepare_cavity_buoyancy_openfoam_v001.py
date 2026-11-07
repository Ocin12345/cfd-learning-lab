"""Prepare the native ideal-gas cavity control; no solver starts."""
from pathlib import Path
import hashlib,json,shutil,sys
from cavity_buoyancy_common_v002 import *
LOCAL=ROOT/'local/natural-convection-cavity-v001';CASE=ROOT/'native/runs/natural-convection-qualification-openfoam-v001/control'
def head(name,kind='dictionary'):return f'FoamFile\n{{ version 2.0; format ascii; class {kind}; object {name}; }}\n'
def save(name,text,kind='dictionary'):
    p=CASE/name;p.parent.mkdir(exist_ok=True,parents=True);p.write_text(head(p.name,kind)+text,encoding='ascii',newline='\n')
def main():
    assert not CASE.exists();CASE.mkdir(parents=True)
    ru=float((LOCAL/'constant-probe/actual-RU.txt').read_text());mw=ru/R
    spec=specification();spec.update(openfoam_universal_gas_constant=ru,openfoam_molecular_weight= mw,molecular_weight_adjustment_reason='Align specific R with installed Fluent constant; native compiled OpenFOAM constant probe retained')
    write(CASE/'declared-settings.json',spec)
    msh=ROOT/'native/meshes/natural-convection-qualification-v001/cavity.msh';shutil.copyfile(msh,CASE/'input-cavity.msh')
    save('system/controlDict','''solver fluid;
startFrom startTime; startTime 0; stopAt endTime; endTime 300; deltaT 1;
writeControl timeStep; writeInterval 100; purgeWrite 0; writeFormat ascii; writePrecision 16;
writeCompression off; timeFormat general; timePrecision 10; runTimeModifiable false;
''')
    save('system/fvSchemes','''ddtSchemes { default steadyState; }
gradSchemes { default leastSquares; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U); div(phi,h) bounded Gauss linearUpwind grad(h); div(phi,K) bounded Gauss linear; div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear orthogonal; }
interpolationSchemes { default linear; }
snGradSchemes { default orthogonal; }
''')
    save('system/fvSolution','''solvers {
 p_rgh { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0.01; }
 p_rghFinal { $p_rgh; relTol 0; }
 "(U|h)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0.01; }
 "(U|h)Final" { $U; relTol 0; }
 "rho.*" { solver diagonal; }
}
PIMPLE { momentumPredictor yes; nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 0; pRefCell 0; pRefValue 0; }
relaxationFactors { fields { p_rgh 0.3; } equations { U 0.7; h 0.7; } }
''')
    save('system/decomposeParDict','numberOfSubdomains 2; method simple; simpleCoeffs { n (2 1 1); delta 0.001; }\n')
    save('constant/physicalProperties',f'''thermoType {{ type heRhoThermo; mixture pureMixture; transport const; thermo hConst; equationOfState perfectGas; specie specie; energy sensibleEnthalpy; }}
mixture {{ specie {{ molWeight {mw:.17g}; }} thermodynamics {{ Cp {CP:.17g}; hf 0; }} transport {{ mu {MU:.17g}; Pr {PR:.17g}; }} }}
dpdt true;
''')
    save('constant/momentumTransport','simulationType laminar;\n')
    save('constant/g',f'dimensions [0 1 -2 0 0 0 0]; value (0 {-G} 0);\n','uniformDimensionedVectorField')
    save('constant/pRef',f'dimensions [1 -1 -2 0 0 0 0]; value {P0};\n','uniformDimensionedScalarField')
    save('constant/hRef',f'dimensions [0 1 0 0 0 0 0]; value {L/2};\n','uniformDimensionedScalarField')
    for name,dim,value,kind in [('U','[0 1 -1 0 0 0 0]','(0 0 0)','volVectorField'),('T','[0 0 0 1 0 0 0]',str(T0),'volScalarField'),('p','[1 -1 -2 0 0 0 0]',str(P0),'volScalarField'),('p_rgh','[1 -1 -2 0 0 0 0]','0','volScalarField')]:
        text=f'dimensions {dim};\ninternalField uniform {value};\nboundaryField\n{{\n'
        for wall in WALLS:
            bc={'U':'type noSlip;','T':f'type fixedValue; value uniform {T0};','p':f'type calculated; value uniform {P0};','p_rgh':'type fixedFluxPressure; value uniform 0;'}[name]
            text+=f' {wall} {{ {bc} }}\n'
        text+=' frontAndBackPlanes { type empty; }\n}\n';save('0/'+name,text,kind)
    shutil.copyfile(__file__,CASE/'executed-preparer.py')
    write(CASE/'prepared-input-hashes.json',{p.relative_to(CASE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in CASE.rglob('*') if p.is_file()})
    print(json.dumps(spec,indent=2))
if __name__=='__main__':main()
