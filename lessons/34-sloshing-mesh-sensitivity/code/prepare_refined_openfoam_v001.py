"""Prepare the matched OpenFOAM 14 sloshing qualification."""
from pathlib import Path
import shutil
import sys

from sloshing_refinement_common_v001 import NX, NY, DX, L, H, MATERIALS, WATER_WIDTH, WATER_HEIGHT, DT, LAST_STEP, FRAME_INTERVAL, INITIAL_TILT, MODE, write_json

DEPTH = 0.01

def main(out):
    if shutil.disk_usage(out.anchor).free < 10 * 2**30:
        raise RuntimeError("Need 10 GiB on E:")
    out.mkdir(parents=True)

    def write(rel, content, kind="dictionary"):
        path = out / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        header = f"FoamFile {{ version 2.0; format ascii; class {kind}; object {path.name}; }}\n"
        path.write_text(header + content.strip() + "\n", encoding="ascii", newline="\n")

    write("system/blockMeshDict", f"""
        scale 1;
        vertices
        (
            (0 0 0) ({L} 0 0) ({L} {H} 0) (0 {H} 0)
            (0 0 {DEPTH}) ({L} 0 {DEPTH}) ({L} {H} {DEPTH}) (0 {H} {DEPTH})
        );
        blocks (hex (0 1 2 3 4 5 6 7) ({NX} {NY} 1) simpleGrading (1 1 1));
        edges ();
        boundary
        (
            tankWalls
            {{
                type wall;
                faces ((0 1 5 4) (1 2 6 5) (2 3 7 6) (3 0 4 7));
            }}
            frontAndBack
            {{
                type empty;
                faces ((0 3 2 1) (4 5 6 7));
            }}
        );
        mergePatchPairs ();
    """)
    write("constant/g", "dimensions [acceleration];\nvalue (0 -9.81 0);", "uniformDimensionedVectorField")
    write("constant/momentumTransport", "simulationType laminar;")
    write("constant/phaseProperties", "phases (water air);\nsigma 0;")
    for phase, material in (("water", "water-liquid"), ("air", "air")):
        density, viscosity = MATERIALS[material]
        write(f"constant/physicalProperties.{phase}",
              f"viscosityModel constant;\nnu {viscosity / density:.16g};\nrho {density};")
    import numpy as np
    x_mid = (np.arange(NX) + 0.5) * (L / NX)
    y_bottom = np.arange(NY) * (H / NY)
    surface = 0.05 + INITIAL_TILT * (x_mid - L / 2)
    alpha = np.clip((surface[None, :] - y_bottom[:, None]) / (H / NY), 0, 1)
    assert abs(alpha.sum() * (L / NX) * (H / NY) - 0.01) < 1e-12
    values = "\n".join(f"{value:.17g}" for value in alpha.ravel())
    write("0/alpha.water", f"""
        dimensions [];
        internalField nonuniform List<scalar> {NX * NY}
        (
        {values}
        );
        boundaryField
        {{
            tankWalls {{ type zeroGradient; }}
            frontAndBack {{ type empty; }}
        }}
    """, "volScalarField")
    write("0/p_rgh", """
        dimensions [pressure];
        internalField uniform 0;
        boundaryField
        {
            tankWalls { type fixedFluxPressure; value uniform 0; }
            frontAndBack { type empty; }
        }
    """, "volScalarField")
    write("0/U", """
        dimensions [velocity];
        internalField uniform (0 0 0);
        boundaryField
        {
            tankWalls { type noSlip; }
            frontAndBack { type empty; }
        }
    """, "volVectorField")
    write("system/controlDict", f"""
        application foamRun;
        solver incompressibleVoF;
        startFrom startTime;
        startTime 0;
        stopAt endTime;
        endTime {LAST_STEP * DT};
        deltaT {DT};
        writeControl timeStep;
        writeInterval {FRAME_INTERVAL};
        purgeWrite 0;
        writeFormat ascii;
        writePrecision 12;
        writeCompression off;
        timeFormat general;
        timePrecision 8;
        runTimeModifiable false;
        adjustTimeStep false;
    """)
    write("system/fvSchemes", """
        ddtSchemes { default Euler; }
        gradSchemes { default Gauss linear; }
        divSchemes
        {
            div(phi,alpha) Gauss interfaceCompression vanLeer 1;
            div(rhoPhi,U) Gauss linearUpwind grad(U);
            div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
        }
        laplacianSchemes { default Gauss linear uncorrected; }
        interpolationSchemes { default linear; }
        snGradSchemes { default uncorrected; }
    """)
    write("system/fvSolution", """
        solvers
        {
            "alpha.water.*"
            {
                nCorrectors 2;
                nSubCycles 1;
                MULESCorr yes;
                MULES { nIter 10; tolerance 1e-2; }
                solver smoothSolver;
                smoother symGaussSeidel;
                tolerance 1e-11;
                relTol 0;
            }
            "pcorr.*"
            {
                solver PCG;
                preconditioner DIC;
                tolerance 1e-5;
                relTol 0;
            }
            p_rgh
            {
                solver PCG;
                preconditioner DIC;
                tolerance 1e-10;
                relTol 0.05;
            }
            p_rghFinal { $p_rgh; relTol 0; }
            "U.*"
            {
                solver smoothSolver;
                smoother symGaussSeidel;
                tolerance 1e-6;
                relTol 0;
                minIter 1;
            }
        }
        PIMPLE
        {
            momentumPredictor yes;
            nOuterCorrectors 8;
            nCorrectors 3;
            nNonOrthogonalCorrectors 0;
            pRefPoint (0.199375 0.099375 0.005);
            pRefValue 0;
        }
        relaxationFactors { equations { ".*" 1; } }
    """)
    write("system/decomposeParDict", """
        numberOfSubdomains 2;
        method simple;
        simpleCoeffs { n (2 1 1); delta 0.001; }
    """)

    write_json(out / "declared-settings.json", dict(case_mode=MODE, initial_tilt=INITIAL_TILT,
        mean_water_height_m=.05,cells=[NX,NY],cell_count=NX*NY,dx_m=DX,dt_s=DT,
        steps=LAST_STEP,frame_interval=FRAME_INTERVAL,duration_s=LAST_STEP*DT,
        outer_correctors=8,momentum_predictor=True,alpha_linear_tolerance=1e-11,
        pressure_linear_tolerance=1e-10,pressure_reference_point_m=[.199375,.099375,.005]))


if __name__ == "__main__":
    main(Path(sys.argv[1]))

