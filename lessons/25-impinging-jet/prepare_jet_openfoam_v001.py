"""Prepare the matched OpenFOAM Foundation 14 impinging-jet case."""
from pathlib import Path
import shutil
import sys

from jet_common_v001 import DEPTH, H, NX, NY, NSLOT, MU, RHO, V_SLOT, W, X

END = 1250
INTERVAL = 125


def main(out):
    if shutil.disk_usage(out.anchor).free < 10 * 2**30:
        raise RuntimeError("Need 10 GiB free on the case drive")
    out.mkdir(parents=True)
    x_breaks = [X[0], -W / 2, W / 2, X[-1]]
    xyz = [(x, y, z) for z in (0., DEPTH) for y in (0., H) for x in x_breaks]
    vertices = "\n".join(f"    ({x:.17g} {y:.17g} {z:.17g})" for x, y, z in xyz)
    side = (NX - NSLOT) // 2

    def write(rel, content, kind="dictionary"):
        path = out / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        header = f"FoamFile {{ version 2.0; format ascii; class {kind}; object {path.name}; }}\n"
        path.write_text(header + content.strip() + "\n", encoding="ascii", newline="\n")

    write("system/blockMeshDict", f"""
        scale 1;
        vertices
        (
{vertices}
        );
        blocks
        (
            hex (0 1 5 4 8 9 13 12) ({side} {NY} 1) simpleGrading (1 1 1)
            hex (1 2 6 5 9 10 14 13) ({NSLOT} {NY} 1) simpleGrading (1 1 1)
            hex (2 3 7 6 10 11 15 14) ({side} {NY} 1) simpleGrading (1 1 1)
        );
        edges ();
        boundary
        (
            left {{ type patch; faces ((4 0 8 12)); }}
            right {{ type patch; faces ((3 7 15 11)); }}
            inlet {{ type patch; faces ((6 5 13 14)); }}
            plate {{ type wall; faces ((0 1 9 8) (1 2 10 9) (2 3 11 10)); }}
            ceiling {{ type wall; faces ((5 4 12 13) (7 6 14 15)); }}
            frontAndBack
            {{
                type empty;
                faces ((0 4 5 1) (1 5 6 2) (2 6 7 3) (8 9 13 12) (9 10 14 13) (10 11 15 14));
            }}
        );
        mergePatchPairs ();
    """)
    profile = "\n".join(f"        (0 {v:.17g} 0)" for v in V_SLOT)
    write("0/U", f"""
        dimensions [0 1 -1 0 0 0 0];
        internalField uniform (0 0 0);
        boundaryField
        {{
            inlet
            {{
                type fixedValue;
                value nonuniform List<vector>
                {NSLOT}
                (
{profile}
                );
            }}
            left {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
            right {{ type pressureInletOutletVelocity; value uniform (0 0 0); }}
            plate {{ type noSlip; }}
            ceiling {{ type noSlip; }}
            frontAndBack {{ type empty; }}
        }}
    """, "volVectorField")
    write("0/p", """
        dimensions [0 2 -2 0 0 0 0];
        internalField uniform 0;
        boundaryField
        {
            inlet { type zeroGradient; }
            left { type fixedValue; value uniform 0; }
            right { type fixedValue; value uniform 0; }
            plate { type zeroGradient; }
            ceiling { type zeroGradient; }
            frontAndBack { type empty; }
        }
    """, "volScalarField")
    write("constant/physicalProperties", f"viscosityModel constant;\nnu {MU / RHO:.17g};")
    write("constant/momentumTransport", "simulationType laminar;")
    write("system/controlDict", f"""
        application foamRun;
        solver incompressibleFluid;
        startFrom startTime;
        startTime 0;
        stopAt endTime;
        endTime {END};
        deltaT 1;
        writeControl timeStep;
        writeInterval {INTERVAL};
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
        relaxationFactors { equations { U 0.9; } }
    """)
    write("system/decomposeParDict", """
        numberOfSubdomains 2;
        method simple;
        simpleCoeffs { n (2 1 1); delta 0.001; }
    """)


if __name__ == "__main__":
    main(Path(sys.argv[1]))

