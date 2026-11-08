"""Prepare a matching transient shockFluid case for OpenFOAM Foundation 14."""
from math import tan
from pathlib import Path
import shutil
import sys

from wedge_common_fine_v001 import CORNER, CP, DEPTH, H, L, MW, NX, NY, P, T, THETA, U

FRAME_DT = 0.00025
END_TIME = 0.008
MAX_CO = 0.3


def main(out):
    if shutil.disk_usage(out.anchor).free < 10 * 2**30:
        raise RuntimeError("Need 10 GiB free")
    out.mkdir(parents=True)
    ramp_height = (L - CORNER) * tan(THETA)
    leading_cells = round(NX * CORNER / L)
    xy = [(0, 0), (CORNER, 0), (CORNER, H), (0, H), (L, ramp_height), (L, H)]
    vertices = "\n".join(f"    ({x:.17g} {y:.17g} {z:.17g})" for z in (0, DEPTH) for x, y in xy)

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
            hex (0 1 2 3 6 7 8 9) ({leading_cells} {NY} 1) simpleGrading (1 1 1)
            hex (1 4 5 2 7 10 11 8) ({NX - leading_cells} {NY} 1) simpleGrading (1 1 1)
        );
        edges ();
        boundary
        (
            inlet {{ type patch; faces ((3 0 6 9)); }}
            outlet {{ type patch; faces ((4 5 11 10)); }}
            top {{ type wall; faces ((2 3 9 8) (5 2 8 11)); }}
            ramp {{ type wall; faces ((0 1 7 6) (1 4 10 7)); }}
            frontAndBack
            {{
                type empty;
                faces ((0 3 2 1) (1 2 5 4) (6 7 8 9) (7 10 11 8));
            }}
        );
        mergePatchPairs ();
    """)
    for name, units, value, kind in (("p", "pressure", f"{P:.16g}", "volScalarField"),
                                     ("T", "temperature", f"{T:.16g}", "volScalarField"),
                                     ("U", "velocity", f"({U:.16g} 0 0)", "volVectorField")):
        wall_type = "slip" if name == "U" else "zeroGradient"
        write(f"0/{name}", f"""
            dimensions [{units}];
            internalField uniform {value};
            boundaryField
            {{
                inlet {{ type fixedValue; value uniform {value}; }}
                outlet {{ type zeroGradient; }}
                ramp {{ type {wall_type}; }}
                top {{ type {wall_type}; }}
                frontAndBack {{ type empty; }}
            }}
        """, kind)
    write("constant/physicalProperties", f"""
        thermoType
        {{
            type hePsiThermo;
            mixture pureMixture;
            transport const;
            thermo hConst;
            equationOfState perfectGas;
            specie specie;
            energy sensibleInternalEnergy;
        }}
        mixture
        {{
            specie {{ molWeight {MW:.16g}; }}
            thermodynamics {{ Cp {CP:.16g}; hf 0; }}
            transport {{ mu 0; Pr 1; }}
        }}
    """)
    write("constant/momentumTransport", "simulationType laminar;")
    write("system/controlDict", f"""
        application foamRun;
        solver shockFluid;
        startFrom startTime;
        startTime 0;
        stopAt endTime;
        endTime {END_TIME};
        deltaT 5e-7;
        writeControl adjustableRunTime;
        writeInterval {FRAME_DT};
        purgeWrite 0;
        writeFormat ascii;
        writePrecision 12;
        writeCompression off;
        timeFormat general;
        timePrecision 9;
        runTimeModifiable false;
        adjustTimeStep yes;
        maxCo {MAX_CO};
        maxDeltaT 2e-6;
    """)
    write("system/fvSchemes", """
        fluxScheme Kurganov;
        ddtSchemes { default Euler; }
        gradSchemes { default Gauss linear; }
        divSchemes { default none; }
        laplacianSchemes { default Gauss linear corrected; }
        interpolationSchemes
        {
            default linear;
            reconstruct(rho) vanLeer;
            reconstruct(U) vanLeerV;
            reconstruct(T) vanLeer;
        }
        snGradSchemes { default corrected; }
    """)
    write("system/fvSolution", """
        solvers
        {
            "(rho|rhoU|rhoE).*" { solver diagonal; }
            "U.*"
            {
                solver smoothSolver;
                smoother GaussSeidel;
                nSweeps 2;
                tolerance 1e-9;
                relTol 0.01;
            }
            "e.*" { $U; tolerance 1e-10; relTol 0; }
        }
        PIMPLE {}
    """)
    write("system/decomposeParDict", """
        numberOfSubdomains 2;
        method simple;
        simpleCoeffs { n (2 1 1); delta 0.001; }
    """)


if __name__ == "__main__":
    main(Path(sys.argv[1]))

