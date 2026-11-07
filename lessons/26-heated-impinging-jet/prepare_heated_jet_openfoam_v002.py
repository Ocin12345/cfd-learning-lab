"""Clone the verified OpenFOAM jet fields and add transient heat transport."""
import hashlib
import json
from pathlib import Path
import shutil
import sys

from heated_jet_common_v001 import ALPHA, COLD, DT, HOT, LAST, STRIDE


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main(source, out):
    assert shutil.disk_usage(out.anchor).free >= 10 * 2**30
    assert (source / "run-status.txt").read_text().strip() == "complete"
    out.mkdir(parents=True)
    hashes = {}
    for rel in ("constant/polyMesh/points", "constant/polyMesh/faces", "constant/polyMesh/owner",
                "constant/polyMesh/neighbour", "constant/polyMesh/boundary",
                "constant/physicalProperties", "constant/momentumTransport",
                "system/decomposeParDict", "1250/U", "1250/p", "1250/phi", "1250/C"):
        p = source / rel
        target = out / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        before = digest(p)
        shutil.copyfile(p, target)
        assert digest(target) == digest(p) == before
        hashes[rel] = before
    (out / "0").mkdir()
    for name in ("U", "p", "phi"):
        shutil.copyfile(out / "1250" / name, out / "0" / name)
        assert digest(out / "0" / name) == hashes[f"1250/{name}"]

    def write(rel, body, kind="dictionary"):
        p = out / rel
        header = f"FoamFile {{ version 2.0; format ascii; class {kind}; object {p.name}; }}\n"
        with p.open("x", encoding="ascii", newline="\n") as f:
            f.write(header + body.strip() + "\n")

    write("0/T", f"""
        dimensions [0 0 0 1 0 0 0];
        internalField uniform {COLD:g};
        boundaryField
        {{
            inlet {{ type fixedValue; value uniform {COLD:g}; }}
            left {{ type inletOutlet; inletValue uniform {COLD:g}; value uniform {COLD:g}; }}
            right {{ type inletOutlet; inletValue uniform {COLD:g}; value uniform {COLD:g}; }}
            plate {{ type fixedValue; value uniform {HOT:g}; }}
            ceiling {{ type zeroGradient; }}
            frontAndBack {{ type empty; }}
        }}
    """, "volScalarField")
    write("system/controlDict", f"""
        solver functions;
        subSolver incompressibleFluid;
        subSolverTime 1250;
        startFrom startTime;
        startTime 0;
        stopAt endTime;
        endTime {LAST * DT};
        deltaT {DT};
        adjustTimeStep false;
        writeControl timeStep;
        writeInterval {STRIDE};
        purgeWrite 0;
        writeFormat ascii;
        writePrecision 14;
        writeCompression off;
        timeFormat general;
        timePrecision 10;
        runTimeModifiable false;
    """)
    write("system/fvSchemes", """
        ddtSchemes { default Euler; }
        gradSchemes { default Gauss linear; }
        divSchemes
        {
            default none;
            div(phi,T) Gauss upwind;
            div(phi,U) bounded Gauss linearUpwind grad(U);
            div((nuEff*dev2(T(grad(U))))) Gauss linear;
        }
        laplacianSchemes { default Gauss linear orthogonal; }
        interpolationSchemes { default linear; }
        snGradSchemes { default orthogonal; }
    """)
    write("system/fvSolution", """
        solvers
        {
            T { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0; nCorrectors 0; }
            TFinal { $T; }
            p { solver GAMG; tolerance 1e-10; relTol 0.01; smoother GaussSeidel; }
            U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01; }
        }
        PIMPLE { nOuterCorrectors 1; nCorrectors 0; nNonOrthogonalCorrectors 0; }
        relaxationFactors { equations { T 1; } }
    """)
    write("system/functions", f"""
        transportTemperature
        {{
            type scalarTransport;
            libs ("libsolverFunctionObjects.so");
            field T;
            phi phi;
            diffusivity constant;
            D {ALPHA:.17g};
            executeControl timeStep;
            executeInterval 1;
            writeControl timeStep;
            writeInterval {STRIDE};
        }}
        previousTemperature
        {{
            type writeObjects;
            libs ("libutilityFunctionObjects.so");
            objects (T);
            writeOption anyWrite;
            writeControl timeStep;
            writeInterval {LAST - 1};
        }}
        frozenRuntimeFields
        {{
            type writeObjects;
            libs ("libutilityFunctionObjects.so");
            objects (U p phi);
            writeOption anyWrite;
            writeControl timeStep;
            writeInterval {LAST};
        }}
    """)
    (out / "frozen-input-sha256.json").write_text(json.dumps(hashes, indent=2) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))

