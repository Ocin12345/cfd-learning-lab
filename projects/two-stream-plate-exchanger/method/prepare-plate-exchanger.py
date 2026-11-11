"""Build a new matched three-region planar heat-exchanger case."""
import argparse
from pathlib import Path

LENGTH = .12
WIDTH = .01
GAP = .002
PLATE = .0005
BULK_SPEED = .015
LAYERS = ("cold", "plate", "hot")
DIM = {"U": "[velocity]", "p": "[pressure]", "T": "[temperature]"}


def setup(case, orientation, level):
    case.mkdir()
    nx, ny, nw = {"base": (120, 16, 4), "fine": (180, 24, 6)}[level]

    def write(name, body, cls="dictionary", obj=None):
        path = case/name
        path.parent.mkdir(parents=True, exist_ok=True)
        header = f"FoamFile {{ version 2.0; format ascii; class {cls}; object {obj or path.name}; }}"
        with path.open("x", encoding="ascii", newline="\n") as stream:
            stream.write(header+"\n"+body.strip()+"\n")

    points = []
    for y in (0, GAP, GAP+PLATE, 2*GAP+PLATE):
        points.extend([(0,y,0), (LENGTH,y,0), (0,y,WIDTH), (LENGTH,y,WIDTH)])
    sections = []
    ends = []
    side_faces = []
    for j, (region, nyj) in enumerate(zip(LAYERS, (ny, nw, ny))):
        a,b,c,d,e,f,g,h = range(j*4,j*4+8)
        sections.append(f"hex ({a} {b} {f} {e} {c} {d} {h} {g}) {region} ({nx} {nyj} 1) simpleGrading (1 1 1)")
        inlet, outlet = f"({a} {c} {g} {e})", f"({b} {f} {h} {d})"
        patch_type = "wall" if region == "plate" else "patch"
        ends += [f"{region}_left {{ type {patch_type}; faces ({inlet}); }}",
                 f"{region}_right {{ type {patch_type}; faces ({outlet}); }}"]
        side_faces += [f"({a} {e} {f} {b})", f"({c} {d} {h} {g})"]
    patches = ends + ["outer_lower { type wall; faces ((0 1 3 2)); }",
                      "outer_upper { type wall; faces ((12 14 15 13)); }",
                      "front_back { type empty; faces ("+" ".join(side_faces)+"); }"]
    vertices = "\n".join("(%g %g %g)" % point for point in points)
    write("system/blockMeshDict", "vertices (\n"+vertices+"\n);\nblocks (\n"+"\n".join(sections)
          +"\n);\nboundary (\n"+"\n".join(patches)+"\n);")
    write("system/controlDict", """
application foamMultiRun;
regionSolvers { cold fluid; plate solid; hot fluid; }
startFrom startTime; startTime 0; stopAt endTime; endTime 900; deltaT 1;
writeControl timeStep; writeInterval 100; purgeWrite 0;
writeFormat ascii; writePrecision 14; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false; adjustTimeStep no;
""")
    write("system/fvSolution", "PIMPLE { nOuterCorrectors 1; }")
    solid_schemes = """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear uncorrected; }
interpolationSchemes { default linear; }
snGradSchemes { default uncorrected; }
"""
    write("system/plate/fvSchemes", solid_schemes)
    write("system/plate/fvSolution", """
solvers { e { solver PCG; preconditioner DIC; tolerance 1e-11; relTol 0; maxIter 400; } eFinal { $e; } }
PIMPLE { nNonOrthogonalCorrectors 0; }
relaxationFactors { equations { e 1; } }
""")
    write("constant/plate/physicalProperties", """
thermoType constSolidThermo;
rho { type uniform; value 7900; }
Cv { type uniform; value 500; }
kappa { type uniform; value 16; }
""")
    fluid_schemes = solid_schemes.replace("divSchemes { default none; }", """
divSchemes {
    default none;
    div(phi,U) bounded Gauss linearUpwind grad(U);
    div(phi,h) bounded Gauss upwind;
    div(phi,K) Gauss upwind;
    div(phi,(p|rho)) Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
wallDist { method meshWave; }
""")
    for region in ("cold", "hot"):
        write(f"system/{region}/fvSchemes", fluid_schemes)
        write(f"system/{region}/fvSolution", """
solvers {
    "(rho|rhoFinal)" { solver diagonal; }
    p { solver GAMG; smoother GaussSeidel; tolerance 1e-11; relTol 0; maxIter 400; }
    pFinal { $p; }
    "(U|h)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-11; relTol 0; maxIter 400; }
    "(U|h)Final" { $U; }
}
PIMPLE { momentumPredictor yes; nCorrectors 1; nNonOrthogonalCorrectors 0; }
relaxationFactors { fields { p 0.5; } equations { U 0.7; h 1; } }
""")
        write(f"constant/{region}/physicalProperties", """
thermoType { type heRhoThermo; mixture pureMixture; transport const; thermo hConst; equationOfState rhoConst; specie specie; energy sensibleEnthalpy; }
mixture { specie { molWeight 18.01528; } equationOfState { rho 1000; } thermodynamics { Cp 4180; hf 0; } transport { mu 0.001; Pr 6.966666666666667; } }
""")
        write(f"constant/{region}/momentumTransport", "simulationType laminar;")
        write(f"constant/{region}/thermophysicalTransport", "laminar { model Fourier; }")

    def field(region, name, initial, boundaries):
        kind = "volVectorField" if name == "U" else "volScalarField"
        body = f"dimensions {DIM[name]};\ninternalField uniform {initial};\nboundaryField\n{{\n"
        body += "\n".join(f"{key} {{ {value} }}" for key, value in boundaries.items())
        write(f"0/{region}/{name}", body+"\n}", kind)

    for region, wall, temp in (("cold", "outer_lower", 300), ("hot", "outer_upper", 340)):
        speed = -BULK_SPEED if region=="cold" and orientation=="counter" else BULK_SPEED
        left, right = f"{region}_left", f"{region}_right"
        inlet, outlet = (left,right) if speed > 0 else (right,left)
        interface = f"{region}_to_plate"
        # Both port face orders follow increasing local wall-normal coordinate.
        # A parabolic profile is symmetric under the alternate patch orientation.
        values = []
        for i in range(ny):
            lo, hi = i/ny, (i+1)/ny
            face_average = speed*ny*(3*(hi**2-lo**2)-2*(hi**3-lo**3))
            values.append(f"({face_average:.14g} 0 0)")
        vel = f"({speed:g} 0 0)"
        profile = "value nonuniform List<vector> "+str(ny)+" ( "+" ".join(values)+" );"
        field(region,"U",vel,{inlet:"type fixedValue; "+profile,
                              outlet:"type pressureInletOutletVelocity; value uniform "+vel+";",
                              wall:"type noSlip;", interface:"type noSlip;", "front_back":"type empty;"})
        field(region,"p","0",{inlet:"type zeroGradient;", outlet:"type fixedValue; value uniform 0;",
                               wall:"type zeroGradient;", interface:"type zeroGradient;", "front_back":"type empty;"})
        field(region,"T",str(temp),{inlet:f"type fixedValue; value uniform {temp};",
                                   outlet:f"type inletOutlet; inletValue uniform {temp}; value uniform {temp};",
                                   wall:"type zeroGradient;", interface:f"type coupledTemperature; value uniform {temp};",
                                   "front_back":"type empty;"})
    field("plate","T","320",{"plate_left":"type zeroGradient;", "plate_right":"type zeroGradient;",
                             "plate_to_cold":"type coupledTemperature; value uniform 320;",
                             "plate_to_hot":"type coupledTemperature; value uniform 320;",
                             "front_back":"type empty;"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("orientation", choices=("counter", "parallel"))
    parser.add_argument("level", choices=("base", "fine"))
    args = parser.parse_args()
    setup(args.case, args.orientation, args.level)
