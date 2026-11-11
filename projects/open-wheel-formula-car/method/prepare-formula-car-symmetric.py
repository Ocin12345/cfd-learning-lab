"""Generate an idealized sports prototype and the corresponding OpenFOAM inputs."""
import argparse
import math
import struct
from collections import Counter
from pathlib import Path

SPEED = 30.
RADIUS = .32
DENSITY = 1.225
DIFFUSER_START = 3.20
LENGTH = 4.55
RIDE_HEIGHT = .18


def chassis():
    # x, half width, roof elevation, flat-floor elevation.
    stations = [(.84,.205,.56,.29),(1.02,.24,.60,.275),(1.38,.31,.63,.23),
                (1.56,.48,.70,.18),(1.82,.59,.98,.18),(2.12,.64,1.11,.18),
                (2.55,.63,.97,.18),(2.87,.62,.72,.18),(3.2,.61,.62,.18),
                (3.9,.53,.59,.18),(4.55,.42,.60,.18)]
    loops = []
    for x,w,roof,z0 in stations:
        z = z0 + max(0,x-DIFFUSER_START)*math.tan(math.radians(8))
        section = [(-w+.07,z),(w-.07,z),(w,z+.065),(w,roof-.08),
                   (w-.12,roof),(-w+.12,roof),(-w,roof-.08),(-w,z+.065)]
        loops.append([(x,y,h) for y,h in section])
    return closed_loft(loops)


def closed_loft(loops):
    n = len(loops[0])
    output = []
    for a,b in zip(loops,loops[1:]):
        for i in range(n):
            j = (i+1)%n
            output.extend([(a[i],a[j],b[j]),(a[i],b[j],b[i])])
    for i in range(1,n-1):
        output.extend([(loops[0][0],loops[0][i+1],loops[0][i]),
                       (loops[-1][0],loops[-1][i],loops[-1][i+1])])
    return output


def axle(x):
    output = []
    for low,high in ((-1.04,-.79),(.79,1.04)):
        sections=[]
        for y,scale in ((low,.94),(low+.018,1.),(high-.018,1.),(high,.94)):
            # z and x form a right-handed cross-section with normal +y.
            sections.append([(x+RADIUS*scale*math.sin(2*math.pi*i/48),y,
                              .38+RADIUS*scale*math.cos(2*math.pi*i/48)) for i in range(48)])
        output.extend(closed_loft(sections))
    return output


def wing(x0, z0, chord, half_span, angle, thickness):
    # Each profile is generated around an inverted, piecewise-parabolic
    # camber line. Cosine spacing resolves the nose and the closed tail.
    stations=[.5*(1-math.cos(math.pi*i/30)) for i in range(31)]
    thickness_values=[5*thickness*(.2969*math.sqrt(x)-.126*x-.3516*x*x
                     +.2843*x**3-.1036*x**4) for x in stations]
    outline=[]
    for i in list(range(31))+list(range(29,0,-1)):
        outline.append((stations[i],thickness_values[i]))
    # z/x order gives an outward +y normal for the right end cap.
    xz=[]
    for i,(x,delta) in enumerate(outline):
        camber=-.045*(2*.4*x-x*x)/.4**2 if x<.4 else -.045*(.2+.8*x-x*x)/.6**2
        d=delta if i<31 else -delta
        z=camber+d
        theta=math.radians(angle)
        xz.append((x0+chord*(x*math.cos(theta)-z*math.sin(theta)),
                   z0+chord*(x*math.sin(theta)+z*math.cos(theta))))
    # Top surface goes towards the trailing edge, giving a +y cap.
    return closed_loft([[(x,y,z) for x,z in xz] for y in (-half_span,half_span)])


def write_stl(path,triangles):
    edges = Counter((a,b) for tri in triangles for a,b in zip(tri,tri[1:]+tri[:1]))
    assert all(n==1 and edges[b,a]==1 for (a,b),n in edges.items())
    with path.open("xb") as f:
        f.write(b"Original generic open-wheel formula car".ljust(80,b" "))
        f.write(struct.pack("<I",len(triangles)))
        for a,b,c in triangles:
            u,v = [y-x for x,y in zip(a,b)],[y-x for x,y in zip(a,c)]
            normal = [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
            length=math.hypot(*normal)
            assert length>0
            f.write(struct.pack("<12fH",*(n/length for n in normal),*a,*b,*c,0))


def setup(root, angle):
    root.mkdir()
    geometry = root/"constant/geometry"
    geometry.mkdir(parents=True)
    parts=(("chassis",chassis()),("frontWheels",axle(1.23)),("rearWheels",axle(3.62)),
           ("frontMain",wing(.035,.238,.44,.89,3,.14)),
           ("frontFlap",wing(.51,.303,.26,.85,angle,.14)),
           ("rearWing",wing(3.78,1.03,.53,.68,11,.13)))
    for name, triangles in parts:
        write_stl(geometry/(name+".stl"),triangles)

    def write(name,body,kind="dictionary"):
        path=root/name
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open("x",encoding="ascii",newline="\n") as f:
            f.write(f"FoamFile {{ version 2.0; format ascii; class {kind}; object {path.name}; }}\n")
            f.write(body.strip()+"\n")

    write("system/blockMeshDict","""
vertices ((-5 0 0) (14 0 0) (14 4 0) (-5 4 0) (-5 0 4.5) (14 0 4.5) (14 4 4.5) (-5 4 4.5));
blocks (hex (0 1 2 3 4 5 6 7) (76 16 18) simpleGrading (1 1 1));
boundary (
    inlet { type patch; faces ((0 4 7 3)); }
    outlet { type patch; faces ((1 2 6 5)); }
    road { type wall; faces ((0 3 2 1)); }
    centrePlane { type symmetryPlane; faces ((0 1 5 4)); }
    farfield { type symmetry; faces ((3 7 6 2) (4 5 6 7)); }
);
""")
    write("system/snappyHexMeshDict","""
castellatedMesh true;
snap true;
addLayers false;
geometry {
    chassis { type triSurface; file "chassis.stl"; }
    frontWheels { type triSurface; file "frontWheels.stl"; }
    rearWheels { type triSurface; file "rearWheels.stl"; }
    frontMain { type triSurface; file "frontMain.stl"; }
    frontFlap { type triSurface; file "frontFlap.stl"; }
    rearWing { type triSurface; file "rearWing.stl"; }
    floorZone { type box; min (-.25 -1.18 -.01); max (4.85 1.18 .75); }
    wakeZone { type box; min (3.6 -1.5 -.01); max (8.1 1.5 1.65); }
}
castellatedMeshControls {
    maxLocalCells 240000; maxGlobalCells 240000;
    minRefinementCells 0; nCellsBetweenLevels 2; maxLoadUnbalance 0.1;
    features ();
    refinementSurfaces {
        chassis { level (3 3); patchInfo { type wall; } }
        frontWheels { level (3 3); patchInfo { type wall; } }
        rearWheels { level (3 3); patchInfo { type wall; } }
        frontMain { level (4 4); patchInfo { type wall; } }
        frontFlap { level (4 4); patchInfo { type wall; } }
        rearWing { level (4 4); patchInfo { type wall; } }
    }
    resolveFeatureAngle 35;
    refinementRegions {
        floorZone { mode inside; level 2; }
        wakeZone { mode inside; level 1; }
    }
    insidePoint (-3.17 .63 1.13);
    allowFreeStandingZoneFaces true;
}
snapControls {
    nSmoothPatch 5; tolerance 2; nSolveIter 40; nRelaxIter 6;
    nFeatureSnapIter 12; implicitFeatureSnap true; explicitFeatureSnap false;
}
addLayersControls {}
meshQualityControls {
    #includeEtc "caseDicts/mesh/generation/meshQualityDict"
    minFaceWeight .03;
    maxBoundarySkewness 3.9;
}
mergeTolerance 1e-6;
""")
    write("system/controlDict","""
solver incompressibleFluid;
startFrom startTime; startTime 0; stopAt endTime; endTime 600; deltaT 1;
writeControl timeStep; writeInterval 50; purgeWrite 0;
writeFormat ascii; writePrecision 14; writeCompression off; runTimeModifiable no;
""")
    write("system/functions","""
vehicleForce
{
    type forces; libs ("libforces.so");
    patches (chassis frontWheels rearWheels frontMain frontFlap rearWing);
    rho rhoInf; rhoInf 1.225; CofR (2.17 0 .18); pitchAxis (0 1 0);
    writeControl timeStep; writeInterval 1;
}
chassisForce
{
    type forces; libs ("libforces.so");
    patches (chassis);
    rho rhoInf; rhoInf 1.225; CofR (2.17 0 .18); pitchAxis (0 1 0);
    writeControl timeStep; writeInterval 1;
}
""")
    with (root/"system/functions").open("a",encoding="ascii",newline="\n") as f:
        for name,group in (("frontWingForce","frontMain frontFlap"),("flapForce","frontFlap"),("rearWingForce","rearWing")):
            f.write(f"""
{name}
{{
    type forces; libs ("libforces.so");
    patches ({group});
    rho rhoInf; rhoInf 1.225; CofR (2.17 0 .18); pitchAxis (0 1 0);
    writeControl timeStep; writeInterval 1;
}}
""")
    write("constant/momentumTransport","simulationType RAS;\nRAS { model kOmegaSST; turbulence on; printCoeffs on; }")
    write("constant/physicalProperties","viscosityModel constant; nu 1.5e-5;")

    walls=("chassis","frontWheels","rearWheels","frontMain","frontFlap","rearWing","road")
    omega= -SPEED/RADIUS
    vel={
        "inlet":"type fixedValue; value uniform (30 0 0);",
        "outlet":"type pressureInletOutletVelocity; value uniform (30 0 0);",
        "farfield":"type symmetry;", "centrePlane":"type symmetryPlane;",
        **{name:"type noSlip;" for name in ("chassis","frontMain","frontFlap","rearWing")},
        "road":"type fixedValue; value uniform (30 0 0);",
    }
    for name,x in (("frontWheels",1.23),("rearWheels",3.62)):
        vel[name]=f"type rotatingWallVelocity; origin ({x} 0 .38); axis (0 1 0); omega constant {omega:g}; value uniform (0 0 0);"
    patches={"inlet":"type zeroGradient;","outlet":"type fixedValue; value uniform 0;",
             "farfield":"type symmetry;", "centrePlane":"type symmetryPlane;",**{wall:"type zeroGradient;" for wall in walls}}
    values={
        "U":("volVectorField","[velocity]","(30 0 0)",vel),
        "p":("volScalarField","[kinematicPressure]","0",patches),
    }
    for name,dim,value,wall_fn in (("k","[0 2 -2 0 0 0 0]",".375","kqRWallFunction"),
                                  ("omega","[0 0 -1 0 0 0 0]","45","omegaWallFunction"),
                                  ("nut","[0 2 -1 0 0 0 0]","0","nutkWallFunction")):
        inlet=f"type fixedValue; value uniform {value};" if name!="nut" else "type calculated; value uniform 0;"
        outlet=f"type inletOutlet; inletValue uniform {value}; value uniform {value};" if name!="nut" else inlet
        conditions={"inlet":inlet,"outlet":outlet,"farfield":"type symmetry;", "centrePlane":"type symmetryPlane;",
                    **{wall:f"type {wall_fn}; value uniform {value};" for wall in walls}}
        values[name]=("volScalarField",dim,value,conditions)
    for name,(kind,dim,value,conditions) in values.items():
        boundary="\n".join(f"{patch} {{ {entry} }}" for patch,entry in conditions.items())
        write("0/"+name,f"dimensions {dim};\ninternalField uniform {value};\nboundaryField\n{{\n{boundary}\n}}",kind)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("angle",type=int,choices=(8,18))
    args=parser.parse_args()
    setup(args.root,args.angle)
