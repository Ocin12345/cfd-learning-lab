"""Generate an idealized sports prototype and the corresponding OpenFOAM inputs."""
import argparse
import math
import struct
from collections import Counter
from pathlib import Path

SPEED = 25.
RADIUS = .32
DENSITY = 1.225
DIFFUSER_START = 3.20
LENGTH = 4.4
RIDE_HEIGHT = .18


def chassis(angle):
    # x, half width, roof elevation, flat-floor elevation.
    stations = [(0,.27,.50,.25),(.33,.42,.54,.195),(.8,.62,.60,RIDE_HEIGHT),
                (1.35,.67,.66,RIDE_HEIGHT),(1.65,.67,.93,RIDE_HEIGHT),
                (2.12,.67,1.07,RIDE_HEIGHT),(2.66,.67,.95,RIDE_HEIGHT),
                (3.2,.67,.68,RIDE_HEIGHT),(3.85,.65,.62,RIDE_HEIGHT),
                (4.4,.57,.60,RIDE_HEIGHT)]
    loops = []
    for x,w,roof,z0 in stations:
        z = z0 + max(0,x-DIFFUSER_START)*math.tan(math.radians(angle))
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


def write_stl(path,triangles):
    edges = Counter((a,b) for tri in triangles for a,b in zip(tri,tri[1:]+tri[:1]))
    assert all(n==1 and edges[b,a]==1 for (a,b),n in edges.items())
    with path.open("xb") as f:
        f.write(b"Original idealized sports prototype".ljust(80,b" "))
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
    for name, triangles in (("chassis",chassis(angle)),("frontWheels",axle(.90)),("rearWheels",axle(3.44))):
        write_stl(geometry/(name+".stl"),triangles)

    def write(name,body,kind="dictionary"):
        path=root/name
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open("x",encoding="ascii",newline="\n") as f:
            f.write(f"FoamFile {{ version 2.0; format ascii; class {kind}; object {path.name}; }}\n")
            f.write(body.strip()+"\n")

    write("system/blockMeshDict","""
vertices ((-5 -4 0) (14 -4 0) (14 4 0) (-5 4 0) (-5 -4 4.5) (14 -4 4.5) (14 4 4.5) (-5 4 4.5));
blocks (hex (0 1 2 3 4 5 6 7) (76 32 18) simpleGrading (1 1 1));
boundary (
    inlet { type patch; faces ((0 4 7 3)); }
    outlet { type patch; faces ((1 2 6 5)); }
    road { type wall; faces ((0 3 2 1)); }
    farfield { type symmetry; faces ((0 1 5 4) (3 7 6 2) (4 5 6 7)); }
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
    patches (chassis frontWheels rearWheels);
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
    write("constant/momentumTransport","simulationType RAS;\nRAS { model kOmegaSST; turbulence on; printCoeffs on; }")
    write("constant/physicalProperties","viscosityModel constant; nu 1.5e-5;")

    walls=("chassis","frontWheels","rearWheels","road")
    omega= -SPEED/RADIUS
    vel={
        "inlet":"type fixedValue; value uniform (25 0 0);",
        "outlet":"type pressureInletOutletVelocity; value uniform (25 0 0);",
        "farfield":"type symmetry;",
        "chassis":"type noSlip;",
        "road":"type fixedValue; value uniform (25 0 0);",
    }
    for name,x in (("frontWheels",.90),("rearWheels",3.44)):
        vel[name]=f"type rotatingWallVelocity; origin ({x} 0 .38); axis (0 1 0); omega constant {omega:g}; value uniform (0 0 0);"
    patches={"inlet":"type zeroGradient;","outlet":"type fixedValue; value uniform 0;",
             "farfield":"type symmetry;",**{wall:"type zeroGradient;" for wall in walls}}
    values={
        "U":("volVectorField","[velocity]","(25 0 0)",vel),
        "p":("volScalarField","[kinematicPressure]","0",patches),
    }
    for name,dim,value,wall_fn in (("k","[0 2 -2 0 0 0 0]",".375","kqRWallFunction"),
                                  ("omega","[0 0 -1 0 0 0 0]","45","omegaWallFunction"),
                                  ("nut","[0 2 -1 0 0 0 0]","0","nutkWallFunction")):
        inlet=f"type fixedValue; value uniform {value};" if name!="nut" else "type calculated; value uniform 0;"
        outlet=f"type inletOutlet; inletValue uniform {value}; value uniform {value};" if name!="nut" else inlet
        conditions={"inlet":inlet,"outlet":outlet,"farfield":"type symmetry;",
                    **{wall:f"type {wall_fn}; value uniform {value};" for wall in walls}}
        values[name]=("volScalarField",dim,value,conditions)
    for name,(kind,dim,value,conditions) in values.items():
        boundary="\n".join(f"{patch} {{ {entry} }}" for patch,entry in conditions.items())
        write("0/"+name,f"dimensions {dim};\ninternalField uniform {value};\nboundaryField\n{{\n{boundary}\n}}",kind)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("angle",type=int,choices=(0,10))
    args=parser.parse_args()
    setup(args.root,args.angle)
