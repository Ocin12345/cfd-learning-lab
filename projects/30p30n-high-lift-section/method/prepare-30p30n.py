"""Extrude the source-locked APC-IV high-lift quadrilateral grid by one cell."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

GRID_SHA="dc3e70dbd11fa3e38a604705b2aa39779f99ebf3da081fd5932a75b4dcedca09"
GRID_URL="https://raw.githubusercontent.com/su2code/VandV/86a7f94ea29949f4d5e7f034d1d04ac3ef04a844/rans/30p30n/2D_L1_coarse_r1.su2.gz"
CHORD=.4572
SPAN=.01
AOA=math.radians(5.5)
SPEED=.17*math.sqrt(1.4*287.058*295.56)
NU=SPEED*CHORD/1.71e6
PARTS=("slat","mainElement","flap")
FAR="farfield"
EMPTY="frontBack"


def source_mesh(path):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==GRID_SHA
    with gzip.open(path,"rt",encoding="ascii") as f:
        def key(name):
            line=f.readline()
            assert line.split("=")[0].strip()==name,line
            return line.split("=")[1].strip()

        assert int(key("NDIME"))==2
        n=int(key("NPOIN"))
        points=[]
        for i in range(n):
            x,y,label=f.readline().split()
            assert int(label)==i
            points.append((float(x),float(y)))
        assert len(points)==64540 and len(set(points))==len(points)
        cells=[]
        for i in range(int(key("NELEM"))):
            kind,*rest=map(int,f.readline().split())
            assert kind==9 and len(rest)==5 and rest[-1]==i
            ids=tuple(rest[:-1])
            assert len(set(ids))==4 and min(ids)>=0 and max(ids)<n
            ring=[points[j] for j in ids]
            a=sum(x*yy-y*xx for (x,y),(xx,yy) in zip(ring,ring[1:]+ring[:1]))
            assert abs(a)>0
            cells.append(ids if a>0 else ids[::-1])
        assert len(cells)==63957 and int(key("NMARK"))==2
        boundaries={}
        for name,count in (("wall",610),(FAR,560)):
            assert key("MARKER_TAG")==name and int(key("MARKER_ELEMS"))==count
            edges=[]
            for _ in range(count):
                type_id,a,b=map(int,f.readline().split())
                assert type_id==3 and 0<=a<n and 0<=b<n and a!=b
                edges.append((a,b))
            boundaries[name]=edges
        assert not f.read().strip()
    return points,cells,boundaries


def separate_walls(points, edges):
    adjacency={}
    for a,b in edges:
        adjacency.setdefault(a,set()).add(b)
        adjacency.setdefault(b,set()).add(a)
    assert all(len(neighbours)==2 for neighbours in adjacency.values())
    remaining=set(adjacency)
    components=[]
    while remaining:
        todo=[min(remaining)]
        found=set()
        while todo:
            node=todo.pop()
            if node not in found:
                found.add(node)
                todo.extend(adjacency[node]-found)
        remaining-=found
        components.append(found)
    assert len(components)==3
    components.sort(key=lambda ids:sum(points[i][0] for i in ids)/len(ids))
    return {name:[(a,b) for a,b in edges if a in ids and b in ids]
            for name,ids in zip(PARTS,components)}


def setup(path, source):
    points,cells,markers=source_mesh(source)
    wall_groups=separate_walls(points,markers["wall"])
    markers={**wall_groups,FAR:markers[FAR]}
    patch_map={}
    for name,edges in markers.items():
        for a,b in edges:
            key=min(a,b),max(a,b)
            assert key not in patch_map
            patch_map[key]=name
    n=len(points)
    halves={}
    inside=[]
    side_planes=[]
    for i,ids in enumerate(cells):
        side_planes.extend(((tuple(reversed(ids)),i),(tuple(j+n for j in ids),i)))
        for a,b in zip(ids,ids[1:]+ids[:1]):
            key=min(a,b),max(a,b)
            if key in halves:
                aa,bb,owner=halves.pop(key)
                assert (aa,bb)==(b,a) and owner<i and key not in patch_map
                inside.append(((aa,bb,bb+n,aa+n),owner,i))
            else:
                halves[key]=(a,b,i)
    assert set(halves)==set(patch_map)
    boundary={name:[] for name in markers}
    for key,(a,b,i) in halves.items():
        boundary[patch_map[key]].append(((a,b,b+n,a+n),i))
    boundary[EMPTY]=side_planes

    path.mkdir()
    def put(name,body,cls="dictionary",loc=""):
        file=path/name
        file.parent.mkdir(parents=True,exist_ok=True)
        with file.open("x",encoding="ascii",newline="\n") as f:
            f.write(f'FoamFile {{ version 2.0; format ascii; class {cls}; location "{loc}"; object {file.name}; }}\n')
            f.write(body.strip()+"\n")

    def foam_list(values):
        return str(len(values))+"\n(\n"+"\n".join(values)+"\n)"
    mesh="constant/polyMesh"
    put(mesh+"/points",foam_list([f"({x:.17g} {y:.17g} {z:g})"
                                 for z in (-SPAN/2,SPAN/2) for x,y in points]),"vectorField",mesh)
    inside.sort(key=lambda record:(record[1],record[2]))
    faces=[face for face,*_ in inside]
    owners=[i for _,i,_ in inside]
    start=len(faces)
    patch_list=[]
    for name,items in boundary.items():
        kind="wall" if name in PARTS else "empty" if name==EMPTY else "patch"
        patch_list.append(f"{name}\n{{\ntype {kind};\nnFaces {len(items)};\nstartFace {start};\n}}")
        start+=len(items)
        faces.extend(face for face,_ in items)
        owners.extend(i for _,i in items)
    put(mesh+"/faces",foam_list([str(len(ids))+"("+" ".join(map(str,ids))+")" for ids in faces]),"faceList",mesh)
    put(mesh+"/owner",foam_list(list(map(str,owners))),"labelList",mesh)
    put(mesh+"/neighbour",foam_list([str(j) for _,_,j in inside]),"labelList",mesh)
    put(mesh+"/boundary",foam_list(patch_list),"polyBoundaryMesh",mesh)
    with (path/"mesh-source.json").open("x",encoding="ascii") as f:
        json.dump({"grid_sha256":GRID_SHA,"grid_url":GRID_URL,"source_points":n,"cells":len(cells),
                   "source_wall_edges":{name:len(items) for name,items in wall_groups.items()},
                   "extrusion_span_m":SPAN,"stowed_chord_m":CHORD,"velocity_m_s":SPEED,
                   "kinematic_viscosity_m2_s":NU},f,indent=2)
        f.write("\n")

    put("constant/physicalProperties",f"viscosityModel constant; nu {NU:.17g};")
    put("constant/momentumTransport","simulationType RAS;\nRAS { model SpalartAllmaras; turbulence on; printCoeffs on; }")
    put("system/controlDict","""
solver incompressibleFluid;
startFrom startTime; startTime 0; stopAt endTime; endTime 1600; deltaT 1;
writeControl timeStep; writeInterval 100; purgeWrite 0;
writeFormat ascii; writePrecision 14; writeCompression off; runTimeModifiable no;
""")
    put("system/fvSolution","""
solvers {
    p { solver GAMG; tolerance 1e-8; relTol 0.02; smoother GaussSeidel; }
    "(U|nuTilda)" { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-10; relTol 0.1; }
}
SIMPLE { nNonOrthogonalCorrectors 1; }
relaxationFactors { fields { p 0.3; } equations { U 0.7; nuTilda 0.7; } }
""")
    put("system/functions","\n".join(f"""
{name}
{{
    type forces; libs ("libforces.so");
    patches ({patches});
    rho rhoInf; rhoInf 1.225; CofR (0 0 0); pitchAxis (0 0 1);
    writeControl timeStep; writeInterval 1;
}}
""" for name,patches in (("allElements"," ".join(PARTS)),*((name+"Force",name) for name in PARTS))))
    velocity=f"({SPEED*math.cos(AOA):.16g} {SPEED*math.sin(AOA):.16g} 0)"
    for key,cls,dim,value,far,wall in (
            ("U","volVectorField","[velocity]",velocity,"freestreamVelocity","noSlip;"),
            ("p","volScalarField","[kinematicPressure]","0","freestreamPressure","zeroGradient;"),
            ("nuTilda","volScalarField","[turbulentViscosity]",f"{3*NU:.17g}","freestream","fixedValue; value uniform 0;"),
            ("nut","volScalarField","[turbulentViscosity]","0","calculated","nutLowReWallFunction; value uniform 0;")):
        field=f"dimensions {dim};\ninternalField uniform {value};\nboundaryField\n{{\n"
        field+=f"{FAR} {{ type {far}; "
        if key!="nut":
            field+=f"freestreamValue uniform {value}; "
        field+=f"value uniform {value}; }}\n"
        field+="\n".join(f"{name} {{ type {wall} }}" for name in PARTS)
        field+=f"\n{EMPTY} {{ type empty; }}\n}}"
        put("0/"+key,field,cls,"0")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path)
    parser.add_argument("source",type=Path)
    args=parser.parse_args()
    setup(args.case,args.source)
