"""Small DFG-inspired Re=100 cylinder/channel mesh; no solver is launched."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import shutil

p = argparse.ArgumentParser()
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
out = a.out.resolve()
if out.exists() or out.drive.upper() != 'E:':
    raise ValueError('Choose a fresh external E: mesh directory')
if shutil.disk_usage(out.anchor).free < 10 * 1024**3:
    raise RuntimeError('Require 10 GiB free disk space')
out.mkdir(parents=True)
import gmsh
gmsh.initialize()
try:
    gmsh.option.setNumber('General.NumThreads', 2)
    gmsh.option.setNumber('Mesh.Algorithm', 6)
    gmsh.option.setNumber('Mesh.ElementOrder', 1)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 1)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.RandomSeed', 1)
    gmsh.model.add('cylinder-channel-Re100')
    geo = gmsh.model.geo
    pts = [geo.addPoint(x, y, 0, .012) for x, y in [(0,0), (2.2,0), (2.2,.41), (0,.41)]]
    outer = [geo.addLine(pts[i], pts[(i+1)%4]) for i in range(4)]
    centre = geo.addPoint(.2, .2, 0)
    circle = [geo.addPoint(.2+.05*math.cos(t), .2+.05*math.sin(t), 0, .0025)
              for t in [0, math.pi/2, math.pi, 3*math.pi/2]]
    arcs = [geo.addCircleArc(circle[i], centre, circle[(i+1)%4]) for i in range(4)]
    surf = geo.addPlaneSurface([geo.addCurveLoop(outer), geo.addCurveLoop(arcs)])
    geo.synchronize()
    for arc in arcs:
        gmsh.model.mesh.setTransfiniteCurve(arc, 33)
    gmsh.model.mesh.setTransfiniteCurve(outer[3], 42)
    gmsh.model.mesh.setTransfiniteCurve(outer[1], 42)
    for tag, name, curves in [(4,'cylinder',arcs), (5,'inlet',[outer[3]]),
                              (6,'outlet',[outer[1]]), (7,'walls',[outer[0],outer[2]])]:
        gmsh.model.addPhysicalGroup(1, curves, tag, name)
    gmsh.model.addPhysicalGroup(2, [surf], 2, 'fluid')
    dist = gmsh.model.mesh.field.add('Distance')
    gmsh.model.mesh.field.setNumbers(dist, 'CurvesList', arcs)
    gmsh.model.mesh.field.setNumber(dist, 'Sampling', 200)
    near = gmsh.model.mesh.field.add('MathEval')
    gmsh.model.mesh.field.setString(near, 'F', f'0.0025 + 0.16*F{dist}')
    wake = gmsh.model.mesh.field.add('Box')
    for key, val in dict(VIn=.009, VOut=.025, XMin=.22, XMax=2.2, YMin=.08, YMax=.33, Thickness=.04).items():
        gmsh.model.mesh.field.setNumber(wake,key,val)
    walls = gmsh.model.mesh.field.add('Distance')
    gmsh.model.mesh.field.setNumbers(walls, 'CurvesList', [outer[0], outer[2]])
    gmsh.model.mesh.field.setNumber(walls, 'Sampling', 250)
    wallsize = gmsh.model.mesh.field.add('MathEval')
    gmsh.model.mesh.field.setString(wallsize, 'F', f'0.007 + 0.4*F{walls}')
    minimum = gmsh.model.mesh.field.add('Min')
    gmsh.model.mesh.field.setNumbers(minimum, 'FieldsList', [near,wake,wallsize])
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize('Relocate2D')
    gmsh.write(str(out / 'source-gmsh.msh'))
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    # Gmsh retains the construction centre, which is not a fluid node: compact
    # only actual triangle nodes so both importers receive zero unused nodes.
    coords = {int(t):(float(xyz[3*i]),float(xyz[3*i+1])) for i,t in enumerate(tags)}
    types, _, blocks = gmsh.model.mesh.getElements(2)
    raw = []
    for kind, block in zip(types,blocks):
        if kind != 2: raise ValueError(f'Nonlinear/nontriangle element: {kind}')
        raw += [list(map(int,block[i:i+3])) for i in range(0,len(block),3)]
    used = sorted({t for tri in raw for t in tri})
    ids = {t:i+1 for i,t in enumerate(used)}
    nodes = [coords[t] for t in used]
    cells, areas = [], []
    for oldtri in raw:
        tri = [ids[t] for t in oldtri]
        aa,bb,cc = [nodes[t-1] for t in tri]
        area=((bb[0]-aa[0])*(cc[1]-aa[1])-(bb[1]-aa[1])*(cc[0]-aa[0]))/2
        if area < 0: tri[1],tri[2]=tri[2],tri[1]
        if area == 0: raise ValueError('Degenerate triangle')
        cells.append(tri);areas.append(abs(area))
    boundary={}
    for zone,curves in [(4,arcs),(5,[outer[3]]),(6,[outer[1]]),(7,[outer[0],outer[2]])]:
        for curve in curves:
            kinds,_,els=gmsh.model.mesh.getElements(1,curve)
            for kind,block in zip(kinds,els):
                assert kind==1
                for i in range(0,len(block),2):
                    boundary[tuple(sorted(ids[int(t)] for t in block[i:i+2]))]=zone
    faces={}
    for cid,tri in enumerate(cells,1):
        for aa,bb in zip(tri,tri[1:]+tri[:1]):
            key=tuple(sorted((aa,bb)))
            if key in faces:
                assert faces[key][:2]==[bb,aa] and not faces[key][3]
                faces[key][3]=cid
            else: faces[key]=[aa,bb,cid,0]
    groups={i:[] for i in (3,4,5,6,7)}
    for key,face in faces.items(): groups[3 if face[3] else boundary[key]].append(face)
    assert len(nodes)-len(faces)+len(cells)==0  # one cylindrical hole
    assert len(groups[4])==128 and len(groups[5])==41 and len(groups[6])==41
    mesh=out/'vortex.msh'
    with mesh.open('w',encoding='ascii',newline='\n') as f:
        f.write('(0 "DFG-inspired Re100 transient cylinder in channel")\n(2 2)\n')
        f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
        for x,y in nodes:f.write(f'{x:.16g} {y:.16g}\n')
        f.write('))\n');f.write(f'(12 (0 1 {len(cells):x} 0))\n(12 (2 1 {len(cells):x} 1 1))\n')
        f.write(f'(13 (0 1 {len(faces):x} 0))\n');first=1
        for zone,bc in [(3,2),(4,3),(5,10),(6,5),(7,3)]:
            rows=groups[zone]
            f.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
            for row in rows:f.write(' '.join(f'{v:x}' for v in row)+'\n')
            f.write('))\n');first+=len(rows)
        for zone,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','cylinder'),(5,'velocity-inlet','inlet'),(6,'pressure-outlet','outlet'),(7,'wall','walls')]:
            f.write(f'(45 ({zone} {kind} {name})())\n')
    inlet=[]
    for n0,n1,_,_ in groups[5]:
        aa,bb=nodes[n0-1],nodes[n1-1]
        x,y=(aa[0]+bb[0])/2,(aa[1]+bb[1])/2
        inlet.append(dict(x=x,y=y,u=6*y*(.41-y)/.41**2,length=math.dist(aa,bb)))
    flux=sum(r['u']*r['length'] for r in inlet)
    (out/'inlet-face-profile.json').write_text(json.dumps(inlet,indent=2)+'\n')
    # Native Fluent point profile: actual inlet face centres, no UDF/compile.
    with (out/'parabolic-inlet.prof').open('w',encoding='ascii',newline='\n') as f:
        f.write(f'((parabolic-inlet point {len(inlet)})\n')
        for key in ('x','y','u'):
            f.write('('+key+'\n'+'\n'.join(f'{r[key]:.16g}' for r in inlet)+'\n)\n')
        f.write(')\n')
    setup=dict(model='DFG-inspired cylinder in channel, Re=100 visual transient pilot',
        cells=len(cells),nodes=len(nodes),mesh_type='linear triangles',diameter_m=.1,
        cylinder_centre_m=[.2,.2],channel_length_m=2.2,channel_height_m=.41,
        density_kg_m3=1,dynamic_viscosity_Pa_s=.001,kinematic_viscosity_m2_s=.001,
        inlet_mean_velocity_m_s=1,inlet_maximum_velocity_m_s=1.5,reference_depth_m=1,
        reference_area_m2=.1,Re_D=100,dt_s=.005,pilot_end_time_s=2,extension_end_time_s=8,
        snapshot_interval_s=.025,cylinder_faces=128,boundary_faces={str(k):len(v) for k,v in groups.items() if k!=3},
        minimum_triangle_area_m2=min(areas),fluid_polygon_area_m2=sum(areas),
        exact_circle_fluid_area_m2=2.2*.41-math.pi*.05**2,
        inlet_midpoint_discrete_flow_m3_s=flux,inlet_exact_flow_m3_s=.41,
        mesh_sha256=hashlib.sha256(mesh.read_bytes()).hexdigest(),gmsh_version=gmsh.__version__,
        unused_nodes=0,mesh_independence=False,physical_validation=False,
        benchmark_note='Geometry/properties/inlet follow DFG 2D-2; practical p0 outlets, mesh, numerics and short startup duration differ. No validation claim.')
    (out/'setup.json').write_text(json.dumps(setup,indent=2)+'\n',encoding='utf-8')
    (out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(setup,indent=2))
finally:
    gmsh.finalize()
