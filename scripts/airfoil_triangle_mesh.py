"""Remesh the existing airfoil polygon without rounding its trailing edge.

Gmsh is a mesher only. Fluent must independently check and solve the resulting
mesh. All boundary segments come from the original Fluent mesh's wall/outer rings.
"""
from pathlib import Path
import argparse,hashlib,json,math

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--scale',type=float,default=1.0)
    p.add_argument('--angle',type=float)
    p.add_argument('--local-edge-sizing',action='store_true')
    p.add_argument('--scale-edge-limits',action='store_true',
                   help='Apply size scale to the original-edge point limits as well as background sizes')
    args=p.parse_args()
    repo=Path(__file__).resolve().parents[1]
    if args.out.resolve().is_relative_to(repo) or args.out.exists():
        raise ValueError('Use a fresh external native output directory')
    if not math.isfinite(args.scale) or args.scale<=0:raise ValueError('Positive scale required')
    source_setup=json.loads((args.source.parent/'setup.json').read_text())
    nt=source_setup['circumferential_cells']
    lines=args.source.read_text().splitlines()
    start=next(i for i,s in enumerate(lines) if s.startswith('(10 (1 '))+1
    stop=next(i for i in range(start,len(lines)) if lines[i]=='))')
    source_nodes=[tuple(map(float,s.split())) for s in lines[start:stop]]
    angle=source_setup['angle_degrees'] if args.angle is None else args.angle
    delta=math.radians(angle-source_setup['angle_degrees'])
    def rotate(xy):
        x,y=xy
        return (.5+(x-.5)*math.cos(delta)+y*math.sin(delta),
                -(x-.5)*math.sin(delta)+y*math.cos(delta))
    wall=list(map(rotate,source_nodes[:nt]))
    outer=list(map(rotate,source_nodes[-nt:]))
    args.out.mkdir(parents=True)
    (args.out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.NumThreads',2)
        gmsh.option.setNumber('Mesh.Algorithm',6)
        gmsh.option.setNumber('Mesh.ElementOrder',1)
        gmsh.option.setNumber('Mesh.MeshSizeFromPoints',int(args.local_edge_sizing))
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',0)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',int(args.local_edge_sizing))
        gmsh.option.setNumber('Mesh.MeshSizeMin',1e-5*args.scale)
        gmsh.option.setNumber('Mesh.MeshSizeMax',1.0*args.scale)
        gmsh.option.setNumber('Mesh.RandomSeed',1)
        gmsh.model.add('preserved-airfoil-polygon')
        def loop(xys,body=False):
            sizes=[0.0]*nt
            if args.local_edge_sizing:
                edges=[math.dist(xys[i],xys[(i+1)%nt]) for i in range(nt)]
                cap=(.008 if body else 1.0)*args.scale
                edge_factor=args.scale if args.scale_edge_limits else 1.0
                sizes=[min(cap,.8*edge_factor*min(edges[i-1],edges[i])) for i in range(nt)]
            points=[gmsh.model.geo.addPoint(x,y,0,sizes[i]) for i,(x,y) in enumerate(xys)]
            curves=[gmsh.model.geo.addLine(points[i],points[(i+1)%nt]) for i in range(nt)]
            return curves,gmsh.model.geo.addCurveLoop(curves)
        body_lines,body_loop=loop(wall,body=True)
        far_lines,far_loop=loop(outer)
        surface=gmsh.model.geo.addPlaneSurface([far_loop,body_loop])
        gmsh.model.geo.synchronize()
        inlet=[line for i,line in enumerate(far_lines) if outer[i][0]+outer[(i+1)%nt][0]<1]
        outlet=[line for line in far_lines if line not in inlet]
        for tag,name,entities in [(4,'airfoil',body_lines),(5,'inlet',inlet),(6,'outlet',outlet)]:
            gmsh.model.addPhysicalGroup(1,entities,tag,name)
        gmsh.model.addPhysicalGroup(2,[surface],2,'fluid')
        distance=gmsh.model.mesh.field.add('Distance')
        gmsh.model.mesh.field.setNumbers(distance,'CurvesList',body_lines)
        gmsh.model.mesh.field.setNumber(distance,'Sampling',20)
        near=gmsh.model.mesh.field.add('MathEval')
        gmsh.model.mesh.field.setString(near,'F',f'{.008*args.scale} + {0.2*args.scale}*F{distance}')
        wake=gmsh.model.mesh.field.add('Box')
        for key,value in dict(VIn=.04*args.scale,VOut=1000,XMin=.5,XMax=10,YMin=-.5,YMax=.5,ZMin=-1,ZMax=1,Thickness=.5).items():
            gmsh.model.mesh.field.setNumber(wake,key,value)
        field=gmsh.model.mesh.field.add('Min')
        gmsh.model.mesh.field.setNumbers(field,'FieldsList',[near,wake])
        gmsh.model.mesh.field.setAsBackgroundMesh(field)
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize('Relocate2D')
        gmsh.write(str(args.out/'source-gmsh.msh'))
        tags,xyz,_=gmsh.model.mesh.getNodes()
        mapping={int(tag):i+1 for i,tag in enumerate(tags)}
        nodes=[(float(xyz[3*i]),float(xyz[3*i+1])) for i in range(len(tags))]
        types,_,blocks=gmsh.model.mesh.getElements(2)
        triangles=[]
        min_area=math.inf
        for kind,block in zip(types,blocks):
            if kind!=2:raise ValueError(f'Unexpected non-linear-triangle element type {kind}')
            for i in range(0,len(block),3):
                tri=[mapping[int(tag)] for tag in block[i:i+3]]
                a,b,c=[nodes[index-1] for index in tri]
                area=((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
                if area<0:tri[1],tri[2]=tri[2],tri[1]
                min_area=min(min_area,abs(area))
                triangles.append(tri)
        assert min_area>0
        boundary={}
        for group,entities in [(4,body_lines),(5,inlet),(6,outlet)]:
            for entity in entities:
                line_types,_,line_blocks=gmsh.model.mesh.getElements(1,entity)
                for kind,block in zip(line_types,line_blocks):
                    assert kind==1
                    for i in range(0,len(block),2):
                        boundary[tuple(sorted(mapping[int(tag)] for tag in block[i:i+2]))]=group
        faces={}
        for cid,tri in enumerate(triangles,1):
            for a,b in zip(tri,tri[1:]+tri[:1]):
                key=tuple(sorted((a,b)))
                if key in faces:
                    assert faces[key][:2]==[b,a] and not faces[key][3]
                    faces[key][3]=cid
                else:faces[key]=[a,b,cid,0]
        groups={3:[],4:[],5:[],6:[]}
        for key,face in faces.items():
            zone=3 if face[3] else boundary[key]
            groups[zone].append(face)
        assert sum(len(groups[z]) for z in [4,5,6])==len(boundary)
        with (args.out/'airfoil.msh').open('w') as f:
            f.write('(0 "Airfoil triangular remesh; original boundary polygon preserved")\n(2 2)\n')
            f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
            for x,y in nodes:f.write(f'{x:.14g} {y:.14g}\n')
            f.write('))\n')
            f.write(f'(12 (0 1 {len(triangles):x} 0))\n(12 (2 1 {len(triangles):x} 1 1))\n')
            f.write(f'(13 (0 1 {len(faces):x} 0))\n')
            first=1
            for zone,bc in [(3,2),(4,3),(5,10),(6,5)]:
                ff=groups[zone]
                f.write(f'(13 ({zone:x} {first:x} {first+len(ff)-1:x} {bc:x} 2)(\n')
                for face in ff:f.write(' '.join(f'{v:x}' for v in face)+'\n')
                f.write('))\n')
                first+=len(ff)
            for zone,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','airfoil'),(5,'velocity-inlet','inlet'),(6,'pressure-outlet','outlet')]:
                f.write(f'(45 ({zone} {kind} {name})())\n')
        coords={(round(x,12),round(y,12)) for x,y in nodes}
        assert all((round(x,12),round(y,12)) in coords for x,y in wall+outer)
        setup=dict(chord_m=1,angle_degrees=angle,geometry='original scaled closed NACA0012 polygon',
                   outer_radius_m=20,density_kg_m3=1,viscosity_Pa_s=.01,inlet_velocity_m_s=1,
                   reynolds_number=100,cells=len(triangles),nodes=len(nodes),mesh_type='triangles',
                   size_scale=args.scale,local_edge_sizing=args.local_edge_sizing,scale_edge_limits=args.scale_edge_limits,
                   preserved_source_polygon_vertices=True,
                   source_mesh_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                   minimum_triangle_area_m2=min_area,gmsh_version=gmsh.__version__)
        (args.out/'setup.json').write_text(json.dumps(setup,indent=2))
        print(json.dumps(setup,indent=2),flush=True)
    finally:
        gmsh.finalize()

if __name__=='__main__':main()
