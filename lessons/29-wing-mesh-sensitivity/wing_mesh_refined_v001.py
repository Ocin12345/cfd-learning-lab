"""Refine the SAME wing using global mesh-size factor 0.9; no CFD launch.

Outputs are fresh E: folders only. The capped profile is NACA0012-like, not an
exact NACA0012. Positive linear tetrahedra are shared in Fluent ASCII, Gmsh2.2
and explicit NumPy connectivity. Root must review the mesh before solving.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import shutil
import tempfile

DOMAIN = [-3., 9., -4., 4., -6., 6.]
ZONES = {2: ('fluid','fluid'), 3: ('interior','interior'), 4: ('wall','wing'),
         5: ('velocity-inlet','inlet'), 6: ('pressure-outlet','outlet'),
         7: ('velocity-inlet','farfield')}


def write_json(path, values):
    path.write_text(json.dumps(values, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):h.update(block)
    return h.hexdigest()


def thickness(x):
    return .6*(.2969*math.sqrt(x)-.1260*x-.3516*x*x+.2843*x**3-.1015*x**4)


def rotate_point(x, y, z=0.):
    a=math.radians(8.)
    return (.25+(x-.25)*math.cos(a)+y*math.sin(a),
            -(x-.25)*math.sin(a)+y*math.cos(a),z)


def profile_data(count=49):
    cut=.96; radius=thickness(cut); scale=1/(cut+radius)
    upper=[]
    for i in range(count):
        x=cut*.5*(1-math.cos(math.pi*i/(count-1)))
        upper.append([x*scale,thickness(x)*scale])
    lower=[[x,-y] for x,y in reversed(upper)]
    cap=[[scale*(cut+radius*math.cos(t)),scale*radius*math.sin(t)]
         for t in [math.pi/2-math.pi*i/32 for i in range(33)]]
    return dict(upper=upper,lower=lower,cap=cap,
        original_truncation_x=cut,original_cap_radius=radius,scale_to_chord_one=scale,
        scaled_cap_center=[cut*scale,0.],scaled_cap_radius=radius*scale,
        profile_xy=upper+cap[1:]+lower[1:],
        construction='Interpolating OCC upper/lower splines; trailing semicircle; total chord scaled to1')


def positive_tets(points, cells):
    import numpy as np
    cells=np.asarray(cells,dtype=np.int64).copy()
    p=points[cells]
    volume=np.einsum('ij,ij->i',np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),p[:,3]-p[:,0])/6
    if not np.isfinite(volume).all() or np.any(volume==0):raise ValueError('Degenerate/nonfinite tetrahedron')
    negative=volume<0
    cells[negative,1],cells[negative,2]=cells[negative,2].copy(),cells[negative,1].copy()
    return cells,np.abs(volume),int(negative.sum())



def native_tet_volume_quality(points,cells):
    """One minus equilateral volume skewness; matched v001 Fluent's .02975645.

    q=V/[(8/(9*sqrt(3)))*R_circum**3]. Gamma/minSICN do not substitute for q.
    Circumcentres are solved relative to node0 to avoid large-origin cancellation.
    """
    import numpy as np
    tetra=points[cells]
    edges=tetra[:,1:]-tetra[:,:1]
    squared=np.einsum('nij,nij->ni',edges,edges)
    centre=np.linalg.solve(2.*edges,squared[:,:,None])[:,:,0]
    radius=np.linalg.norm(centre,axis=1)
    volume=np.abs(np.einsum('ij,ij->i',np.cross(edges[:,0],edges[:,1]),edges[:,2]))/6.
    quality=volume/((8./(9.*math.sqrt(3.)))*radius**3)
    if not np.isfinite(quality).all() or np.any(quality<=0) or np.any(quality>1.+1e-8):
        raise ValueError('Invalid native-equivalent tetrahedron volume quality')
    return quality,radius


def face_connectivity(points,cells,boundary):
    """ASCII Fluent requires the right-hand node normal toward c0 (not CFF order)."""
    import numpy as np
    faces={}
    for cid,(a,b,c,d) in enumerate(cells):
        for tri in ((b,d,c),(a,c,d),(a,d,b),(a,b,c)):
            key=tuple(sorted(map(int,tri)))
            if key in faces:
                if faces[key][4]!=-1:raise ValueError('Nonmanifold shared face')
                faces[key][4]=cid
            else:faces[key]=[*map(int,tri),cid,-1]
    groups={zone:[] for zone in (3,4,5,6,7)}
    for key,face in faces.items():
        zone=3 if face[4]>=0 else boundary.get(key)
        if zone is None:raise ValueError('Unclassified exterior tetrahedron face')
        groups[zone].append(face)
    if sum(len(groups[z]) for z in (4,5,6,7))!=len(boundary):raise ValueError('Boundary/facet count mismatch')
    if any(not groups[zone] for zone in groups):raise ValueError('Empty face group')
    centres=points[cells].mean(axis=1)
    minimum=1.
    for zone,rows in groups.items():
        values=np.asarray(rows,dtype=np.int64)
        tri=points[values[:,:3]]; centre=tri.mean(axis=1)
        normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        norm=np.linalg.norm(normal,axis=1)
        toward_c0=centres[values[:,3]]-centre
        dot=np.einsum('ij,ij->i',normal,toward_c0)
        if np.any(dot<=0):raise ValueError('Fluent ASCII face orientation fails c0 normal test')
        minimum=min(minimum,float(np.min(dot/(norm*np.linalg.norm(toward_c0,axis=1)))))
        if zone==3:
            delta=centres[values[:,3]]-centres[values[:,4]]
            metric=np.einsum('ij,ij->i',normal,delta)/(norm*np.linalg.norm(delta,axis=1))
            if np.any(metric<=0):raise ValueError('Interior face orientation/centres invalid')
            minimum=min(minimum,float(metric.min()))
    return groups,centres,minimum


def write_fluent(path,points,cells,groups):
    with path.open('w',encoding='ascii',newline='\n') as f:
        f.write('(0 "AR4 rounded-profile finite wing, steady laminar Re400")\n(2 3)\n')
        f.write(f'(10 (0 1 {len(points):x} 0))\n(10 (1 1 {len(points):x} 1 3)(\n')
        for row in points:f.write(' '.join(f'{v:.16g}' for v in row)+'\n')
        f.write('))\n')
        f.write(f'(12 (0 1 {len(cells):x} 0))\n(12 (2 1 {len(cells):x} 1 2))\n')
        f.write(f'(13 (0 1 {sum(map(len,groups.values())):x} 0))\n');first=1
        for zone,bc in ((3,2),(4,3),(5,10),(6,5),(7,10)):
            rows=groups[zone]
            f.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 3)(\n')
            for a,b,c,c0,c1 in rows:f.write(' '.join(f'{v:x}' for v in (a+1,b+1,c+1,c0+1,c1+1))+'\n')
            f.write('))\n');first+=len(rows)
        for zone,(kind,name) in ZONES.items():f.write(f'(45 ({zone} {kind} {name})())\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();out=a.out.resolve()
    if out.exists() or out.drive.upper()!='E:':raise ValueError('Use a fresh E: output folder')
    if shutil.disk_usage(out.anchor).free<10*1024**3:raise RuntimeError('E: free space below10GiB')
    out.mkdir(parents=True)
    temp=out/'temp';temp.mkdir()
    for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temp)
    tempfile.tempdir=str(temp)
    (out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
    profile=profile_data();write_json(out/'profile.json',profile)
    import numpy as np
    import gmsh
    gmsh.initialize()
    gmsh.logger.start()
    try:
        gmsh.option.setNumber('General.NumThreads',2)
        for key,value in {'Mesh.Algorithm':6,'Mesh.Algorithm3D':1,'Mesh.ElementOrder':1,
            'Mesh.MeshSizeFromPoints':0,'Mesh.MeshSizeFromCurvature':0,
            'Mesh.MeshSizeExtendFromBoundary':0,'Mesh.MeshSizeMin':.003,
            'Mesh.MeshSizeMax':1.3,'Mesh.MeshSizeFactor':0.9,'Mesh.Optimize':1,'Mesh.OptimizeThreshold':.35,
            'Mesh.RandomSeed':7,'Mesh.MshFileVersion':2.2,'Mesh.Binary':0}.items():
            gmsh.option.setNumber(key,value)
        gmsh.model.add('finite-wing-rounded-Re400')
        occ=gmsh.model.occ
        upper=[occ.addPoint(x,y,-2.) for x,y in profile['upper']]
        # Reuse shared leading/trailing points rather than coincident vertices.
        lower=[upper[-1]]  # placeholder replaced with lower trailing point
        lower=[occ.addPoint(profile['lower'][0][0],profile['lower'][0][1],-2.)]+[
            occ.addPoint(x,y,-2.) for x,y in profile['lower'][1:-1]]+[upper[0]]
        centre=occ.addPoint(*profile['scaled_cap_center'],-2.)
        tip=occ.addPoint(1.,0.,-2.)
        curves=[occ.addSpline(upper),occ.addCircleArc(upper[-1],centre,tip),
                occ.addCircleArc(tip,centre,lower[0]),occ.addSpline(lower)]
        surface=occ.addPlaneSurface([occ.addCurveLoop(curves)])
        extrusion=occ.extrude([(2,surface)],0,0,4.)
        wing_volume=next(tag for dim,tag in extrusion if dim==3)
        occ.rotate([(3,wing_volume)],.25,0,0,0,0,1,-math.radians(8.))
        box=occ.addBox(-3,-4,-6,12,8,12)
        fluid,_=occ.cut([(3,box)],[(3,wing_volume)],removeObject=True,removeTool=True)
        occ.synchronize()
        if len(fluid)!=1:raise ValueError('Expected one connected fluid CAD volume')
        volume=fluid[0][1]
        entities=gmsh.model.getBoundary([(3,volume)],oriented=False,recursive=False)
        surfaces={4:[],5:[],6:[],7:[]}
        for dim,tag in entities:
            x0,y0,z0,x1,y1,z1=gmsh.model.getBoundingBox(dim,tag)
            tol=1e-5
            if abs(x0+3)<tol and abs(x1+3)<tol:zone=5
            elif abs(x0-9)<tol and abs(x1-9)<tol:zone=6
            elif ((abs(y0+4)<tol and abs(y1+4)<tol) or (abs(y0-4)<tol and abs(y1-4)<tol)
                 or (abs(z0+6)<tol and abs(z1+6)<tol) or (abs(z0-6)<tol and abs(z1-6)<tol)):zone=7
            else:zone=4
            surfaces[zone].append(tag)
        if len(surfaces[5])!=1 or len(surfaces[6])!=1 or len(surfaces[7])!=4:raise ValueError('Outer CAD boundary classification failed')
        for zone,tags in surfaces.items():gmsh.model.addPhysicalGroup(2,tags,zone,ZONES[zone][1])
        gmsh.model.addPhysicalGroup(3,[volume],2,'fluid')
        distance=gmsh.model.mesh.field.add('Distance')
        gmsh.model.mesh.field.setNumbers(distance,'SurfacesList',surfaces[4])
        gmsh.model.mesh.field.setNumber(distance,'Sampling',80)
        near=gmsh.model.mesh.field.add('MathEval')
        gmsh.model.mesh.field.setString(near,'F',f'0.035+0.40*F{distance}')
        fields=[near]
        # The v001 .035 surface size under-resolved the thin round cap across
        # the full span. Resolve the SAME two cylindrical cap surfaces locally;
        # other wing, tip-wake and far-domain size fields remain unchanged.
        cap_surfaces=[tag for tag in surfaces[4] if gmsh.model.getType(2,tag)=='Cylinder']
        if len(cap_surfaces)!=2:raise ValueError('Expected the two unchanged round-cap CAD surfaces')
        cap_distance=gmsh.model.mesh.field.add('Distance')
        gmsh.model.mesh.field.setNumbers(cap_distance,'SurfacesList',cap_surfaces)
        gmsh.model.mesh.field.setNumber(cap_distance,'Sampling',400)
        cap_size=gmsh.model.mesh.field.add('MathEval')
        gmsh.model.mesh.field.setString(cap_size,'F',f'0.0075+0.32*F{cap_distance}')
        fields.append(cap_size)
        for values in [dict(VIn=.18,VOut=1.3,XMin=.1,XMax=6.,YMin=-.8,YMax=.45,ZMin=-2.5,ZMax=2.5,Thickness=.4),
                       dict(VIn=.13,VOut=1.3,XMin=.2,XMax=6.,YMin=-.7,YMax=.35,ZMin=1.6,ZMax=2.5,Thickness=.3),
                       dict(VIn=.13,VOut=1.3,XMin=.2,XMax=6.,YMin=-.7,YMax=.35,ZMin=-2.5,ZMax=-1.6,Thickness=.3)]:
            field=gmsh.model.mesh.field.add('Box')
            for key,value in values.items():gmsh.model.mesh.field.setNumber(field,key,value)
            fields.append(field)
        minimum=gmsh.model.mesh.field.add('Min')
        gmsh.model.mesh.field.setNumbers(minimum,'FieldsList',fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(minimum)
        gmsh.write(str(out/'geometry.brep'))
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize('Netgen')
        kinds,element_tags,blocks=gmsh.model.mesh.getElements(3)
        if list(kinds)!=[4]:raise ValueError('Expected only linear tetrahedra')
        raw=np.asarray(blocks[0],dtype=np.int64).reshape(-1,4)
        if not 1<=len(raw)<=250000:raise ValueError(f'Tetrahedron hard cell cap exceeded: {len(raw)}')
        tags,xyz,_=gmsh.model.mesh.getNodes()
        lookup={int(tag):xyz[3*i:3*i+3] for i,tag in enumerate(tags)}
        used=np.unique(raw);ids={int(tag):i for i,tag in enumerate(used)}
        points=np.asarray([lookup[int(tag)] for tag in used],dtype=float)
        cells=np.asarray([[ids[int(tag)] for tag in row] for row in raw],dtype=np.int64)
        cells,volumes,flipped=positive_tets(points,cells)
        boundary={}
        for zone,tags in surfaces.items():
            for tag in tags:
                kinds,_,blocks=gmsh.model.mesh.getElements(2,tag)
                for kind,block in zip(kinds,blocks):
                    if kind!=2:raise ValueError('Nonlinear/nontriangular boundary element')
                    for tri in np.asarray(block).reshape(-1,3):
                        key=tuple(sorted(ids[int(v)] for v in tri))
                        if key in boundary:raise ValueError('Duplicated exterior facet')
                        boundary[key]=zone
        groups,centres,proxy=face_connectivity(points,cells,boundary)
        native_quality,circumradius=native_tet_volume_quality(points,cells)
        candidate_floor=.105  # .10 native gate plus an absolute .005 margin.
        quality_pass=bool(float(native_quality.min())>=candidate_floor and proxy>=candidate_floor)
        quality_receipt=dict(definition='V / [(8/(9sqrt3))*R_circum^3], one minus equilateral volume skewness',
            minimum=float(native_quality.min()),percentile1=float(np.percentile(native_quality,1)),
            median=float(np.median(native_quality)),count_below_native_0p10=int((native_quality<.10).sum()),
            count_below_candidate_0p105=int((native_quality<candidate_floor).sum()),
            face_centre_orthogonality_proxy_minimum=proxy,required_candidate_floor=candidate_floor,
            all_candidate_tets_checked=len(cells),candidate_gate_pass=quality_pass,
            matched_v001_native_minimum=.029756453880275018,native_fluent_quality_check_still_required=True,
            cap_surfaces=cap_surfaces,cap_surface_size_m=.0075,cap_growth=.32,cap_distance_sampling=400,
            geometry_changed=False,repair_scope='Same CAD and size fields, all mesh lengths multiplied by 0.9 for mesh sensitivity')
        write_json(out/'candidate-native-quality.json',quality_receipt)
        mesh=out/'wing.msh';write_fluent(mesh,points,cells,groups)
        gmsh.write(str(out/'source-gmsh.msh'))
        qualities={name:np.asarray(gmsh.model.mesh.getElementQualities(element_tags[0],name)) for name in ('gamma','minSICN')}
        wing=np.asarray(groups[4],dtype=np.int64)
        tri=points[wing[:,:3]];normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        # ASCII normals point into fluid c0, therefore out of the wing solid.
        normals/=np.linalg.norm(normals,axis=1)[:,None]
        np.savez_compressed(out/'geometry.npz',points_xyz_m=points,nodes=points,tets=cells,
            tet_centres=centres,tet_volumes=volumes,tet_native_volume_quality=native_quality,tet_circumradius=circumradius,wing_face_nodes=wing[:,:3],wing_faces=wing[:,:3],
            wing_face_owner=wing[:,3],wing_face_owner_cells=wing[:,3],wing_normals=normals,
            profile_xy=np.asarray(profile['profile_xy']),profile_xyz_incidence=np.asarray([rotate_point(x,y) for x,y in profile['profile_xy']]))
        with (out/'wing-surface.stl').open('w',encoding='ascii') as f:
            f.write('solid wing\n')
            for normal,vertices in zip(normals,tri):
                f.write(' facet normal '+' '.join(map(str,normal))+'\n  outer loop\n')
                for vertex in vertices:f.write('   vertex '+' '.join(map(str,vertex))+'\n')
                f.write('  endloop\n endfacet\n')
            f.write('endsolid wing\n')
        cad_volume=gmsh.model.occ.getMass(3,volume)
        setup=dict(kind='finite-wing-steady-Re400',geometry='Rounded capped NACA0012-like rectangular finite wing; not exact NACA0012',
            cells=len(cells),nodes=len(points),unused_nodes=0,mesh_type='positive linear tetrahedra',
            target_cells=[180000,250000],hard_maximum_cells=250000,target_cell_count_met=180000<=len(cells)<=250000,
            chord_m=1.,span_m=4.,reference_area_m2=4.,aspect_ratio=4.,incidence_degrees_clockwise=8.,
            rotation_axis='z',rotation_origin_m=[.25,0.,0.],span_limits_m=[-2.,2.],domain_bounds_m=DOMAIN,
            density_kg_m3=1.,dynamic_viscosity_Pa_s=.0025,kinematic_viscosity_m2_s=.0025,
            inlet_velocity_m_s=[1.,0.,0.],Re_c=400,Cp_reference_gauge_pressure_Pa=0.,
            profile_metadata='profile.json',geometry_npz='geometry.npz',npz_index_base=0,
            wing_normal_definition='unit outward normal from wing solid into fluid; Fluent ASCII right-hand normal toward c0',
            boundary_zones={str(zone):dict(name=ZONES[zone][1],native_type=ZONES[zone][0],faces=len(groups[zone]),cad_surfaces=surfaces.get(zone,[])) for zone in groups},
            minimum_tet_volume_m3=float(volumes.min()),maximum_tet_volume_m3=float(volumes.max()),
            mesh_fluid_volume_m3=float(volumes.sum()),cad_fluid_volume_m3=cad_volume,
            mesh_solid_volume_m3=1152.-float(volumes.sum()),cad_solid_volume_m3=1152.-cad_volume,
            orientation_repaired_tets=flipped,face_centre_orthogonality_proxy_minimum=proxy,
            gmsh_quality={name:dict(minimum=float(v.min()),percentile1=float(np.percentile(v,1)),median=float(np.median(v))) for name,v in qualities.items()},
            fluent_native_orthogonal_quality_floor=.10,fluent_native_check_pending=True,
            mesh_version='sensitivity-v001-factor0p9',same_geometry_as_v001=True,mesh_size_factor=0.9,candidate_quality_gate=quality_receipt,
            gmsh_version=gmsh.__version__,mesh_sha256=sha256(mesh),gmsh_mesh_sha256=sha256(out/'source-gmsh.msh'),
            geometry_npz_sha256=sha256(out/'geometry.npz'),mesher_source_sha256=sha256(Path(__file__)),
            mesh_independence=False,physical_validation=False,steady_visualization_only=True,
            geometry_sources=['https://gmsh.info/doc/texinfo/gmsh.html',
                'https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/tgd_user_format_gridsect.html'])
        write_json(out/'setup.json',setup)
        (out/'gmsh-log.txt').write_text('\n'.join(gmsh.logger.get())+'\n',encoding='utf-8')
        write_json(out/'mesh-source-manifest.json',dict(files=[dict(name=p.name,sha256=sha256(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()]))
        print(json.dumps(setup,indent=2),flush=True)
        return 0 if quality_pass else 2
    finally:
        gmsh.logger.stop();gmsh.finalize()


if __name__=='__main__':raise SystemExit(main())
