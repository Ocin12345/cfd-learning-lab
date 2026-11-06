"""Fresh source-only rotating-cylinder O-grid; no Fluent, WSL or CFD launch."""
from pathlib import Path
import argparse,ctypes,hashlib,json,shutil
import numpy as np

RUN='2026-11-02-magnus-cylinder-v001'
TARGET_SHA='87dfdeb12300ba0d6dfaefdedf58d0c08ef984673522f4432b1fe90807f0dbbe'
PLAN_SHA='0e811eb689652575cdb036560b65a2a3f3c0cc2abd604b0c6073bac8a390d49d'
NT=128
NR=64
RI=.05
RO=2.
ZONES=(3,4,5,6)
ZONE_NAMES={3:'interior',4:'cylinder',5:'inlet',6:'outlet'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(value,message):
    if not value:raise ValueError(message)


def json_new(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def components(adjacency):
    unseen=set(adjacency);groups=[]
    while unseen:
        start=min(unseen);unseen.remove(start);seen={start};stack=[start]
        while stack:
            for other in adjacency[stack.pop()]:
                if other not in seen:seen.add(other);unseen.discard(other);stack.append(other)
        groups.append(seen)
    return groups


def available_ram():
    class Memory(ctypes.Structure):
        _fields_=[('dwLength',ctypes.c_ulong),('dwMemoryLoad',ctypes.c_ulong)]+[(name,ctypes.c_ulonglong) for name in ('totalPhys','availPhys','totalPage','availPage','totalVirtual','availVirtual','availExtended')]
    m=Memory();m.dwLength=ctypes.sizeof(m);require(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))!=0,'Read actual available RAM')
    return int(m.availPhys)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    source=Path(__file__).resolve();lab=source.parents[3];private=lab/'local/magnus-cylinder-v001';out=args.out.resolve()
    require(sha(private/'TARGETS.md')==TARGET_SHA and sha(private/'EXECUTION-PLAN-v001.md')==PLAN_SHA,'Exact root targets and plan required')
    freeze=json.loads((private/'mesh-source-freeze-v002.json').read_text(encoding='utf-8-sig'));require(freeze['mesher_source_sha256']==sha(source) and freeze['scope']=='source_mesh_only_no_native_preflight_or_CFD','Mesher source frozen before mesh generation')
    require(out==(lab/'native/runs'/RUN/'mesh').resolve(),'Only exact fresh dedicated Magnus mesh output')
    require(out.drive.upper()=='E:' and not out.exists(),'Fresh external E mesh directory required')
    free_disk=shutil.disk_usage(out.anchor).free;free_ram=available_ram()
    require(free_disk>=10*1024**3,'Require10GiB free E:')
    require(free_ram>=512*1024**2,'Lightweight source-mesh preparation requires512MiB available RAM; no CFD is launched')
    theta_edges=np.linspace(0.,2*np.pi,NT+1);theta=theta_edges[:-1]
    radius_edges=RI+(RO-RI)*np.expm1(5*np.arange(NR+1)/NR)/np.expm1(5.)
    require(radius_edges[0]==RI and radius_edges[-1]==RO and np.all(np.diff(radius_edges)>0),'Exact increasing frozen radial grading')
    nodes=np.array([(r*np.cos(t),r*np.sin(t)) for r in radius_edges for t in theta])
    def nid(i,j):return j*NT+(i%NT)
    quads=np.array([[nid(i,j),nid(i,j+1),nid(i+1,j+1),nid(i+1,j)] for j in range(NR) for i in range(NT)],dtype=np.int64)
    polygons=nodes[quads];origin=polygons[:,:1].copy();local=polygons-origin;after=np.roll(local,-1,axis=1)
    cross=local[:,:,0]*after[:,:,1]-after[:,:,0]*local[:,:,1];areas=.5*cross.sum(axis=1)
    centres=origin[:,0]+((local+after)*cross[:,:,None]).sum(axis=1)/(6*areas[:,None])
    sector_index=np.tile(np.arange(NT,dtype=np.int64),NR);radial_index=np.repeat(np.arange(NR,dtype=np.int64),NT)
    require(nodes.shape==(8320,2) and quads.shape==(8192,4) and np.isfinite(nodes).all(),'Exact finite node/cell count')
    require(np.all(areas>0) and np.isfinite(centres).all(),'Positive polygons and finite exact polygon centroids')
    require(np.array_equal(np.unique(quads),np.arange(8320)) and len(np.unique(np.round(nodes,14),axis=0))==8320,'Every node used without duplicate seam')
    require(np.max(abs(np.linalg.norm(nodes,axis=1)-np.repeat(radius_edges,NT)))<2e-15,'Exact concentric node radii')
    polygon_area=.5*NT*np.sin(2*np.pi/NT)*(RO**2-RI**2);circle_area=np.pi*(RO**2-RI**2)
    require(abs(float(areas.sum())-polygon_area)<2e-13,'Exact inscribed polygon annulus area')
    faces={}
    for cid,q in enumerate(quads):
        for a,b in zip(q,np.roll(q,-1)):
            key=tuple(sorted((int(a),int(b))))
            if key in faces:
                row=faces[key];require(row[3]==-1 and row[0]==b and row[1]==a,'Exactly two opposite shared-edge incidences');row[3]=cid
            else:faces[key]=[int(a),int(b),cid,-1]
    groups={z:[] for z in ZONES};orth=[];cell_adj={i:set() for i in range(8192)};boundary_adj={}
    for a,b,owner,neighbour in faces.values():
        edge=nodes[b]-nodes[a];mid=.5*(nodes[a]+nodes[b]);normal=np.array([edge[1],-edge[0]])
        if neighbour>=0:
            zone=3;delta=centres[neighbour]-centres[owner];cell_adj[owner].add(neighbour);cell_adj[neighbour].add(owner)
        else:
            require(a//NT==b//NT and a//NT in (0,NR),'Only inner/outer rings have boundary edges')
            if a//NT==0:zone=4;require(np.dot(normal,mid)<0,'Cylinder fluid-outward normal points into hole')
            else:
                require(mid[0]!=0 and np.dot(normal,mid)>0,'Outer midpoint split and outward normal unambiguous')
                zone=5 if mid[0]<0 else 6
            delta=mid-centres[owner];boundary_adj.setdefault(a,set()).add(b);boundary_adj.setdefault(b,set()).add(a)
        cosine=float(np.dot(normal,delta)/(np.linalg.norm(normal)*np.linalg.norm(delta)))
        require(0<cosine<=1+1e-12,'Every normal points from owner outward');orth.append(min(cosine,1.))
        groups[zone].append((a,b,owner,neighbour))
    require(len(faces)==16512 and [len(groups[z]) for z in ZONES]==[16256,128,64,64],'Exact face and patch counts')
    require(len(nodes)-len(faces)+len(quads)==0 and len(components(cell_adj))==1,'Euler0 and one connected fluid component')
    loops=components(boundary_adj);require(len(loops)==2 and sorted(map(len,loops))==[128,128] and all(len(v)==2 for v in boundary_adj.values()),'Two simple boundary loops and one genuine hole')
    ordered=np.array([row for z in ZONES for row in groups[z]],dtype=np.int64)
    face_nodes=ordered[:,:2];face_owner=ordered[:,2];face_neighbour=ordered[:,3]
    face_zone=np.concatenate([np.full(len(groups[z]),z,dtype=np.int64) for z in ZONES])
    face_centres=nodes[face_nodes].mean(axis=1);edges=nodes[face_nodes[:,1]]-nodes[face_nodes[:,0]]
    face_area_vectors=np.column_stack([edges[:,1],-edges[:,0]]);face_lengths=np.linalg.norm(face_area_vectors,axis=1);face_normals=face_area_vectors/face_lengths[:,None]
    incidence=np.bincount(face_owner,minlength=8192)+np.bincount(face_neighbour[face_neighbour>=0],minlength=8192);require(np.all(incidence==4),'Every cell has four incident faces')
    closure=np.zeros((8192,2));np.add.at(closure,face_owner,face_area_vectors);m=face_neighbour>=0;np.add.at(closure,face_neighbour[m],-face_area_vectors[m]);require(abs(closure).max()<1e-14,'Actual cell face-area closure')
    inlet=face_zone==5;outlet=face_zone==6;body=face_zone==4;boundary=np.flatnonzero(face_neighbour<0)
    inlet_flow=float((face_area_vectors[inlet]@np.array([.1,0.])).sum());outlet_flow=float((face_area_vectors[outlet]@np.array([.1,0.])).sum())
    require(abs(inlet_flow+.4)<1e-14 and abs(outlet_flow-.4)<1e-14,'Exact polygon projected uniform-flow boundary integrals')
    omega=-2.;wall_velocity=omega*np.column_stack([-face_centres[body,1],face_centres[body,0]])
    wall_velocity-=face_normals[body]*np.sum(wall_velocity*face_normals[body],axis=1)[:,None]
    speeds=np.linalg.norm(wall_velocity,axis=1);expected_speed=abs(omega)*RI*np.cos(np.pi/NT)
    require(np.allclose(speeds,expected_speed,atol=1e-16,rtol=0) and abs(np.sum(wall_velocity*face_normals[body],axis=1)).max()<1e-16,'Projected clockwise facet-centre velocity geometry only')
    out.mkdir(parents=True);mesh=out/'magnus-cylinder.msh';ranges={};first=1
    with mesh.open('x',encoding='ascii',newline='\n') as f:
        f.write('(0 "Rotating-cylinder fixed O-grid;8192graded quads;source geometry only, no computed flow")\n(2 2)\n')
        f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
        for x,y in nodes:f.write(f'{x:.17g} {y:.17g}\n')
        f.write(f'))\n(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
        for z in ZONES:
            rows=groups[z];bc={3:2,4:3,5:10,6:5}[z];ranges[ZONE_NAMES[z]]=[first,first+len(rows)-1]
            f.write(f'(13 ({z:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
            for a,b,o,n in rows:f.write(f'{a+1:x} {b+1:x} {o+1:x} {n+1:x}\n')
            f.write('))\n');first+=len(rows)
        for z,kind,name in ((2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','cylinder'),(5,'velocity-inlet','inlet'),(6,'pressure-outlet','outlet')):f.write(f'(45 ({z} {kind} {name})())\n')
    np.savez_compressed(out/'geometry.npz',nodes=nodes,quads=quads,centres=centres,areas=areas,grid_shape=np.array([NR,NT]),ix=sector_index,iy=radial_index,
        sector_index=sector_index,radial_index=radial_index,radius_edges=radius_edges,theta_edges=theta_edges,r_centres=np.linalg.norm(centres,axis=1),theta_centres=np.mod(np.arctan2(centres[:,1],centres[:,0]),2*np.pi),
        face_nodes=face_nodes,face_owner=face_owner,face_neighbour=face_neighbour,face_zone=face_zone,face_centres=face_centres,face_normals=face_normals,face_area_vectors=face_area_vectors,face_lengths=face_lengths,
        face_source_ids=np.arange(1,16513,dtype=np.int64),boundary_face_source_ids=boundary+1,boundary_face_owner=face_owner[boundary],boundary_face_zone=face_zone[boundary],boundary_face_nodes=face_nodes[boundary],
        boundary_face_centres=face_centres[boundary],boundary_face_normals=face_normals[boundary],boundary_face_area_vectors=face_area_vectors[boundary],boundary_face_lengths=face_lengths[boundary],
        boundary_face_sector_index=sector_index[face_owner[boundary]],cylinder_face_source_ids=np.flatnonzero(body)+1,cylinder_projected_prescribed_wall_velocity=wall_velocity)
    proof=dict(cells=8192,nodes=8320,faces=16512,interior_faces=16256,physical_boundary_faces=256,boundary_face_counts={'cylinder':128,'inlet':64,'outlet':64},fluent_face_ranges_1based=ranges,
        positive_areas=True,minimum_area_m2=float(areas.min()),maximum_area_m2=float(areas.max()),total_polygon_area_m2=float(areas.sum()),exact_polygon_formula_area_m2=float(polygon_area),
        circular_annulus_area_m2=float(circle_area),polygon_vs_circle_relative_area_error=float((areas.sum()-circle_area)/circle_area),all_nodes_used=True,duplicate_seam_nodes=False,
        all_shared_edges_oppositely_oriented=True,all_cells_four_faces=True,connected_fluid_components=1,boundary_closed_loops=2,boundary_loop_nodes=[128,128],holes=1,euler_characteristic=0,
        orthogonality_proxy_min=float(min(orth)),maximum_nonorthogonality_proxy_deg=float(np.degrees(np.arccos(min(orth)))),maximum_cell_area_vector_closure_m=float(abs(closure).max()),
        minimum_edge_length_m=float(face_lengths.min()),maximum_edge_length_m=float(face_lengths.max()),uniform_inlet_outward_flow_per_depth_m2_s=inlet_flow,uniform_outlet_outward_flow_per_depth_m2_s=outlet_flow,
        prescribed_cylinder_facet_centre_speed_min_m_s=float(speeds.min()),prescribed_cylinder_facet_centre_speed_max_m_s=float(speeds.max()),prescribed_wall_velocity_is_geometry_only_not_native_CFD=True,
        native_solver_mesh_quality_pending=True,expected_openfoam_extrusion=dict(points=16640,faces=32896,internal_faces=16256,physical_faces=256,empty_faces=16384,depth_m=.01))
    setup=dict(model='Viscous laminar rotating-cylinder startup, generic incompressible fluid',cells=8192,nx=NT,ny=NR,ntheta=NT,nr=NR,
        diameter_m=.1,cylinder_radius_m=RI,inner_radius_m=RI,outer_radius_m=RO,radial_grading='0.05+1.95*expm1(5*j/64)/expm1(5)',
        density_kg_m3=1.,dynamic_viscosity_Pa_s=.0005,kinematic_viscosity_m2_s=.0005,inlet_velocity_m_s=[.1,0.],outlet_gauge_pressure_Pa=0.,operating_pressure_Pa=0.,
        Re_D=20.,Re_body=20.,reynolds_number=20.,cylinder_omega_rad_s=-2.,omega_rad_s=-2.,rotation_axis_origin_m=[0.,0.,0.],rotation_axis_direction=[0.,0.,1.],nominal_wall_speed_m_s=.1,
        nominal_spin_ratio=1.,prescribed_facet_centre_wall_speed_m_s=float(expected_speed),prescribed_wall_motion='Absolute rotational moving wall, fixed stationary cell zone/mesh',
        reference_length_m=.1,reference_velocity_m_s=.1,reference_area_m2=.1,reference_depth_m=1.,openfoam_extrusion_depth_m=.01,openfoam_z_bounds_m=[-.005,.005],
        force_coefficient_denominator_N=.0005,openfoam_force_coefficient_denominator_N=.000005,
        initial_velocity_m_s=[0.,0.],initial_pressure_Pa=0.,initial_state_is_prescribed_not_conservative_solved_flow=True,
        time_step_s=.02,pilot_time_steps=5,pilot_final_time_s=.1,time_steps=500,final_time_s=10.,saved_state_count=501,animation_step_indices=list(range(0,501,5)),
        boundary_names=['cylinder','inlet','outlet'],boundary_types={'cylinder':'wall','inlet':'velocity-inlet','outlet':'pressure-outlet'},boundary_zone_ids={ZONE_NAMES[z]:z for z in ZONES},boundary_face_counts=proof['boundary_face_counts'],
        source_cell_order='radial_index*128+sector_index; node=radial_edge_index*128+sector_index with angular wrap',face_order='interior,cylinder,inlet,outlet',
        geometry_index_convention='NPZ node/cell/owner indices0-based; face_neighbour=-1boundary; *_source_ids are Fluent1-based faceids',
        energy=False,turbulence=False,gravity=False,MRF=False,moving_mesh=False,physical_validation=False,mesh_independence=False,time_step_independence=False,
        geometry_checks=proof,mesh_sha256=sha(mesh),geometry_npz_sha256=sha(out/'geometry.npz'),mesher_source_sha256=sha(source),target_document_sha256=TARGET_SHA,execution_plan_sha256=PLAN_SHA,
        source_freeze_sha256=sha(private/'mesh-source-freeze-v002.json'),generated_by='Codex; source mesh only; no native preflight or CFD')
    json_new(out/'mesh-proof.json',proof);json_new(out/'setup.json',setup)
    json_new(out/'mesh-resource-receipt.json',dict(available_ram_bytes_before_mesh=free_ram,mesh_only_ram_floor_bytes=512*1024**2,free_E_bytes_before_mesh=free_disk,disk_floor_bytes=10*1024**3,
        CFD_RAM_floor_bytes_from_TARGETS=4*1024**3,CFD_or_native_preflight_launched=False,mesh_only=True,threaded_solver_launched=False))
    (out/'executed-mesher.py').write_bytes(source.read_bytes());require(sha(source)==sha(out/'executed-mesher.py'),'Executed source archive exact')
    json_new(out/'source-manifest.json',dict(kind='Frozen source mesh only, no native preflight or CFD',source_path=str(source),source_sha256=sha(source),target_document_sha256=TARGET_SHA,execution_plan_sha256=PLAN_SHA,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}))
    print(json.dumps(dict(out=str(out),source_sha256=sha(source),mesh_sha256=sha(mesh),geometry_npz_sha256=sha(out/'geometry.npz'),setup_sha256=sha(out/'setup.json'),proof_sha256=sha(out/'mesh-proof.json'),geometry_checks=proof),indent=2))


if __name__=='__main__':main()