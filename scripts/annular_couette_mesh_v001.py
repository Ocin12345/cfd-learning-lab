"""Fresh shared planar circular Couette mesh and source geometry; never launches CFD."""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np

RUN='2026-11-02-annular-couette-v001'
NT=64
NR=24
RI=.05
RO=.1


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


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
                if other not in seen:
                    seen.add(other);unseen.discard(other);stack.append(other)
        groups.append(seen)
    return groups


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    source=Path(__file__).resolve();lab=source.parents[3];out=a.out.resolve()
    require(out==(lab/'native/runs'/RUN/'mesh').resolve(),'Only exact fresh dedicated annular mesh directory allowed')
    require(out.drive.upper()=='E:' and not out.exists(),'Fresh external E: directory required')
    require(shutil.disk_usage(out.anchor).free>=10*1024**3,'Require10GiB free E:')
    radius_edges=np.linspace(RI,RO,NR+1);theta_edges=np.linspace(0.,2*np.pi,NT+1)
    theta=theta_edges[:-1];nodes=np.array([(r*np.cos(t),r*np.sin(t)) for r in radius_edges for t in theta])
    def nid(i,j):return j*NT+(i%NT)
    quads=np.array([[nid(i,j),nid(i,j+1),nid(i+1,j+1),nid(i+1,j)] for j in range(NR) for i in range(NT)],dtype=np.int64)
    polygons=nodes[quads];after=np.roll(polygons,-1,axis=1)
    cross=polygons[:,:,0]*after[:,:,1]-after[:,:,0]*polygons[:,:,1]
    areas=.5*cross.sum(axis=1);centres=((polygons+after)*cross[:,:,None]).sum(axis=1)/(6*areas[:,None])
    sector_index=np.tile(np.arange(NT,dtype=np.int64),NR);radial_index=np.repeat(np.arange(NR,dtype=np.int64),NT)
    require(nodes.shape==(1600,2) and quads.shape==(1536,4) and np.isfinite(nodes).all(),'Exact finite node/cell count')
    require(np.all(areas>0) and np.isfinite(centres).all(),'Positive oriented polygons and finite centroids')
    require(np.array_equal(np.unique(quads),np.arange(1600)),'Every node used without duplicate seam indices')
    require(len(np.unique(np.round(nodes,14),axis=0))==1600,'No duplicate seam nodes')
    require(np.max(abs(np.linalg.norm(nodes,axis=1)-np.repeat(radius_edges,NT)))<1e-15,'Exact concentric node radii')
    expected_area=.5*NT*np.sin(2*np.pi/NT)*(RO**2-RI**2);circle_area=np.pi*(RO**2-RI**2)
    require(abs(float(areas.sum())-expected_area)<1e-16,'Exact polygon annulus area')
    faces={}
    for cid,q in enumerate(quads):
        require(len(set(map(int,q)))==4,'Four unique vertices per source cell')
        for u,v in zip(q,np.roll(q,-1)):
            key=tuple(sorted((int(u),int(v))))
            if key in faces:
                row=faces[key];require(row[3]==-1 and row[0]==v and row[1]==u,'Shared edge occurs exactly twice with opposite directions');row[3]=cid
            else:faces[key]=[int(u),int(v),cid,-1]
    groups={3:[],4:[],5:[]};orth=[];cell_adj={i:set() for i in range(1536)};boundary_adj={}
    for u,v,owner,neighbour in faces.values():
        edge=nodes[v]-nodes[u];mid=.5*(nodes[u]+nodes[v]);normal=np.array([edge[1],-edge[0]])
        if neighbour>=0:
            zone=3;delta=centres[neighbour]-centres[owner];cell_adj[owner].add(neighbour);cell_adj[neighbour].add(owner)
        else:
            require(u//NT==v//NT and u//NT in (0,NR),'Only complete inner or outer ring edges are boundaries')
            zone=4 if u//NT==0 else 5;delta=mid-centres[owner]
            boundary_adj.setdefault(u,set()).add(v);boundary_adj.setdefault(v,set()).add(u)
            require((np.dot(normal,mid)<0) if zone==4 else (np.dot(normal,mid)>0),'Boundary normal points outward from fluid')
        cosine=float(np.dot(normal,delta)/(np.linalg.norm(normal)*np.linalg.norm(delta)))
        require(0<cosine<=1+1e-12,'Every face normal points owner to neighbour or outside fluid');orth.append(min(cosine,1.))
        groups[zone].append((u,v,owner,neighbour))
    require(len(faces)==3136 and [len(groups[z]) for z in (3,4,5)]==[3008,64,64],'Exact edge/zone coverage')
    require(len(nodes)-len(faces)+len(quads)==0,'Annular Euler characteristic0')
    require(len(components(cell_adj))==1,'One connected fluid component')
    loops=components(boundary_adj);require(len(loops)==2 and sorted(map(len,loops))==[64,64] and all(len(a)==2 for a in boundary_adj.values()),'Two simple closed boundary loops establish one hole')
    ordered=np.array([row for zone in (3,4,5) for row in groups[zone]],dtype=np.int64)
    face_nodes=ordered[:,:2];face_owner=ordered[:,2];face_neighbour=ordered[:,3]
    face_zone=np.concatenate([np.full(len(groups[z]),z,dtype=np.int64) for z in (3,4,5)])
    face_centres=nodes[face_nodes].mean(axis=1);edges=nodes[face_nodes[:,1]]-nodes[face_nodes[:,0]]
    face_area_vectors=np.column_stack([edges[:,1],-edges[:,0]]);face_lengths=np.linalg.norm(face_area_vectors,axis=1);face_normals=face_area_vectors/face_lengths[:,None]
    incidence=np.bincount(face_owner,minlength=1536)+np.bincount(face_neighbour[face_neighbour>=0],minlength=1536)
    require(np.all(incidence==4),'All1536 cells have exactly four incident faces')
    closure=np.zeros((1536,2));np.add.at(closure,face_owner,face_area_vectors);mask=face_neighbour>=0;np.add.at(closure,face_neighbour[mask],-face_area_vectors[mask])
    require(float(abs(closure).max())<1e-16,'Actual oriented face vectors close every polygon')
    boundary=np.flatnonzero(face_neighbour<0);reference=centres[0]
    require(RI<float(np.linalg.norm(reference))<RO,'Pressure reference source cell0 is in fluid, outside central hole')
    out.mkdir(parents=True);mesh=out/'annular-couette.msh'
    with mesh.open('x',encoding='ascii',newline='\n') as f:
        f.write('(0 "Planar annular circular Couette;64angular x24radial quads; source mesh only")\n(2 2)\n')
        f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
        for x,y in nodes:f.write(f'{x:.17g} {y:.17g}\n')
        f.write(f'))\n(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
        first=1
        for zone in (3,4,5):
            rows=groups[zone];bc=2 if zone==3 else 3
            f.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
            for u,v,o,n in rows:f.write(f'{u+1:x} {v+1:x} {o+1:x} {n+1:x}\n')
            f.write('))\n');first+=len(rows)
        for zone,kind,name in ((2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','inner'),(5,'wall','outer')):f.write(f'(45 ({zone} {kind} {name})())\n')
    np.savez_compressed(out/'geometry.npz',nodes=nodes,quads=quads,centres=centres,areas=areas,grid_shape=np.array([NR,NT]),
        ix=sector_index,iy=radial_index,sector_index=sector_index,radial_index=radial_index,radius_edges=radius_edges,theta_edges=theta_edges,
        r_centres=np.linalg.norm(centres,axis=1),theta_centres=np.mod(np.arctan2(centres[:,1],centres[:,0]),2*np.pi),
        face_nodes=face_nodes,face_owner=face_owner,face_neighbour=face_neighbour,face_zone=face_zone,
        face_centres=face_centres,face_normals=face_normals,face_area_vectors=face_area_vectors,face_lengths=face_lengths,
        face_source_ids=np.arange(1,3137,dtype=np.int64),boundary_face_source_ids=boundary+1,boundary_face_owner=face_owner[boundary],
        boundary_face_zone=face_zone[boundary],boundary_face_nodes=face_nodes[boundary],boundary_face_centres=face_centres[boundary],
        boundary_face_normals=face_normals[boundary],boundary_face_area_vectors=face_area_vectors[boundary],boundary_face_lengths=face_lengths[boundary],
        boundary_face_sector_index=sector_index[face_owner[boundary]],pressure_reference_source_cell=np.array(0),pressure_reference_xy_m=reference)
    proof=dict(nodes=1600,cells=1536,edges=3136,interior_faces=3008,wall_faces=dict(inner=64,outer=64),positive_areas=True,
        all_nodes_used=True,duplicate_seam_nodes=False,all_shared_edges_oppositely_oriented=True,all_cells_four_faces=True,
        connected_fluid_components=1,boundary_closed_loops=2,boundary_loop_nodes=[64,64],holes=1,euler_characteristic=0,
        total_polygon_area_m2=float(areas.sum()),exact_polygon_formula_area_m2=float(expected_area),circular_annulus_area_m2=float(circle_area),
        polygon_vs_circle_relative_area_error=float((areas.sum()-circle_area)/circle_area),minimum_area_m2=float(areas.min()),maximum_area_m2=float(areas.max()),
        orthogonality_proxy_min=float(min(orth)),maximum_nonorthogonality_proxy_deg=float(np.degrees(np.arccos(min(orth)))),
        maximum_cell_area_vector_closure_m=float(abs(closure).max()),solver_mesh_quality_not_yet_evaluated=True)
    geometry_path=out/'geometry.npz';source_sha=sha(source)
    setup=dict(model='Planar circular Couette spin-up; fixed inertial mesh; generic incompressible laminar fluid',nx=NT,ny=NR,ntheta=NT,nr=NR,cells=1536,
        inner_radius_m=RI,outer_radius_m=RO,inner_omega_rad_s=1.,omega_rad_s=1.,outer_omega_rad_s=0.,rotation_axis_origin_m=[0.,0.,0.],rotation_axis_direction=[0.,0.,1.],
        inner_wall_speed_m_s=.05,density_kg_m3=1.,dynamic_viscosity_Pa_s=.000125,kinematic_viscosity_m2_s=.000125,Re_gap=20.,
        reference_depth_m=1.,openfoam_extrusion_depth_m=.01,initial_velocity_m_s=[0.,0.],initial_pressure_Pa=0.,operating_pressure_Pa=0.,
        pressure_reference_source_cell=0,pressure_reference_xy_m=reference.tolist(),pressure_reference_xyz_m=[float(reference[0]),float(reference[1]),0.],
        time_step_s=.05,pilot_time_steps=5,pilot_final_time_s=.25,time_steps=400,final_time_s=20.,saved_state_count=401,animation_step_indices=list(range(0,401,4)),
        boundary_names=['inner','outer'],boundary_zone_ids={'interior':3,'inner':4,'outer':5},source_cell_order='radial_index*64+sector_index',
        source_node_order='radial_edge_index*64+sector_index, sector wrap; no duplicated seam nodes',face_order='interior, inner, outer; Fluent source ids1-based',
        geometry_index_convention='All NPZ node/cell/face-owner indices0-based; face_neighbour=-1 atboundary; *_source_ids are Fluent1-based face ids',
        energy=False,gravity=False,turbulence=False,moving_mesh=False,MRF=False,physical_validation=False,mesh_independence=False,time_step_independence=False,
        geometry_checks=proof,mesh_sha256=sha(mesh),geometry_npz_sha256=sha(geometry_path),mesher_source_sha256=source_sha,generated_by='Codex; no computed flow')
    json_new(out/'setup.json',setup);json_new(out/'mesh-proof.json',proof)
    (out/'executed-mesher.py').write_bytes(source.read_bytes());require(sha(out/'executed-mesher.py')==source_sha,'Exact executed source archive')
    manifest={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    json_new(out/'source-manifest.json',dict(kind='Source mesh/geometry only; no simulation',source_path=str(source),source_sha256=source_sha,files=manifest))
    print(json.dumps(dict(out=str(out),cells=1536,nodes=1600,faces=3136,mesh_sha256=sha(mesh),geometry_sha256=sha(geometry_path),mesher_sha256=source_sha,
        setup_sha256=sha(out/'setup.json'),pressure_reference_xy_m=reference.tolist(),geometry_checks=proof),indent=2))


if __name__=='__main__':main()