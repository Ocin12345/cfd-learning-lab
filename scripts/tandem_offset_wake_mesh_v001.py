"""Shared fixed-grid offset tandem-square source mesh; no solver or WSL launch."""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np

RUN='2026-11-02-tandem-offset-wake-v001'
TARGET_SHA='426448b7e6909e75dbf08f065d93f90699b98f12be007891631d9ed9b10f71d7'
NX=160
NY=60
DX=.0125
DY=.0125
Y0=-.375
BOXES=((36,44,24,32),(68,76,28,36))
ZONES=(3,4,5,6,7,8,9)
ZONE_NAMES={3:'interior',4:'inlet',5:'outlet',6:'top',7:'bottom',8:'body1',9:'body2'}


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


def on_box_edge(a,b,box):
    x0,x1,y0,y1=box
    return ((a[0]==b[0] and a[0] in (x0,x1) and y0<=a[1]<=y1 and y0<=b[1]<=y1)
        or (a[1]==b[1] and a[1] in (y0,y1) and x0<=a[0]<=x1 and x0<=b[0]<=x1))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    source=Path(__file__).resolve();lab=source.parents[3];out=args.out.resolve();targets=lab/'local/tandem-offset-wake-v001/TARGETS.md'
    require(sha(targets)==TARGET_SHA,'Exact root frozen targets required before mesh generation')
    require(out==(lab/'native/runs'/RUN/'mesh').resolve(),'Only exact dedicated tandem mesh output allowed')
    require(out.drive.upper()=='E:' and not out.exists(),'Fresh external E mesh directory required')
    require(shutil.disk_usage(out.anchor).free>=10*1024**3,'Require10GiB free E:')
    body1_mask=np.zeros((NY,NX),bool);body2_mask=body1_mask.copy()
    for mask,(x0,x1,y0,y1) in zip((body1_mask,body2_mask),BOXES):mask[y0:y1,x0:x1]=True
    require(not np.any(body1_mask&body2_mask),'Two distinct non-overlapping holes')
    solid_mask=body1_mask|body2_mask;fluid_mask=~solid_mask
    cells=[(i,j) for j in range(NY) for i in range(NX) if fluid_mask[j,i]]
    vertices=sorted({v for i,j in cells for v in ((i,j),(i+1,j),(i+1,j+1),(i,j+1))},key=lambda v:(v[1],v[0]))
    numbering={v:k for k,v in enumerate(vertices)};grid_vertices=np.asarray(vertices,dtype=np.int64)
    nodes=np.column_stack([grid_vertices[:,0]*DX,Y0+grid_vertices[:,1]*DY])
    quads=np.array([[numbering[v] for v in ((i,j),(i+1,j),(i+1,j+1),(i,j+1))] for i,j in cells],dtype=np.int64)
    polygons=nodes[quads];origin=polygons[:,:1].copy();local=polygons-origin;after=np.roll(local,-1,axis=1)
    cross=local[:,:,0]*after[:,:,1]-after[:,:,0]*local[:,:,1];areas=.5*cross.sum(axis=1)
    centres=origin[:,0]+((local+after)*cross[:,:,None]).sum(axis=1)/(6*areas[:,None])
    slots=np.asarray(cells,dtype=np.int64);ix,iy=slots.T
    require(nodes.shape==(9723,2) and quads.shape==(9472,4) and np.isfinite(nodes).all(),'Exact finite node/cell count')
    require(np.all(areas>0) and np.isfinite(centres).all(),'Positive oriented source polygons')
    require(np.allclose(areas,DX*DY,atol=1e-16,rtol=0) and abs(float(areas.sum())-1.48)<1e-13,'Uniform cells and exact two-hole fluid area')
    require(np.allclose(centres,np.column_stack([(ix+.5)*DX,Y0+(iy+.5)*DY]),rtol=0,atol=1e-13),'Exact slot-to-centroid mapping')
    require(np.array_equal(np.unique(quads),np.arange(9723)) and len(np.unique(grid_vertices,axis=0))==9723,'Every unique grid node used; no seam duplicates')
    slot_to_source_cell=np.full((NY,NX),-1,dtype=np.int64);slot_to_source_cell[iy,ix]=np.arange(9472)
    source_cell_to_slot=iy*NX+ix
    require(np.array_equal(slot_to_source_cell.ravel()[source_cell_to_slot],np.arange(9472)),'Source-slot bijection')
    require(np.count_nonzero(body1_mask)==64 and np.count_nonzero(body2_mask)==64 and np.count_nonzero(slot_to_source_cell<0)==128,'Exactly two64-slot holes')
    faces={}
    for cid,q in enumerate(quads):
        for a,b in zip(q,np.roll(q,-1)):
            key=tuple(sorted((int(a),int(b))))
            if key in faces:
                row=faces[key];require(row[3]==-1 and row[0]==b and row[1]==a,'Every shared edge has exactly two opposite incidences');row[3]=cid
            else:faces[key]=[int(a),int(b),cid,-1]
    groups={z:[] for z in ZONES};orth=[];cell_adj={i:set() for i in range(9472)};boundary_adj={}
    for a,b,owner,neighbour in faces.values():
        edge=nodes[b]-nodes[a];mid=.5*(nodes[a]+nodes[b]);normal=np.array([edge[1],-edge[0]])
        if neighbour>=0:
            zone=3;delta=centres[neighbour]-centres[owner];cell_adj[owner].add(neighbour);cell_adj[neighbour].add(owner)
        else:
            delta=mid-centres[owner];a_grid,b_grid=vertices[a],vertices[b]
            if a_grid[0]==b_grid[0]==0:zone=4
            elif a_grid[0]==b_grid[0]==NX:zone=5
            elif a_grid[1]==b_grid[1]==NY:zone=6
            elif a_grid[1]==b_grid[1]==0:zone=7
            elif on_box_edge(a_grid,b_grid,BOXES[0]):zone=8
            elif on_box_edge(a_grid,b_grid,BOXES[1]):zone=9
            else:raise ValueError('Unclassified physical boundary edge')
            boundary_adj.setdefault(a,set()).add(b);boundary_adj.setdefault(b,set()).add(a)
        cosine=float(np.dot(normal,delta)/(np.linalg.norm(normal)*np.linalg.norm(delta)))
        require(0<cosine<=1+1e-12,'Owner-outward face orientations');orth.append(min(cosine,1.))
        groups[zone].append((a,b,owner,neighbour))
    require(len(faces)==19196 and [len(groups[z]) for z in ZONES]==[18692,60,60,160,160,32,32],'Exact global/per-zone face counts')
    require(len(nodes)-len(faces)+len(quads)==-1 and len(components(cell_adj))==1,'One connected fluid domain with Euler-1')
    loops=components(boundary_adj);require(len(loops)==3 and sorted(map(len,loops))==[32,32,440] and all(len(v)==2 for v in boundary_adj.values()),'Three simple closed boundary loops/two holes')
    ordered=np.array([row for z in ZONES for row in groups[z]],dtype=np.int64)
    face_nodes=ordered[:,:2];face_owner=ordered[:,2];face_neighbour=ordered[:,3]
    face_zone=np.concatenate([np.full(len(groups[z]),z,dtype=np.int64) for z in ZONES])
    face_centres=nodes[face_nodes].mean(axis=1);edges=nodes[face_nodes[:,1]]-nodes[face_nodes[:,0]]
    face_area_vectors=np.column_stack([edges[:,1],-edges[:,0]]);face_lengths=np.linalg.norm(face_area_vectors,axis=1);face_normals=face_area_vectors/face_lengths[:,None]
    require(np.allclose(face_lengths,DX,atol=1e-15,rtol=0),'Every source face edge retains full.0125m length')
    incidence=np.bincount(face_owner,minlength=9472)+np.bincount(face_neighbour[face_neighbour>=0],minlength=9472)
    require(np.all(incidence==4),'Exactly four incident faces on every cell')
    closure=np.zeros((9472,2));np.add.at(closure,face_owner,face_area_vectors);m=face_neighbour>=0;np.add.at(closure,face_neighbour[m],-face_area_vectors[m])
    require(float(abs(closure).max())<1e-15,'Each polygon closes its oriented face-area vectors')
    boundary=np.flatnonzero(face_neighbour<0);x_edges=np.arange(NX+1)*DX;y_edges=Y0+np.arange(NY+1)*DY
    out.mkdir(parents=True);mesh=out/'tandem-offset-wake.msh';ranges={};first=1
    with mesh.open('x',encoding='ascii',newline='\n') as f:
        f.write('(0 "Offset tandem-square source mesh;9472uniform quads;two genuine holes;no computed flow")\n(2 2)\n')
        f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
        for x,y in nodes:f.write(f'{x:.17g} {y:.17g}\n')
        f.write(f'))\n(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
        for z in ZONES:
            rows=groups[z];bc={3:2,4:10,5:5,6:7,7:7,8:3,9:3}[z]
            ranges[ZONE_NAMES[z]]=[first,first+len(rows)-1]
            f.write(f'(13 ({z:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
            for a,b,o,n in rows:f.write(f'{a+1:x} {b+1:x} {o+1:x} {n+1:x}\n')
            f.write('))\n');first+=len(rows)
        for z,kind,name in ((2,'fluid','fluid'),(3,'interior','interior'),(4,'velocity-inlet','inlet'),(5,'pressure-outlet','outlet'),(6,'symmetry','top'),(7,'symmetry','bottom'),(8,'wall','body1'),(9,'wall','body2')):f.write(f'(45 ({z} {kind} {name})())\n')
    np.savez_compressed(out/'geometry.npz',nodes=nodes,quads=quads,centres=centres,areas=areas,grid_shape=np.array([NY,NX]),ix=ix,iy=iy,x_edges=x_edges,y_edges=y_edges,
        grid_node_coordinates=grid_vertices,fluid_mask=fluid_mask,solid_mask=solid_mask,body1_mask=body1_mask,body2_mask=body2_mask,
        slot_to_source_cell=slot_to_source_cell,source_cell_to_slot=source_cell_to_slot,source_slot_coordinates=slots,
        face_nodes=face_nodes,face_owner=face_owner,face_neighbour=face_neighbour,face_zone=face_zone,face_centres=face_centres,face_normals=face_normals,
        face_area_vectors=face_area_vectors,face_lengths=face_lengths,face_source_ids=np.arange(1,19197,dtype=np.int64),
        boundary_face_source_ids=boundary+1,boundary_face_owner=face_owner[boundary],boundary_face_zone=face_zone[boundary],boundary_face_nodes=face_nodes[boundary],
        boundary_face_centres=face_centres[boundary],boundary_face_normals=face_normals[boundary],boundary_face_area_vectors=face_area_vectors[boundary],boundary_face_lengths=face_lengths[boundary])
    proof=dict(cells=9472,nodes=9723,faces=19196,interior_faces=18692,physical_boundary_faces=504,
        boundary_face_counts={ZONE_NAMES[z]:len(groups[z]) for z in ZONES[1:]},fluent_face_ranges_1based=ranges,
        positive_areas=True,minimum_area_m2=float(areas.min()),maximum_area_m2=float(areas.max()),total_area_m2=float(areas.sum()),
        minimum_edge_length_m=float(face_lengths.min()),maximum_edge_length_m=float(face_lengths.max()),
        all_nodes_used=True,duplicate_nodes=False,all_shared_edges_oppositely_oriented=True,all_cells_four_faces=True,
        connected_fluid_components=1,boundary_closed_loops=3,boundary_loop_nodes=[440,32,32],holes=2,euler_characteristic=-1,
        root_grid_slots=9600,removed_solid_slots=128,body1_slots=64,body2_slots=64,source_slot_bijection=True,
        orthogonality_proxy_min=float(min(orth)),maximum_nonorthogonality_proxy_deg=float(np.degrees(np.arccos(min(orth)))),maximum_cell_area_vector_closure_m=float(abs(closure).max()),
        native_solver_mesh_quality_pending=True,expected_openfoam_extrusion=dict(points=19446,faces=38140,internal_faces=18692,physical_faces=504,empty_faces=18944,depth_m=.01))
    setup=dict(model='Transient offset tandem-square wake interaction, generic incompressible laminar flow',cells=9472,nx=NX,ny=NY,dx_m=DX,dy_m=DY,
        density_kg_m3=1.,dynamic_viscosity_Pa_s=.0001,kinematic_viscosity_m2_s=.0001,inlet_velocity_m_s=[.1,0.],outlet_gauge_pressure_Pa=0.,
        Re_body=100.,Re_D=100.,body_side_m=.1,body_centres_m=[[.5,-.025],[.9,.025]],body1_centre_m=[.5,-.025],body2_centre_m=[.9,.025],
        body1_x_bounds_m=[.45,.55],body1_y_bounds_m=[-.075,.025],body2_x_bounds_m=[.85,.95],body2_y_bounds_m=[-.025,.075],
        body_slot_bounds_ij=[[36,44,24,32],[68,76,28,36]],domain_x_bounds_m=[0.,2.],domain_y_bounds_m=[-.375,.375],length_m=2.,height_m=.75,
        streamwise_clear_gap_m=.3,transverse_centre_offset_m=.05,reference_length_m=.1,reference_velocity_m_s=.1,reference_area_m2=.1,reference_depth_m=1.,openfoam_extrusion_depth_m=.01,
        openfoam_z_bounds_m=[-.005,.005],operating_pressure_Pa=0.,initial_velocity_m_s=[0.,0.],initial_pressure_Pa=0.,initial_state_is_prescribed_not_conservative_solved_flow=True,
        expected_flow_per_unit_depth_m2_s=dict(inlet=.075,outlet=.075),boundary_names=['inlet','outlet','top','bottom','body1','body2'],
        boundary_types=dict(inlet='velocity-inlet',outlet='pressure-outlet',top='symmetry',bottom='symmetry',body1='wall',body2='wall'),
        boundary_zone_ids={ZONE_NAMES[z]:z for z in ZONES},boundary_face_counts=proof['boundary_face_counts'],body_wall_motion='Stationary Wall',body_shear_condition='No Slip',
        time_step_s=.05,pilot_time_steps=5,pilot_final_time_s=.25,time_steps=800,final_time_s=40.,saved_state_count=801,animation_step_indices=list(range(0,801,8)),
        initial_physical_time_s=0.,initial_physical_step_index=0,energy=False,turbulence=False,gravity=False,MRF=False,moving_mesh=False,
        source_cell_order='j then i over fluid slots; source nodes sorted j then i; no hole-interior unused nodes',
        face_order='interior,inlet,outlet,top,bottom,body1,body2',geometry_index_convention='NPZ node/cell/owner indices0-based; face_neighbour=-1atboundary; *_source_ids Fluent1-based faceids',
        geometry_checks=proof,mesh_sha256=sha(mesh),geometry_npz_sha256=sha(out/'geometry.npz'),mesher_source_sha256=sha(source),target_document_sha256=TARGET_SHA,
        generated_by='Codex; mesh only, no computed flow',physical_validation=False,mesh_independence=False,time_step_independence=False,transient_shedding_established=False)
    json_new(out/'mesh-proof.json',proof);json_new(out/'setup.json',setup)
    (out/'executed-mesher.py').write_bytes(source.read_bytes());require(sha(source)==sha(out/'executed-mesher.py'),'Executed source archive exact')
    json_new(out/'source-manifest.json',dict(kind='Source mesh only; native solver quality pending',source_path=str(source),source_sha256=sha(source),target_document_sha256=TARGET_SHA,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}))
    print(json.dumps(dict(out=str(out),mesh_sha256=sha(mesh),geometry_npz_sha256=sha(out/'geometry.npz'),setup_sha256=sha(out/'setup.json'),proof_sha256=sha(out/'mesh-proof.json'),source_sha256=sha(source),geometry_checks=proof),indent=2))


if __name__=='__main__':main()