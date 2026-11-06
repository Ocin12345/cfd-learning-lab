"""Bounded interior-node repair of the preserved v002 wing mesh; no remeshing/CFD.

Only free interior nodes in stars of failing tetrahedra may move. All boundary
coordinates, cells, profile and CAD stay unchanged. One fresh candidate is saved
with a complete witness even on failure; no target reduction or automatic retry.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import shutil
import time
import numpy as np

FLOOR=.105


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def digest(values):
    return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def boundary_faces(path):
    import re
    faces={};capture=False;zone=None
    for line in path.read_text(encoding='ascii').splitlines():
        match=re.match(r'^\(13 \(([0-9a-f]+) [0-9a-f]+ [0-9a-f]+ [0-9a-f]+ 3\)\($',line)
        if match:zone=int(match[1],16);capture=True;continue
        if capture and line=='))':capture=False;continue
        if capture and zone!=3:
            nodes=[int(value,16)-1 for value in line.split()[:3]]
            key=tuple(sorted(nodes))
            if key in faces:raise ValueError('Duplicate source exterior face')
            faces[key]=zone
    if not faces or set(faces.values())!={4,5,6,7}:raise ValueError('Wrong source patch contract')
    return faces


def update_gmsh(source,destination,old_points,new_points,moved):
    """Map source node tags by a verified coordinate bijection, never by ordering."""
    lines=source.read_text(encoding='ascii').splitlines(keepends=True)
    start=next(i for i,line in enumerate(lines) if line.strip()=='$Nodes');count=int(lines[start+1])
    if count!=len(old_points):raise ValueError('Gmsh/NPZ node count mismatch')
    lookup={}
    for i,row in enumerate(old_points):
        key=tuple(np.round(row,11))
        if key in lookup:raise ValueError('Source coordinate key is not unique')
        lookup[key]=i
    seen=set();changed=[]
    for offset in range(count):
        index=start+2+offset;fields=lines[index].split();tag=int(fields[0]);xyz=np.asarray(fields[1:],float)
        node=lookup.get(tuple(np.round(xyz,11)))
        if node is None or node in seen or not np.allclose(xyz,old_points[node],rtol=0,atol=1e-10):
            raise ValueError('Source Gmsh node does not bijectively match NPZ coordinates')
        seen.add(node)
        if node in moved:
            lines[index]=str(tag)+' '+' '.join(f'{v:.16g}' for v in new_points[node])+'\n'
            changed.append(dict(gmsh_tag=tag,npz_node_zero_based=int(node)))
    if len(seen)!=len(old_points) or len(changed)!=len(moved):raise ValueError('Incomplete source Gmsh bijection')
    destination.write_text(''.join(lines),encoding='ascii',newline='')
    return changed


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();source=a.source.resolve();out=a.out.resolve()
    if out.exists() or out.drive.upper()!='E:' or not source.is_dir():raise ValueError('Use preserved source and fresh E: output')
    if shutil.disk_usage(out.anchor).free<10*1024**3:raise RuntimeError('Require10GiBfreeE:')
    helper_path=Path(__file__).with_name('wing_mesh_v002.py')
    spec=importlib.util.spec_from_file_location('wing_helper_v002',helper_path);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    source_paths=[source/name for name in ('geometry.npz','setup.json','wing.msh','source-gmsh.msh','profile.json','geometry.brep','wing-surface.stl','candidate-native-quality.json','executed-mesher.py')]
    before={str(path):sha256(path) for path in source_paths}
    g=np.load(source/'geometry.npz');arrays={key:g[key].copy() for key in g.files};g.close()
    old_points=arrays['points_xyz_m'].copy();points=old_points.copy();cells=arrays['tets'].copy()
    boundary=boundary_faces(source/'wing.msh');fixed=np.asarray(sorted({n for tri in boundary for n in tri}),dtype=np.int64)
    fixed_set=set(map(int,fixed));initial_q,_=helper.native_tet_volume_quality(points,cells);quality=initial_q.copy()
    stars=[[] for _ in range(len(points))]
    for cid,tet in enumerate(cells):
        for node in tet:stars[node].append(cid)
    directions=np.asarray([v for v in itertools.product((-1.,0.,1.),repeat=3) if any(v)])
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    attempts=0;accepted=[];moved=set();start=time.monotonic();reason=''
    bad_initial=np.where(quality<FLOOR)[0].tolist()
    for sweep in range(12):
        if float(quality.min())>=FLOOR+1e-7:reason='global floor achieved';break
        if time.monotonic()-start>115:reason='bounded optimization time exhausted';break
        worst=int(np.argmin(quality));best=None
        for node in cells[worst]:
            node=int(node)
            if node in fixed_set:continue
            star=np.asarray(stars[node],dtype=np.int64);old=points[node].copy()
            neighbors=np.unique(cells[star]);neighbors=neighbors[neighbors!=node]
            h=float(np.median(np.linalg.norm(old_points[neighbors]-old_points[node],axis=1)))
            centre=points[neighbors].mean(axis=0)-old
            local_directions=[v for v in directions]
            if np.linalg.norm(centre)>0:local_directions.append(centre/np.linalg.norm(centre))
            protected=quality[star]>=FLOOR
            current_min=float(quality[star].min())
            for fraction in (.015,.03,.06,.10,.15):
                for direction in local_directions:
                    if time.monotonic()-start>115:break
                    candidate=old+h*fraction*direction;attempts+=1
                    if np.linalg.norm(candidate-old_points[node])>.20*h:continue
                    points[node]=candidate
                    tetra=points[cells[star]];edges=tetra[:,1:]-tetra[:,:1]
                    signed=np.einsum('ij,ij->i',np.cross(edges[:,0],edges[:,1]),edges[:,2])/6.
                    if np.any(signed<=0):points[node]=old;continue
                    try:q,_=helper.native_tet_volume_quality(points,cells[star])
                    except (ValueError,np.linalg.LinAlgError):points[node]=old;continue
                    points[node]=old
                    if np.any(q[protected]<FLOOR):continue
                    score=float(q.min())
                    if score<=current_min+1e-7:continue
                    if best is None or score>best['score']:
                        best=dict(node=node,star=star,q=q,point=candidate.copy(),score=score,
                                  before_star_minimum=current_min,h=h,old=old.copy())
            points[node]=old
        if best is None:reason='no small feasible local improvement';break
        node=best['node'];points[node]=best['point'];quality[best['star']]=best['q'];moved.add(node)
        accepted.append(dict(sweep=sweep,node_zero_based=node,old=best['old'].tolist(),new=points[node].tolist(),
            node_star_tets=len(best['star']),before_star_minimum=best['before_star_minimum'],after_star_minimum=best['score'],
            displacement_from_source_m=float(np.linalg.norm(points[node]-old_points[node])),original_neighbor_median_edge_m=best['h']))
    elapsed=time.monotonic()-start
    q,radius=helper.native_tet_volume_quality(points,cells)
    signed_cells,volumes,flipped=helper.positive_tets(points,cells)
    if flipped or not np.array_equal(signed_cells,cells):raise ValueError('Interior relocation changed tetrahedron orientation')
    if not np.array_equal(points[fixed],old_points[fixed]):raise ValueError('Boundary coordinate changed')
    groups,centres,proxy=helper.face_connectivity(points,cells,boundary)
    passed=bool(float(q.min())>=FLOOR and proxy>=FLOOR and np.all(q[initial_q>=FLOOR]>=FLOOR))
    out.mkdir(parents=True)
    for name in ('profile.json','geometry.brep','wing-surface.stl'):shutil.copyfile(source/name,out/name)
    (out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
    shutil.copyfile(source/'executed-mesher.py',out/'executed-parent-mesher.py')
    shutil.copyfile(helper_path,out/'executed-mesh-helper.py')
    helper.write_fluent(out/'wing.msh',points,cells,groups)
    tag_map=update_gmsh(source/'source-gmsh.msh',out/'source-gmsh.msh',old_points,points,moved)
    arrays.update(points_xyz_m=points,nodes=points,tet_centres=centres,tet_volumes=volumes,
        tet_native_volume_quality=q,tet_circumradius=radius)
    np.savez_compressed(out/'geometry.npz',**arrays)
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.NumThreads',2);gmsh.open(str(out/'source-gmsh.msh'))
        kinds,tags,_=gmsh.model.mesh.getElements(3)
        if list(kinds)!=[4] or len(tags[0])!=len(cells):raise ValueError('Exported Gmsh cell count/type mismatch')
        gmsh_quality={name:np.asarray(gmsh.model.mesh.getElementQualities(tags[0],name)) for name in ('gamma','minSICN')}
    finally:gmsh.finalize()
    witness=dict(scope='Interior-node relocation only; no remeshing or boundary/CAD/profile change',
        initial_failing_cells=bad_initial,initial_minimum=float(initial_q.min()),final_minimum=float(q.min()),
        final_below_native_0p10=int((q<.10).sum()),final_below_margin_0p105=int((q<FLOOR).sum()),
        required_floor=FLOOR,face_centre_orthogonality_proxy_minimum=proxy,candidate_gate_pass=passed,
        optimization_seconds=elapsed,candidate_offsets_tested=attempts,stop_reason=reason,moved_node_count=len(moved),accepted_moves=accepted,
        moved_nodes=[dict(node_zero_based=node,before=old_points[node].tolist(),after=points[node].tolist()) for node in sorted(moved)],
        fixed_boundary_node_count=len(fixed),fixed_boundary_node_ids_sha256=digest(fixed),
        boundary_coordinates_before_sha256=digest(old_points[fixed]),boundary_coordinates_after_sha256=digest(points[fixed]),
        boundary_coordinates_byte_identical=np.array_equal(points[fixed],old_points[fixed]),
        connectivity_before_sha256=digest(cells),connectivity_after_sha256=digest(arrays['tets']),
        cell_count_unchanged=len(cells),node_count_unchanged=len(points),positive_volumes=True,
        formerly_passing_tets_remain_above_floor=bool(np.all(q[initial_q>=FLOOR]>=FLOOR)),
        source_gmsh_coordinate_bijection_verified=True,changed_gmsh_node_tags=tag_map,
        profile_before_sha256=sha256(source/'profile.json'),profile_after_sha256=sha256(out/'profile.json'),
        cad_before_sha256=sha256(source/'geometry.brep'),cad_after_sha256=sha256(out/'geometry.brep'),
        wing_surface_before_sha256=sha256(source/'wing-surface.stl'),wing_surface_after_sha256=sha256(out/'wing-surface.stl'),
        parent_source_hashes=before,relocator_source_sha256=sha256(Path(__file__)),helper_source_sha256=sha256(helper_path),
        native_fluent_check_still_required=True,CFD_launched=False)
    write_json(out/'node-relocation-witness.json',witness)
    receipt=dict(definition='V / [(8/(9sqrt3))*R_circum^3], one minus equilateral volume skewness',
        minimum=float(q.min()),percentile1=float(np.percentile(q,1)),median=float(np.median(q)),
        count_below_native_0p10=int((q<.10).sum()),count_below_candidate_0p105=int((q<FLOOR).sum()),
        face_centre_orthogonality_proxy_minimum=proxy,required_candidate_floor=FLOOR,
        all_candidate_tets_checked=len(cells),candidate_gate_pass=passed,native_fluent_quality_check_still_required=True,
        geometry_changed=False,repair_scope='Interior relocation of only failing tetrahedron-star free nodes',witness='node-relocation-witness.json')
    write_json(out/'candidate-native-quality.json',receipt)
    setup=json.loads((source/'setup.json').read_text(encoding='utf-8'))
    setup.update(mesh_version='v003',same_geometry_as_v001=True,candidate_quality_gate=receipt,
        minimum_tet_volume_m3=float(volumes.min()),maximum_tet_volume_m3=float(volumes.max()),
        mesh_fluid_volume_m3=float(volumes.sum()),mesh_solid_volume_m3=1152.-float(volumes.sum()),
        face_centre_orthogonality_proxy_minimum=proxy,
        gmsh_quality={name:dict(minimum=float(v.min()),percentile1=float(np.percentile(v,1)),median=float(np.median(v))) for name,v in gmsh_quality.items()},
        mesh_sha256=sha256(out/'wing.msh'),gmsh_mesh_sha256=sha256(out/'source-gmsh.msh'),
        geometry_npz_sha256=sha256(out/'geometry.npz'),mesher_source_sha256=sha256(Path(__file__)),
        parent_mesh_sha256=sha256(source/'wing.msh'),interior_relocation_witness='node-relocation-witness.json')
    write_json(out/'setup.json',setup)
    unchanged=all(sha256(path)==value for path,value in before.items())
    if not unchanged:raise ValueError('Preserved v002 input changed during relocation')
    write_json(out/'mesh-source-manifest.json',dict(parent_inputs_unchanged=True,
        files=[dict(name=path.name,sha256=sha256(path),bytes=path.stat().st_size) for path in sorted(out.iterdir()) if path.is_file()]))
    print(json.dumps(dict(candidate_gate_pass=passed,quality=receipt,witness=witness,setup=setup),indent=2),flush=True)
    return 0 if passed else 2


if __name__=='__main__':raise SystemExit(main())
