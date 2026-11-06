"""Read-only accepted-flow geometry adapter. Never reads guessed/new UDS data or runs CFD."""
from pathlib import Path
import ast,hashlib,importlib.util,json,sys
import h5py
import numpy as np
ROOT=Path('E:/CFD-Learning-Lab');PRIVATE=ROOT/'local/prospective-tracer-v001';QUAL=PRIVATE/'native-flow-qualification-v001.json';QUAL_SHA='dccbb17e813f2552298332dc3fd7ff713e30df71a70508316a30831995ccd312'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def require(v,m):
 if not v:raise ValueError(m)
def qualified_geometry():
 require(sha(QUAL)==QUAL_SHA,'Pinned prospective qualification');q=json.loads(QUAL.read_text(encoding='utf-8'));pins={Path(r['path']).resolve():r['sha256_after'] for r in q['source_manifest']};before={QUAL:QUAL_SHA,Path(__file__).resolve():sha(__file__)}
 def touch(p):
  p=Path(p).resolve();require(p in pins and sha(p)==pins[p],'Qualification input hash '+str(p));before[p]=pins[p];return p
 def module(name,p):
  p=touch(p);require(name not in sys.modules,'Fresh checked import namespace '+name);s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
 # Authenticate every imported dependency BEFORE executing any imported source.
 helper=ROOT/'repo/cfd-learning-lab/scripts/analyse_nozzle_fluent_v005.py';fp=ROOT/'local/tjunction-v001/fluent-checks-v001/executed-analyser.py';op=ROOT/'local/tjunction-v001/openfoam-checks-v001/executed-analyse_tjunction_openfoam_v001.py';prep=ROOT/'repo/cfd-learning-lab/scripts/prepare_openfoam_tjunction_v001.py'
 for p in (helper,fp,op,prep):touch(p)
 fhelper=module('analyse_nozzle_fluent_v005',helper);f=module('qualified_original_tjunction_fluent',fp);module('prepare_openfoam_tjunction_v001',prep);o=module('qualified_original_tjunction_openfoam',op)
 # Local source is explicit code dependency, not a native dataset selector.
 local=PRIVATE/'scalar_fv_audit_v002.py';before[local]=sha(local);spec=importlib.util.spec_from_file_location('qualified_scalar_fv',local);fv=importlib.util.module_from_spec(spec);spec.loader.exec_module(fv)
 N=ROOT/'native/runs/2026-11-01-tjunction-v001';gp=touch(N/'mesh/geometry.npz')
 with np.load(gp,allow_pickle=False) as z:g={k:z[k] for k in ('nodes','quads','centres','areas')}
 edges={}
 for cid,quad in enumerate(g['quads']):
  for a,b in zip(quad,np.roll(quad,-1)):edges.setdefault(tuple(sorted((int(a),int(b)))),[]).append(cid)
 require(len(edges)==9040 and all(len(v) in (1,2) for v in edges.values()),'Actual canonical 9040-edge incidence')
 lookup={tuple(np.round(x,10)):i for i,x in enumerate(g['nodes'])}
 def point_map(points):
  out=np.array([lookup[tuple(np.round(p[:2],10))] for p in points]);require(np.max(abs(g['nodes'][out]-np.asarray(points)[:,:2]))<1e-10,'Actual native/source vertices');return out
 def edge_check(stencil,pairs,mapping):
  seen={};normal={}
  for i,pair in enumerate(pairs):
   k=tuple(sorted(map(int,pair)));require(k in edges and k not in seen,'Native active edge unique/source identity');adj=[int(mapping[stencil['owner'][i]])]
   if stencil['neighbour'][i]>=0:adj.append(int(mapping[stencil['neighbour'][i]]))
   require(sorted(adj)==sorted(edges[k]),'Native actual owner/neighbour/source incidence');seen[k]=i;normal[k]=stencil['area_vectors'][i]/stencil['depth_m']*(1 if adj[0]==min(adj) else -1)
  require(set(seen)==set(edges),'All 9040 source edges covered');return seen,normal
 case=touch(N/'fluent/tjunction-final.cas.h5');data=touch(N/'fluent/tjunction-final.dat.h5');c0,c1,phi,zones=f.native_faces(case,data);mapping=f.common_mapping(case,g);require(len(phi)==9040,'Actual Fluent9040 faces')
 with h5py.File(case,'r') as h:
  r=h['meshes/1'];points=fhelper.chunk_array(r['nodes/coords'],len(g['nodes']),(2,));pairs=np.empty((9040,2),int);seen=np.zeros(9040,bool)
  for block in r['faces/nodes'].values():
   lo=int(np.asarray(block.attrs['minId']).ravel()[0])-1;hi=int(np.asarray(block.attrs['maxId']).ravel()[0]);require(not seen[lo:hi].any() and np.all(np.asarray(block['nnodes'][:])==2),'Actual face-node chunks');pairs[lo:hi]=np.asarray(block['nodes'][:]).reshape(-1,2).astype(int)-1;seen[lo:hi]=True
  require(seen.all(),'Full Fluent face coverage')
 vertices=[set() for _ in range(4400)]
 for i,pair in enumerate(pairs):
  for cid in (c0[i],c1[i]):
   if cid:vertices[cid-1].update(map(int,pair))
 centres=np.empty((4400,2));area=np.empty(4400)
 for i,ids in enumerate(vertices):
  require(len(ids)==4,'Native cell quad vertices');xy=points[sorted(ids)];mid=xy.mean(axis=0);xy=xy[np.argsort(np.arctan2(xy[:,1]-mid[1],xy[:,0]-mid[0]))];nxt=np.roll(xy,-1,axis=0);cross=xy[:,0]*nxt[:,1]-xy[:,1]*nxt[:,0];a=.5*cross.sum();require(a>0,'Actual positive native polygon area');area[i]=a;centres[i]=np.sum((xy+nxt)*cross[:,None],axis=0)/(6*a)
 require(np.max(abs(area-g['areas'][mapping]))<1e-14 and np.max(abs(centres-g['centres'][mapping]))<1e-12,'Fluent actual polygon area/centroid source identity');fc=points[pairs].mean(axis=1);edge=points[pairs[:,1]]-points[pairs[:,0]];sf=np.column_stack((-edge[:,1],edge[:,0]));direction=fc-centres[c0-1];inside=c1>0;direction[inside]=centres[c1[inside]-1]-centres[c0[inside]-1];sf[np.sum(sf*direction,axis=1)<0]*=-1;patch=np.full(9040,'',dtype='U32')
 for name,z in zones.items():patch[z['first']:z['last']]=name
 F=fv.make_stencil(owner=c0-1,neighbour=c1-1,cell_centres=centres,cell_volumes=area,face_centres=fc,area_vectors=sf,patch=patch,depth_m=1.,native_face_ids=np.arange(1,9041));fe,fn=edge_check(F,point_map(points)[pairs],mapping)
 with h5py.File(data,'r') as h:rho=fhelper.chunk_array(h['results/1/phase-1/cells/SV_DENSITY'],4400)
 run=N/'openfoam';poly=o.native_poly(run)
 for name in ('points','faces','owner','neighbour','boundary'):touch(poly/name)
 geom=o.native_geometry(poly);pa=o.native_patches(poly);C=o.field(touch(run/'2000/C').read_text(),4400,3,[0,1,0,0,0,0,0]);V=np.asarray(geom['volumes']);points3=np.asarray(geom['points']);faces=np.asarray(geom['faces']);xyz=points3[faces];Sf=.5*np.sum(np.cross(xyz,np.roll(xyz,-1,axis=1)),axis=1);Cf=xyz.mean(axis=1);active=np.flatnonzero(np.linalg.norm(Sf[:,:2],axis=1)>1e-12);require(len(faces)==17840 and len(active)==9040,'Actual OF17840 faces/9040 vertical');require(np.max(abs(Sf[active,2]))<1e-12,'Planar vertical area vectors')
 raw=touch(run/'2000/phi').read_text();allphi=np.zeros(17840);allphi[:8560]=o.field(raw,8560,1,[0,3,-1,0,0,0,0]);blocks=o.boundary_blocks(raw);names=np.full(17840,'interior',dtype='U32')
 for name in o.PATCHES:
  z=pa[name];ids=slice(z['start'],z['start']+z['count']);allphi[ids]=o.parse_value(blocks[name],z['count'],1,'value');names[ids]=name
 z=pa[o.EMPTY];names[z['start']:z['start']+z['count']]='excluded_empty';require(z['count']==8800 and not np.any(names[active]=='excluded_empty'),'Explicit8800 empty planes excluded; not invented observed scalar flux')
 neigh=np.full(17840,-1);neigh[:8560]=geom['neighbour'];omap=o.cell_mapping(C,gp,V,native_points=points3)[0]['native_cell_to_source_cell'];omap=np.asarray(omap);O=fv.make_stencil(owner=np.asarray(geom['owner'])[active],neighbour=neigh[active],cell_centres=C[:,:2],cell_volumes=V,face_centres=Cf[active,:2],area_vectors=Sf[active,:2],patch=names[active],depth_m=.01,native_face_ids=active)
 pm=point_map(points3);opairs=np.array([list(set(pm[face])) for face in faces[active]]);oe,on=edge_check(O,opairs,omap);require(max(np.max(abs(fn[k]-on[k])) for k in edges)<1e-12,'Independent native F/OF matched oriented per-depth edge vectors')
 require(all(sha(p)==digest for p,digest in before.items()),'Consumed old inputs/source unchanged')
 report=dict(source_preparation_only=True,scalar_fields_read=False,scalar_case_created=False,CFD_started=False,scalar_acceptance_claimed=False,native_cells_each=4400,fluent_faces=9040,openfoam_faces=17840,openfoam_active_faces=9040,openfoam_excluded_empty_faces=8800,matched_source_edge_count=9040,actual_domain_area_m2=float(area.sum()),actual_orthogonal_distances_m=dict(fluent_min=float(F['distance'].min()),fluent_max=float(F['distance'].max()),openfoam_min=float(O['distance'].min()),openfoam_max=float(O['distance'].max())),max_Co_dt0p025=dict(fluent=float((.025*.5*(np.bincount(F['owner'],weights=abs(phi),minlength=4400)+np.bincount(F['neighbour'][F['inside']],weights=abs(phi[F['inside']]),minlength=4400))/(rho*F['volumes'])).max()),openfoam=float((.025*.5*(np.bincount(O['owner'],weights=abs(allphi[active]),minlength=4400)+np.bincount(O['neighbour'][O['inside']],weights=abs(allphi[active][O['inside']]),minlength=4400))/V).max())),consumed_inputs=[dict(path=str(p),sha256=digest) for p,digest in before.items()],all_consumed_inputs_unchanged=True)
 return dict(fluent=dict(stencil=F,native_mass_phi=phi,rho=rho,native_cell_to_source_cell=mapping),openfoam=dict(stencil=O,native_mass_phi=allphi[active],rho=np.ones(4400),native_cell_to_source_cell=omap)),report
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args();require(a.out.parent.resolve()==PRIVATE.resolve() and not a.out.exists(),'Fresh private geometry-only receipt');_,report=qualified_geometry()
 with a.out.open('x',encoding='utf-8') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps(dict(report=str(a.out),sha256=sha(a.out),native_cells_each=4400,matched_active_edges=9040,scalar_fields_read=False)))
