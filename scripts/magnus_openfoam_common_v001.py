"""Parameterized native ASCII mesh/rank/force readers. Imports require no private archives or native runtime."""
from pathlib import Path
from dataclasses import dataclass
import hashlib,json,re,math,shutil
FLOAT=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
@dataclass(frozen=True)
class MeshSpec:
 cells:int
 points:int
 faces:int
 internal:int
 patch_counts:dict
 depth:float=.01
 empty:str='frontAndBackPlanes'

def require(v,m):
 if not v:raise ValueError(m)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def clean(s):return re.sub(r'/\*.*?\*/|//[^\n]*','',s,flags=re.S)
def write_new(p,s):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x',encoding='ascii',newline='\n') as f:f.write(s.rstrip()+'\n')
def json_new(p,v):write_new(p,json.dumps(v,indent=2,allow_nan=False))
def copied(a,b):
 a=Path(a);b=Path(b);require(not b.exists(),'Preserve '+str(b));b.parent.mkdir(parents=True,exist_ok=True);h=sha(a);shutil.copy2(a,b);require(sha(a)==sha(b)==h,'Copy hash differs')
def header(n,c='dictionary'):return 'FoamFile\n{\n version 2.0;\n format ascii;\n class '+c+';\n object '+n+';\n}\n\n'
def counted(p):
 t=clean(Path(p).read_text());h=re.search(r'FoamFile\s*\{([^{}]*)\}',t);require(h is not None and re.search(r'\bformat\s+ascii\s*;',h[1]),'Native ASCII list');m=re.fullmatch(r'\s*(\d+)\s*\((.*)\)\s*',t[h.end():],re.S);require(m is not None,'Counted native list '+str(p));return int(m[1]),m[2]
def values(t,n,c,label='internalField'):
 m=re.search(r'\b'+label+r'\s+uniform\s+([^;]+);',t)
 if m:
  v=list(map(float,m[1].replace('(',' ').replace(')',' ').split()));require(len(v)==c,'Uniform components');return [v[:] for _ in range(n)] if c>1 else [v[0]]*n
 kind='vector' if c==3 else 'scalar';m=re.search(r'\b'+label+r'\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',t,re.S);require(m is not None and int(m[1])==n,'Native list count');v=list(map(float,m[2].replace('(',' ').replace(')',' ').split()));require(len(v)==n*c,'Native component length');return [v[i:i+c] for i in range(0,len(v),c)] if c>1 else v
def field(raw,n,c,dims):
 t=clean(raw);m=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',t);require(m is not None and list(map(int,m[1].split()))==dims and re.search(r'\bformat\s+ascii\s*;',t),'Actual field dimensions/ASCII');return values(t,n,c)
def blocks(raw):
 t=clean(raw);require('boundaryField' in t,'Boundary field missing');i=t.index('{',t.index('boundaryField'))+1;out={}
 while i<len(t):
  m=re.match(r'\s*("[^"]+"|[^\s{}]+)\s*\{',t[i:])
  if not m:break
  name=m[1].strip('"');a=i+m.end();j=a;level=1
  while level:
   require(j<len(t),'Unclosed patch dictionary');level+=(t[j]=='{')-(t[j]=='}');j+=1
  require(name not in out,'Duplicate patch');out[name]=t[a:j-1];i=j
 return out
def top(t):
 out=[];level=0
 for c in t:
  if c=='{':level+=1
  elif c=='}':level-=1
  elif level==0:out.append(c)
 return ''.join(out)
def kind(t):
 m=re.search(r'\btype\s+(\w+)\s*;',top(t));require(m is not None,'Actual patch field type');return m[1]
def patch_table(poly):
 n,t=counted(poly/'boundary');rows={name.strip('"'):dict(kind=re.search(r'\btype\s+(\w+)\s*;',b)[1],count=int(re.search(r'\bnFaces\s+(\d+)\s*;',b)[1]),start=int(re.search(r'\bstartFace\s+(\d+)\s*;',b)[1])) for name,b in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',t)};require(len(rows)==n,'Native mesh patch count');return rows
def native_poly(case):
 p=[case/n for n in ('0/polyMesh','constant/polyMesh') if (case/n/'points').is_file()];require(len(p)==1,'Exactly one initial global native mesh');return p[0]
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def native_geometry(poly,spec):
 n,b=counted(poly/'points');points=[list(map(float,s.split())) for s in re.findall(r'\(([^()]*)\)',b)];require(n==spec.points==len(points),'Native point count');n,b=counted(poly/'faces');faces=[list(map(int,s.split())) for k,s in re.findall(r'(\d+)\s*\(([^()]*)\)',b)];require(n==spec.faces==len(faces) and all(len(f)==4 and len(set(f))==4 for f in faces),'Native quad faces')
 def labels(name,count):
  n,b=counted(poly/name);v=list(map(int,b.split()));require(n==count==len(v) and all(0<=i<spec.cells for i in v),'Native owner/neighbour');return v
 owner=labels('owner',spec.faces);neighbour=labels('neighbour',spec.internal);V=[0.]*spec.cells;Sf=[];Cf=[];vertices=[set() for _ in range(spec.cells)]
 for i,ids in enumerate(faces):
  v=[points[j] for j in ids];sf=[0.,0.,0.];volume=0.
  for a,b in zip(v,v[1:]+v[:1]):sf=[x+.5*y for x,y in zip(sf,cross(a,b))]
  for k in (1,2):volume+=dot(v[0],cross(v[k],v[k+1]))/6
  V[owner[i]]+=volume;vertices[owner[i]].update(ids)
  if i<spec.internal:V[neighbour[i]]-=volume;vertices[neighbour[i]].update(ids)
  Sf.append(sf);Cf.append([sum(q[j] for q in v)/4 for j in range(3)])
 require(all(math.isfinite(v) and v>0 for v in V) and all(len(v)==8 for v in vertices),'Positive oriented native hexa volumes');return dict(points=points,faces=faces,owner=owner,neighbour=neighbour,volumes=V,Sf=Sf,Cf=Cf,cell_vertices=[sorted(v) for v in vertices])
def source_mapping(g,src,C,V,spec,patches,zones):
 lookup={tuple(round(x,10) for x in p):i for i,p in enumerate(src['nodes'])};pm=[];plane=[]
 for x,y,z in g['points']:
  key=(round(x,10),round(y,10));require(key in lookup and min(abs(z-spec.depth/2),abs(z+spec.depth/2))<=1e-12,'Source/native projected vertex/extrusion');pm.append(lookup[key]);plane.append(1 if z>0 else -1)
 require(len(set(zip(pm,plane)))==spec.points==2*len(src['nodes']),'Both native planes/source nodes bijective');clook={frozenset(q):i for i,q in enumerate(src['quads'])};n2s=[]
 for vertices in g['cell_vertices']:
  key=frozenset(pm[j] for j in vertices);require(key in clook and len({pm[j] for j in vertices if plane[j]>0})==4 and len({pm[j] for j in vertices if plane[j]<0})==4,'Native8vertex/sourcequad');n2s.append(clook[key])
 require(sorted(n2s)==list(range(spec.cells)),'Native/source cell bijection');s2n=[0]*spec.cells
 for n,s in enumerate(n2s):s2n[s]=n
 require(all(math.dist(C[n][:2],src['centres'][s])<=1e-10 and abs(C[n][2])<=1e-12 and abs(V[n]-g['volumes'][n])<=1e-13 and abs(V[n]/spec.depth-src['areas'][s])<=1e-12 for n,s in enumerate(n2s)),'Actual native C/Vc/source areas')
 edge={frozenset(f):i for i,f in enumerate(src['face_nodes'])};fmap=[-1]*spec.faces;signs=[0]*spec.faces;seen=[]
 for f,ids in enumerate(g['faces']):
  nodes=frozenset(pm[j] for j in ids)
  if len(nodes)!=2:continue
  require(nodes in edge,'Native verticalface/sourceedge');s=edge[nodes];fmap[f]=s;expected=[x*spec.depth for x in src['face_area_vectors'][s]]+[0.];sign=1 if dot(g['Sf'][f],expected)>0 else -1;signs[f]=sign;require(math.dist(g['Sf'][f],[sign*x for x in expected])<=1e-12,'Native/source face area vector');o=src['face_owner'][s];n=src['face_neighbour'][s];require(n2s[g['owner'][f]]==(o if sign>0 else n),'Native/source face owner')
  if f<spec.internal:require(n2s[g['neighbour'][f]]==(n if sign>0 else o),'Native/source face neighbour')
  else:require(n==-1 and sign==1,'Native source physical wall orientation')
  seen.append(s)
 require(sorted(seen)==list(range(len(src['face_nodes']))),'Each canonical source edge once')
 for name,zone in zones.items():
  q=patches[name];require(q['count']==spec.patch_counts[name] and all(src['face_zone'][fmap[f]]==zone for f in range(q['start'],q['start']+q['count'])),'Actual physical patch/source zone identity '+name)
 return dict(native_cell_to_source_cell=n2s,source_native_cell_index=s2n,native_point_to_source_node=pm,native_point_z_sign=plane,native_face_to_source_edge=fmap,native_face_source_orientation=signs)
class Evidence:
 def __init__(self):self.before={}
 def touch(self,p):
  p=Path(p).resolve();require(p.is_file(),'Missing native evidence '+str(p))
  if p not in self.before:self.before[p]=sha(p)
  return p
 def text(self,p):return self.touch(p).read_text(encoding='utf-8-sig')
 def json(self,p):return json.loads(self.text(p))
 def finish(self):
  rows=[dict(path=str(p),bytes=p.stat().st_size,sha256_before=h,sha256_after=sha(p)) for p,h in self.before.items()];require(all(r['sha256_before']==r['sha256_after'] for r in rows),'Consumed input changed during read-only analysis');return dict(sources=rows,inputs_unchanged_during_analysis=True,solver_started_or_modified=False)
def finite_json(v):
 if isinstance(v,dict):return {k:finite_json(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [finite_json(x) for x in v]
 if isinstance(v,bool):return v
 if isinstance(v,float):return v if math.isfinite(v) else None
 return v
class Ranks:
 def __init__(self,e,case,g,patches,spec):
  import numpy as np
  self.e=e;self.case=Path(case);self.g=g;self.spec=spec;self.patches=patches;self.rows=[];self.mesh_paths=[];cells_seen=[];faces_seen=np.zeros(spec.faces,dtype=int);physical_seen=[];points=np.asarray(g['points']);faces=np.asarray(g['faces']);owner=np.asarray(g['owner']);neighbour=np.asarray(g['neighbour']);Sf=np.asarray(g['Sf'])
  for rank in (0,1):
   poly=self.case/f'processor{rank}/0/polyMesh';self.mesh_paths.append({n:str(poly/n) for n in ('points','faces','owner','neighbour','boundary','cellProcAddressing','faceProcAddressing','pointProcAddressing')})
   def labels(n):
    count,b=counted(e.touch(poly/n));a=np.array(list(map(int,b.split())),dtype=np.int64);require(len(a)==count,'Rank label count');return a
   cells=labels('cellProcAddressing');addr=labels('faceProcAddressing');own=labels('owner');nei=labels('neighbour');pm=labels('pointProcAddressing');cells_seen.extend(cells.tolist());require(np.all((cells>=0)&(cells<spec.cells)) and len(set(cells))==len(cells),'Rank cell addressing')
   n,b=counted(e.touch(poly/'points'));rp=np.array([list(map(float,s.split())) for s in re.findall(r'\(([^()]*)\)',b)]);require(rp.shape==(n,3) and len(pm)==n and np.all((pm>=0)&(pm<spec.points)) and np.max(abs(rp-points[pm]))<=1e-12,'Rank projected native point mapping');n,b=counted(e.touch(poly/'faces'));rf=np.array([list(map(int,s.split())) for size,s in re.findall(r'(\d+)\s*\(([^()]*)\)',b)]);require(rf.shape==(n,4) and len(own)==len(addr)==n and np.all(addr!=0) and np.all(abs(addr)<=spec.faces) and len(set(abs(addr)))==len(addr),'Rank face addressing')
   for lf,a in enumerate(addr):
    fi=abs(int(a))-1;sign=1 if a>0 else -1;require(set(pm[rf[lf]].tolist())==set(faces[fi].tolist()),'Every rank face vertex set');v=rp[rf[lf]];sf=.5*np.sum(np.cross(v,np.roll(v,-1,axis=0)),axis=0);require(np.max(abs(sf-sign*Sf[fi]))<=1e-13,'Rank signed face area orientation');require(cells[own[lf]]==(owner[fi] if sign>0 else neighbour[fi] if fi<spec.internal else -1),'Rank owner orientation')
    if lf<len(nei):require(cells[nei[lf]]==(neighbour[fi] if sign>0 else owner[fi]),'Rank neighbour orientation')
   faces_seen[abs(addr)-1]+=1;e.touch(poly/'boundary');bp=patch_table(poly);maps={}
   for name,qg in patches.items():
    q=bp[name];local=np.arange(q['start'],q['start']+q['count']);ids=abs(addr[local])-1;require(np.all(addr[local]>0) and set(ids).issubset(set(range(qg['start'],qg['start']+qg['count']))),'Rank physical/empty patch membership');physical_seen.extend(ids.tolist());maps[name]=dict(patch=q,local=local,global_faces=ids,local_owners=own[local])
   self.rows.append(dict(rank=rank,cells=cells,addr=addr,ninternal=len(nei),bp=bp,maps=maps))
  require(sorted(cells_seen)==list(range(spec.cells)) and np.all((faces_seen>=1)&(faces_seen<=2)) and np.all(faces_seen[spec.internal:]==1),'Complete disjoint cell and signed face cover');require(sorted(physical_seen)==list(range(spec.internal,spec.faces)),'Every physical and empty boundary face exactly once');self.report=dict(rank_cell_counts=[len(r['cells']) for r in self.rows],all_cells_disjoint_complete=True,all_faces_signed_addressed=True,all_boundary_faces_once=True,all_face_vertices_area_vectors_owners_neighbours_proved=True)
 def read(self,label):
  import numpy as np
  s=self.spec;U=np.full((s.cells,3),np.nan);P=np.full(s.cells,np.nan);W=np.full((s.cells,3),np.nan);Phi=np.full(s.faces,np.nan);seen=np.zeros(s.faces,dtype=bool);duplicate=0.;finite=True;boundary={n:dict(phi=np.full(q['count'],np.nan),U=np.full((q['count'],3),np.nan),p=np.full(q['count'],np.nan),shear=np.full((q['count'],3),np.nan)) for n,q in self.patches.items() if n!=s.empty}
  for r in self.rows:
   folder=self.case/f'processor{r["rank"]}'/label;nc=len(r['cells']);raw={n:self.e.text(folder/n) for n in ('U','p','phi','vorticity','wallShearStress')};u=np.asarray(field(raw['U'],nc,3,[0,1,-1,0,0,0,0])).reshape(nc,3);p=np.asarray(field(raw['p'],nc,1,[0,2,-2,0,0,0,0]));w=np.asarray(field(raw['vorticity'],nc,3,[0,0,-1,0,0,0,0])).reshape(nc,3);stress=np.asarray(field(raw['wallShearStress'],nc,3,[0,2,-2,0,0,0,0])).reshape(nc,3);U[r['cells']]=u;P[r['cells']]=p;W[r['cells']]=w;f=np.zeros(len(r['addr']));f[:r['ninternal']]=field(raw['phi'],r['ninternal'],1,[0,3,-1,0,0,0,0]);bb={n:blocks(t) for n,t in raw.items()}
   for name,q in r['bp'].items():
    body=top(bb['phi'][name]);k=kind(body)
    if k=='empty':continue
    # Foundation14 symmetryFvsPatchField::write serializes value; require actual surface data.
    f[q['start']:q['start']+q['count']]=values(body,q['count'],1,'value')
   finite=finite and bool(all(np.isfinite(x).all() for x in (u,p,w,stress,f)));ids=abs(r['addr'])-1;oriented=np.sign(r['addr'])*f;overlap=seen[ids]
   if overlap.any():duplicate=max(duplicate,float(abs(Phi[ids[overlap]]-oriented[overlap]).max()))
   Phi[ids[~overlap]]=oriented[~overlap];seen[ids]=True
   for name,m in r['maps'].items():
    if name==s.empty:require(kind(bb['U'][name])==kind(bb['p'][name])=='empty','Actual empty field types');continue
    q=m['patch'];count=q['count'];slots=m['global_faces']-self.patches[name]['start'];uv=top(bb['U'][name]);pv=top(bb['p'][name]);uk=kind(uv);pk=kind(pv);boundary[name]['phi'][slots]=f[m['local']]
    expected=('fixedValue','zeroGradient') if name=='inlet' else ('zeroGradient','fixedValue') if name=='outlet' else ('rotatingWallVelocity','zeroGradient');require((uk,pk)==expected,'Actual inlet/outlet/body U/p types '+name)
    if name=='cylinder':rotating_wall_configuration(bb['U'][name])
    if re.search(r'\bvalue\s+',uv):ub=np.asarray(values(uv,count,3,'value')).reshape(count,3)
    elif uk=='zeroGradient':ub=u[m['local_owners']]
    else:raise ValueError('Native actual inlet U value required')
    pb=np.asarray(values(pv,count,1,'value')) if re.search(r'\bvalue\s+',pv) else p[m['local_owners']] if pk=='zeroGradient' else None;require(pb is not None,'Native actual outlet p value');boundary[name]['U'][slots]=ub;boundary[name]['p'][slots]=pb
    if pk=='zeroGradient':require(not count or np.max(abs(pb-p[m['local_owners']]))<=1e-10,'Actual zeroGradient p semantics')
    if name=='cylinder':boundary[name]['shear'][slots]=np.asarray(values(top(bb['wallShearStress'][name]),count,3,'value')).reshape(count,3)
  require(seen.all(),'Native all-face coverage');return dict(U=U,p=P,vorticity=W,phi=Phi,boundary=boundary,finite=finite,interface_phi_max_difference=duplicate)

def time_metadata(e,folder,index,dt=.05):
 raw=e.text(folder/'uniform/time');d={n:float(re.search(r'\b'+n+r'\s+('+FLOAT+r')\s*;',raw)[1]) for n in ('value','deltaT','index')};require(abs(d['value']-index*dt)<=1e-10 and abs(d['deltaT']-dt)<=1e-14 and d['index']==index,'Actual native value/dt/index');return d

def residual_history(e,case,first,count,dt=.05,max_outer=30):
 raw=e.text(Path(case)/'solver.log');require(re.search(r'^End\s*$',raw,re.M) and 'Finalising parallel run' in raw and re.search(r'\bnProcs\s*:\s*2\b',raw) and 'FOAM FATAL' not in raw,'Actual completed native MPI2 log');parts=re.split(r'^Time = ('+FLOAT+r')s?\s*$',raw,flags=re.M);require(len(parts)==2*count+1,'Bounded actual time count');rows=[];final=[];counts=[];co=[];names=set();num='(?:'+FLOAT+'|[-+]?(?:nan|inf))'
 for j in range(count):
  i=first+j;t=float(parts[j*2+1]);require(abs(t-i*dt)<=1e-10,'Actual physical time, no pilot replay');body=parts[j*2+2];outer=re.split(r'PIMPLE:\s+Iteration\s+(\d+)\s*',body);seq=[int(outer[k]) for k in range(1,len(outer),2)];require(seq==list(range(1,len(seq)+1)) and len(seq)==max_outer,'Actual fixed outer budget');counts.append(len(seq));last=[]
  for k in range(1,len(outer),2):
   o=int(outer[k]);samples=re.findall(r'Solving for (p|Ux|Uy|Uz), Initial residual = ('+num+r'), Final residual = ('+num+r'), No Iterations (\d+)',outer[k+1],re.I);require(all(any(x[0]==n for x in samples) for n in ('p','Ux','Uy')),'Every actual p/U component correction')
   for q,(n,a,b,it) in enumerate(samples):
    names.add(n);r=dict(step_index=i,time_s=t,outer_iteration=o,solve_index_in_outer=q+1,component=n,initial_residual=float(a),final_residual=float(b),linear_iterations=int(it));rows.append(r)
    if o==seq[-1]:last.append(r)
  final.append(last);c=re.findall(r'Courant Number mean:\s*('+num+r')\s+max:\s*('+num+r')',body,re.I);co.append(float(c[-1][1]) if c else float('nan'))
 return rows,final,counts,co,names

def force_rows(e,path,columns):
 import numpy as np
 raw=e.text(path);rows=[]
 for line in raw.splitlines():
  if not line.strip() or line.lstrip().startswith('#'):continue
  v=list(map(float,line.replace('(',' ').replace(')',' ').split()));require(len(v)==columns,'Native force history shape');rows.append(v)
 require(rows,'Actual native force rows required');a=np.asarray(rows);require(len(np.unique(a[:,0]))==len(a),'Duplicate times within native force file');return a

def body_field_forces(state,g,patches,rho=1.):
 import numpy as np
 out={}
 for name in ('cylinder',):
  q=patches[name];ids=np.arange(q['start'],q['start']+q['count']);sf=np.asarray(g['Sf'])[ids];area=np.linalg.norm(sf,axis=1);b=state['boundary'][name];require(np.all(area>0) and np.linalg.norm(sf.sum(axis=0))<=1e-13,'Closed body actual surface vectors');out[name]=dict(pressure_force_N=np.sum(rho*b['p'][:,None]*sf,axis=0),viscous_force_N=-rho*np.sum(b['shear']*area[:,None],axis=0),finite=bool(np.isfinite(b['p']).all() and np.isfinite(b['shear']).all() and np.isfinite(b['U']).all()),maximum_native_wall_velocity_error=float(abs(b['U']-rotating_wall_velocity(np.asarray(g['Cf'])[ids],sf)).max()),native_wall_speed_min=float(np.linalg.norm(b['U'],axis=1).min()),native_wall_speed_max=float(np.linalg.norm(b['U'],axis=1).max()))
 return out
def rotating_wall_velocity(centres,area_vectors,omega=-2.):
 import numpy as np
 centres=np.asarray(centres);sf=np.asarray(area_vectors);n=sf/np.linalg.norm(sf,axis=1)[:,None];up=np.cross(np.tile([0.,0.,omega],(len(centres),1)),centres);return up-n*np.sum(n*up,axis=1)[:,None]

def rotating_wall_configuration(block):
 body=top(block)
 def vector(name):
  m=re.search(r'\b'+name+r'\s*\(([^()]*)\)\s*;',body);require(m is not None,'Actual rotating-wall '+name);return list(map(float,m[1].split()))
 require(vector('origin')==[0.,0.,0.] and vector('axis')==[0.,0.,1.],'Actual native rotation origin/+z axis')
 m=re.search(r'\bomega\s*\{([^{}]*)\}',block,re.S)
 if m:
  require(re.search(r'\btype\s+constant\s*;',m[1]),'Actual constant omega Function1');v=re.search(r'\bvalue\s+('+FLOAT+r')\s*;',m[1]);require(v is not None,'Actual constant omega value');w=float(v[1])
 else:
  m=re.search(r'\bomega\s+('+FLOAT+r')\s*\[rad/s\]\s*;',body);require(m is not None,'Explicit angular-velocity SI unit');w=float(m[1])
 require(w==-2.,'Actual clockwise omega-2rad/s');return w
