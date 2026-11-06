"""Prepare/check annular Couette pilot files; never launches CFD."""
from pathlib import Path
import argparse,hashlib,json,math,re,shutil
RUN='2026-11-02-annular-couette-v001';PATCHES=('inner','outer');EMPTY='frontAndBackPlanes'
CELLS=1536;POINTS=3200;FACES=6208;INTERNAL=3008;DEPTH=.01;DT=.05
TARGET_SHA='4c3f6bdb49ffb208230a7c89f574d9c030fd0d97de800e256f791ecf64abfff9'
MESH_SHA='13d68cc96c8ea8f55f82f212315245c6e7a377223f3fcc87813328912d13f752';NPZ_SHA='d9fed1a699d1558e32de2fce7cb81776bf7294da8e4f405274d930c2625b1a97'
SETUP_SHA='bcb50febd503aaaf9f951fce98405a3cfa97ea47c2454d76e05fdd503a4cbfee'
RESEARCH_SHA='3197de4c5b8c20e83c45946dacdaf929f22a38fa9c5b9e26a2e9b995301d2168'
FLOAT=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def require(v,m):
 if not v:raise ValueError(m)
def root():return Path(__file__).resolve().parents[3]
def run_root():return root()/'native/runs'/RUN
def private():return root()/'local/annular-couette-v001'
def checked_case(p):
 p=Path(p).resolve();require(p==(run_root()/'openfoam-pilot-v001').resolve(),'Dedicated annular pilot only');return p
INSTALLED_SOURCES=tuple(r['path'] for r in json.loads((private()/'official-source-evidence-v001/source-manifest.json').read_text()))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write_new(p,s):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x',encoding='ascii',newline='\n') as f:f.write(s.rstrip()+'\n')
def json_new(p,v):write_new(p,json.dumps(v,indent=2,allow_nan=False))
def copied(a,b):
 require(not b.exists(),'Preserve existing '+str(b));b.parent.mkdir(parents=True,exist_ok=True);h=sha(a);shutil.copy2(a,b);require(sha(a)==sha(b)==h,'Copy hash mismatch')
def clean(s):return re.sub(r'/\*.*?\*/|//[^\n]*','',s,flags=re.S)
def header(n,c='dictionary'):return 'FoamFile\n{\n version 2.0;\n format ascii;\n class '+c+';\n object '+n+';\n}\n\n'
def dictionaries():
 d={'system/controlDict':'''solver incompressibleFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 0.25;
deltaT 0.05;
adjustTimeStep false;
writeControl timeStep;
writeInterval 1;
purgeWrite 0;
writeFormat ascii;
writePrecision 14;
writeCompression off;
timeFormat general;
timePrecision 10;
runTimeModifiable false;
''','system/fvSchemes':'''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss upwind; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''','system/fvSolution':'''solvers
{
 p { solver GAMG; tolerance 1e-10; relTol 0; smoother GaussSeidel; }
 pFinal { $p; relTol 0; }
 U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0; }
 UFinal { $U; relTol 0; }
}
PIMPLE
{
 nOuterCorrectors 10;
 nCorrectors 2;
 nNonOrthogonalCorrectors 0;
 momentumPredictor yes;
 consistent yes;
 pRefPoint (0.05092584598549641 0.00250182638511251 0);
 pRefValue 0;
}
''','system/decomposeParDict':'numberOfSubdomains 2;\nmethod simple;\nsimpleCoeffs { n (2 1 1); delta 0.001; }\n','constant/physicalProperties':'viscosityModel constant;\nnu 0.000125;\n','constant/momentumTransport':'simulationType laminar;\n','0/U':'''dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
 inner { type rotatingWallVelocity; origin (0 0 0); axis (0 0 1); omega 1 [rad/s]; }
 outer { type noSlip; }
 frontAndBackPlanes { type empty; }
}
''','0/p':'''dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
 inner { type fixedFluxPressure; }
 outer { type fixedFluxPressure; }
 frontAndBackPlanes { type empty; }
}
'''}
 vortex='''nativeVorticity
{
 type vorticity;
 libs ("libfieldFunctionObjects.so");
 field U;
 result vorticity;
 executeControl timeStep;
 executeInterval 1;
 writeControl timeStep;
 writeInterval 1;
}
'''
 d['system/functions']=vortex
 d['system/preflightDiagnostics']='''nativeCentres { type writeCellCentres; libs ("libfieldFunctionObjects.so"); }
nativeVolumes { type writeCellVolumes; libs ("libfieldFunctionObjects.so"); }
nativeInitialFields { type writeObjects; libs ("libutilityFunctionObjects.so"); objects (U p phi); writeOption anyWrite; }
'''+vortex
 return {n:header(Path(n).name,'volVectorField' if n=='0/U' else 'volScalarField' if n=='0/p' else 'dictionary')+v for n,v in d.items()}
def counted(p):
 t=clean(Path(p).read_text());h=re.search(r'FoamFile\s*\{([^{}]*)\}',t);require(h is not None and re.search(r'\bformat\s+ascii\s*;',h[1]),'ASCII native list required');m=re.fullmatch(r'\s*(\d+)\s*\((.*)\)\s*',t[h.end():],re.S);require(m is not None,'Invalid counted native list '+str(p));return int(m[1]),m[2]
def values(t,count,components,label='internalField'):
 m=re.search(r'\b'+label+r'\s+uniform\s+([^;]+);',t)
 if m:
  row=list(map(float,m[1].replace('(',' ').replace(')',' ').split()));require(len(row)==components,'Uniform component count');return [row[:] for _ in range(count)] if components>1 else [row[0]]*count
 kind='vector' if components==3 else 'scalar';m=re.search(r'\b'+label+r'\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',t,re.S);require(m is not None and int(m[1])==count,'Native list count');a=list(map(float,m[2].replace('(',' ').replace(')',' ').split()));require(len(a)==count*components,'Native value length');return [a[i:i+components] for i in range(0,len(a),components)] if components>1 else a
def field(p,count,components,dims):
 t=clean(Path(p).read_text());m=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',t);require(m is not None and list(map(int,m[1].split()))==dims,'Field dimensions');return values(t,count,components)
def blocks(p):
 t=clean(Path(p).read_text());require('boundaryField' in t,'Missing boundaryField');return dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',t.split('boundaryField',1)[1]))
def native_poly(c):
 a=[c/n for n in ('0/polyMesh','constant/polyMesh') if (c/n/'points').is_file()];require(len(a)==1,'Exactly one initial mesh');return a[0]
def patches(poly):
 nb,b=counted(poly/'boundary');out={n:dict(kind=re.search(r'\btype\s+(\w+)\s*;',v)[1],count=int(re.search(r'\bnFaces\s+(\d+)\s*;',v)[1]),start=int(re.search(r'\bstartFace\s+(\d+)\s*;',v)[1])) for n,v in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',b)}
 require(nb==3 and set(out)==set(PATCHES)|{EMPTY},'Annular physical patch set');require(all(out[n]['kind']=='wall' and out[n]['count']==64 for n in PATCHES),'64native faces perwall');require(out[EMPTY]['kind']=='empty' and out[EMPTY]['count']==3072,'Empty face count');return out
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def native_geometry(poly):
 n,b=counted(poly/'points');pts=[list(map(float,x.split())) for x in re.findall(r'\(([^()]*)\)',b)];require(n==POINTS==len(pts),'Native point count');n,b=counted(poly/'faces');faces=[list(map(int,x.split())) for k,x in re.findall(r'(\d+)\s*\(([^()]*)\)',b)];require(n==FACES==len(faces) and all(len(x)==4 and len(set(x))==4 for x in faces),'Native quadrilateral faces')
 def labels(name,expected):
  n,b=counted(poly/name);a=list(map(int,b.split()));require(n==expected==len(a) and all(0<=v<CELLS for v in a),'Native '+name+' labels');return a
 own=labels('owner',FACES);nei=labels('neighbour',INTERNAL);V=[0.]*CELLS;Sf=[];Cf=[];vertices=[set() for _ in range(CELLS)]
 for i,ids in enumerate(faces):
  v=[pts[j] for j in ids];sf=[0.,0.,0.];volume=0.
  for a,b in zip(v,v[1:]+v[:1]):sf=[x+.5*y for x,y in zip(sf,cross(a,b))]
  for j in range(1,3):volume+=dot(v[0],cross(v[j],v[j+1]))/6
  V[own[i]]+=volume;vertices[own[i]].update(ids)
  if i<len(nei):V[nei[i]]-=volume;vertices[nei[i]].update(ids)
  Sf.append(sf);Cf.append([sum(q[j] for q in v)/4 for j in range(3)])
 require(all(math.isfinite(v) and v>0 for v in V) and all(len(v)==8 for v in vertices),'Positive native eight-vertex cells');return dict(points=pts,faces=faces,owner=own,neighbour=nei,volumes=V,Sf=Sf,Cf=Cf,cell_vertices=[sorted(v) for v in vertices])
def verify_inputs(c):
 require(sha(private()/'official-source-evidence-v001/source-manifest.json')==RESEARCH_SHA,'Frozen primary research manifest changed')
 m=json.loads((c/'provenance/prepared-input-sha256.json').read_text());require(all(sha(c/n)==d for n,d in m.items()),'Prepared immutable input changed');return m
def prepare(mesh,case):
 import numpy as np
 case=checked_case(case);mesh=Path(mesh).resolve();require(not case.exists(),'Fresh case required');require(mesh==run_root()/'mesh/annular-couette.msh' and sha(mesh)==MESH_SHA and sha(mesh.parent/'geometry.npz')==NPZ_SHA,'Selected mesh/NPZ hash');require(sha(mesh.parent/'setup.json')==SETUP_SHA,'Exact shared physical setup hash');setup=json.loads((mesh.parent/'setup.json').read_text());require(setup['cells']==CELLS and setup['kinematic_viscosity_m2_s']==.000125 and setup['density_kg_m3']==1 and setup['time_step_s']==DT and setup['pilot_time_steps']==5 and setup['time_steps']==400 and setup['rotation_axis_direction']==[0,0,1],'Physics/time setup');require(sha(private()/'TARGETS.md')==TARGET_SHA,'Frozen target hash');require(sha(private()/'official-source-evidence-v001/source-manifest.json')==RESEARCH_SHA,'Official research manifest hash');case.mkdir(parents=True)
 for n,t in dictionaries().items():write_new(case/n,t)
 for n in ('U','p'):write_new(case/'provenance/initial-fields'/n,dictionaries()['0/'+n])
 copied(mesh,case/'input-annular-couette.msh')
 for p in (mesh.parent/'setup.json',mesh.parent/'geometry.npz',private()/'TARGETS.md',Path(__file__),private()/'run-openfoam-v001.sh',private()/'Launch-OpenFOAM-v001.ps1'):copied(p,case/'provenance'/p.name)
 copied(private()/'official-source-evidence-v001/source-manifest.json',case/'provenance/official-research-source-manifest.json')
 with np.load(mesh.parent/'geometry.npz',allow_pickle=False) as f:g={n:np.asarray(f[n]).tolist() for n in ('nodes','quads','centres','areas','grid_shape','ix','iy','face_nodes','face_owner','face_neighbour','face_zone','face_centres','face_area_vectors')}
 require(len(g['nodes'])==1600 and len(g['quads'])==CELLS and g['grid_shape']==[24,64],'Actual source array shape');json_new(case/'provenance/source-geometry.json',g)
 json_new(case/'declared-comparison.json',dict(phase='five_step_transient_pilot',physical_setup=setup,initial_cell_U=[0,0,0],initial_p=0,mesh_fixed=True,frame='inertial',solver='Foundation14 incompressibleFluid',ddt='Euler',momentum_convection='Gauss upwind',dt_s=DT,pilot_steps=5,pilot_end_s=.25,prospective_total_steps=400,production_admitted=False,outer_correctors=10,pressure_correctors=2,last_outer_initial_residual_target=1e-6,last_outer_final_linear_residual_target=1e-8,local_divergence_scale_s=1,local_divergence_max_target=1e-5,local_divergence_RMS_target=1e-6,absolute_wall_phi_target_m3_s=1e-12,inner_wall_U_target_m_s=1e-10,outer_wall_U_target_m_s=1e-12,pressure_reference_point_m=setup['pressure_reference_xyz_m'],pressure_reference_value=0,pressure_gauge_clarification='pRefPoint identifies sourcecell0 acrossMPI; pRefCell denotes master-local only',target_document_sha256=TARGET_SHA,physical_validation=False))
 json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file() and p.relative_to(case).as_posix() not in ('0/U','0/p')});json_new(case/'run-status.json',dict(stage='prepared_only',solver_started=False));verify_inputs(case);return dict(case=str(case),prepared=True,solver_started=False)
def validate_native(case):
 c=checked_case(case);verify_inputs(c);require(not (c/'solver.log').exists(),'Preflight before CFD');log=(c/'mesh-check.log').read_text();require('Mesh OK.' in log and 'checkMesh -allGeometry -allTopology' in log and not re.search(r'FOAM FATAL|Failed\s+\d+\s+mesh checks',log),'Full checkMesh');nonorth=float(re.search(r'Mesh non-orthogonality Max:\s*('+FLOAT+')',log)[1]);skew=float(re.search(r'Max skewness\s*=\s*('+FLOAT+')',log)[1]);require(nonorth<=75 and skew<=4,'Native quality')
 poly=native_poly(c);pa=patches(poly);g=native_geometry(poly);src=json.loads((c/'provenance/source-geometry.json').read_text());lookup={tuple(round(v,10) for v in q):i for i,q in enumerate(src['nodes'])};pointmap=[];planes=[]
 for x,y,z in g['points']:
  key=(round(x,10),round(y,10));require(key in lookup,'Nativepoint absent source');i=lookup[key];require(math.hypot(x-src['nodes'][i][0],y-src['nodes'][i][1])<=1e-10 and min(abs(z-.005),abs(z+.005))<=1e-12,'Sourcepoint/extrusion');pointmap.append(i);planes.append(1 if z>0 else -1)
 require(len(set(zip(pointmap,planes)))==3200,'Eachsourcevertexbothplanes');celllookup={frozenset(q):i for i,q in enumerate(src['quads'])};n2s=[]
 for ids in g['cell_vertices']:
  key=frozenset(pointmap[j] for j in ids);require(key in celllookup,'Nativecell sourcevertexset');n2s.append(celllookup[key]);require(len({pointmap[j] for j in ids if planes[j]>0})==4 and len({pointmap[j] for j in ids if planes[j]<0})==4,'Nativecell bothplanes')
 require(sorted(n2s)==list(range(CELLS)),'Native/sourcecell bijection');inverse=[0]*CELLS
 for n,s in enumerate(n2s):inverse[s]=n
 C=field(c/'0/C',CELLS,3,[0,1,0,0,0,0,0]);Vc=field(c/'0/Vc',CELLS,1,[0,3,0,0,0,0,0]);require(all(math.dist(v[:2],src['centres'][n2s[i]])<=1e-10 and abs(v[2])<=1e-12 for i,v in enumerate(C)),'NativeC/sourcecentres');require(all(v>0 and abs(v-g['volumes'][i])<=1e-14 and abs(v/DEPTH-src['areas'][n2s[i]])<=1e-12 for i,v in enumerate(Vc)),'NativeVc/oriented/sourceareas')
 edge_lookup={frozenset(e):i for i,e in enumerate(src['face_nodes'])};fmap=[-1]*FACES;fsign=[0]*FACES;seen=[]
 for i,ids in enumerate(g['faces']):
  key=frozenset(pointmap[j] for j in ids)
  if len(key)==2:
   require(key in edge_lookup,'Nativeface missing sourceedge');s=edge_lookup[key];fmap[i]=s;expected=[v*DEPTH for v in src['face_area_vectors'][s]]+[0.];sign=1 if dot(g['Sf'][i],expected)>0 else -1;fsign[i]=sign;require(math.dist(g['Sf'][i],[sign*v for v in expected])<=1e-12,'Source/nativeface Sf');so=src['face_owner'][s];sn=src['face_neighbour'][s];require(n2s[g['owner'][i]]==(so if sign>0 else sn),'Sourceface owner')
   if i<INTERNAL:require(n2s[g['neighbour'][i]]==(sn if sign>0 else so),'Sourceface neighbour')
   else:require(sn==-1 and sign==1,'Nativewall orientation')
   seen.append(s)
 require(sorted(seen)==list(range(3136)),'Eachsourceedge one nativeverticalface');U=field(c/'0/U',CELLS,3,[0,1,-1,0,0,0,0]);P=field(c/'0/p',CELLS,1,[0,2,-2,0,0,0,0]);require(max(abs(v) for r in U for v in r)<=1e-12 and max(abs(v) for v in P)<=1e-12,'Actualinitial cellrest/p0');ub=blocks(c/'0/U');pb=blocks(c/'0/p');bc={}
 for name in PATCHES:
  q=pa[name];ids=range(q['start'],q['start']+q['count']);kind='rotatingWallVelocity' if name=='inner' else 'noSlip';require(re.search(r'\btype\s+'+kind+r'\s*;',ub[name]) and re.search(r'\btype\s+fixedFluxPressure\s*;',pb[name]),'Native wall/pressure types');actual=values(ub[name],64,3,'value') if re.search(r'\bvalue\s+',ub[name]) else [[0.,0.,0.]]*64;error=0.
  for j,fi in enumerate(ids):
   cf=g['Cf'][fi];sf=g['Sf'][fi];a=math.sqrt(dot(sf,sf));n=[v/a for v in sf];u=[-cf[1],cf[0],0.] if name=='inner' else [0.,0.,0.];target=[x-y*dot(n,u) for x,y in zip(u,n)];error=max(error,max(abs(x-y) for x,y in zip(actual[j],target)))
  require(error<=(1e-10 if name=='inner' else 1e-12),'Actual wallvelocity projection');bc[name]=dict(max_abs_velocity_error_m_s=error,actual_faces=64)
 owning=[]
 for rank in (0,1):
  n,b=counted(c/f'processor{rank}/0/polyMesh/cellProcAddressing');a=list(map(int,b.split()));require(n==len(a),'Rankcelladdress count')
  if inverse[0] in a:owning.append(dict(rank=rank,local_cell=a.index(inverse[0]),global_native_cell=inverse[0],source_cell=0))
 require(len(owning)==1,'pRefPoint uniqueMPIowner');mapping=dict(native_cell_to_source_cell=n2s,source_native_cell_index=inverse,native_point_to_source_node=pointmap,native_point_z_sign=planes,native_face_to_source_edge=fmap,native_face_source_orientation=fsign);json_new(c/'provenance/native-mapping.json',mapping)
 receipt=dict(passed=True,solver_started=False,cells=CELLS,points=POINTS,faces=FACES,internal_faces=INTERNAL,physical_faces=128,empty_faces=3072,maximum_nonorthogonality_degrees=nonorth,maximum_skewness=skew,min_positive_volume_m3=min(Vc),total_volume_m3=sum(Vc),native_mesh_sha256={n:sha(poly/n) for n in ('points','faces','owner','neighbour','boundary')},actual_native_C_Vc_sha256={n:sha(c/'0'/n) for n in ('C','Vc')},actual_initial_cell_rest=True,actual_boundary_velocity_proof=bc,source_cell0_pressure_reference_unique_rank=owning[0],runtime_parallel_pressureReference_constructed_successfully=True,full_source_cell_point_face_vertex_mapping=True);json_new(c/'provenance/preflight-validation.json',receipt);return receipt

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,required=True);g=ap.add_mutually_exclusive_group(required=True);g.add_argument('--mesh',type=Path);g.add_argument('--validate-native',action='store_true');g.add_argument('--verify-inputs',action='store_true');a=ap.parse_args();print(json.dumps(prepare(a.mesh,a.case) if a.mesh else validate_native(a.case) if a.validate_native else dict(verified_files=len(verify_inputs(checked_case(a.case)))),indent=2))
if __name__=='__main__':main()
