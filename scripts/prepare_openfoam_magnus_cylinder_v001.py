"""Prepare/check the frozen rotating-cylinder transient pilot; never invokes native programs or CFD."""
from pathlib import Path
import argparse,json,re,math
import magnus_openfoam_common_v001 as common
from magnus_openfoam_common_v001 import MeshSpec,require,sha,write_new,json_new,copied,header,field,blocks,top,kind,values,counted,native_poly,patch_table,native_geometry,source_mapping
COMMON_SHA='587991b8159edce429e2133170938bbc7c40ece58a0197f1ef63ba0147c3490f'
TARGET_SHA='87dfdeb12300ba0d6dfaefdedf58d0c08ef984673522f4432b1fe90807f0dbbe'
PLAN_SHA='0e811eb689652575cdb036560b65a2a3f3c0cc2abd604b0c6073bac8a390d49d'
MESH_SHA='b6db8e780813a5e56f73896a4b0a462149ad4e351144946c74d5cc330a640501';NPZ_SHA='7490cd34096043a789a8a7abdc01ea806f4f533a9b9d437ac63a8158c190092a';SETUP_SHA='936c01fe314266f1acd9c3c26640801b7c5f784ed2fb6337fe990ef33f1711f8'
RUN='2026-11-02-magnus-cylinder-v001';PATCH_COUNTS=dict(inlet=64,outlet=64,cylinder=128);ZONES=dict(inlet=5,outlet=6,cylinder=4);SPEC=MeshSpec(8192,16640,32896,16256,PATCH_COUNTS)
INSTALLED_SOURCES=(
'applications/modules/incompressibleFluid/incompressibleFluid.C','applications/modules/incompressibleFluid/correctPressure.C',
'src/finiteVolume/fields/fvPatchFields/derived/rotatingWallVelocity/rotatingWallVelocityFvPatchVectorField.C',
'src/finiteVolume/fields/fvsPatchFields/constraint/symmetry/symmetryFvsPatchField.C',
'src/finiteVolume/cfdTools/general/solutionControl/pimpleControl/pimpleLoop/pimpleLoop.C','src/finiteVolume/cfdTools/general/solutionControl/pimpleControl/pimpleNoLoopControl/pimpleNoLoopControl.C','src/finiteVolume/cfdTools/general/solutionControl/pisoControl/pisoControl.C','src/finiteVolume/fvMatrices/fvMatrix/fvMatrixSolve.C',
'src/functionObjects/forces/forces/forces.C','src/functionObjects/forces/forcesBase/forcesBase.C','src/functionObjects/forces/forceCoeffs/forceCoeffs.C','src/functionObjects/field/wallShearStress/wallShearStress.C',
'src/OpenFOAM/db/functionObjects/functionObject/functionObject.C','src/OpenFOAM/db/functionObjects/functionObjectList/functionObjectList.C','src/functionObjects/field/fieldExpression/fieldExpression.C',
'src/functionObjects/utilities/writeObjects/writeObjects.C','src/functionObjects/field/writeCellCentres/writeCellCentres.C','src/functionObjects/field/writeCellVolumes/writeCellVolumes.C','src/functionObjects/field/vorticity/vorticity.C','src/finiteVolume/finiteVolume/fvc/fvcCurl.C','applications/utilities/postProcessing/foamPostProcess/foamPostProcess.C')
def root():return Path(__file__).resolve().parents[3]
def run_root():return root()/'native/runs'/RUN
def private():return root()/'local/magnus-cylinder-v001'
def checked_case(p):
 p=Path(p).resolve();require(p==(run_root()/'openfoam-pilot-v001').resolve(),'Dedicated frozen Magnus pilot only');return p

def force_function(name,typ,body,centre):
 s=f'''{name}
{{
 type {typ};
 libs ("libforces.so");
 patches ({body});
 rho rhoInf; rhoInf 1; p p; U U; pRef 0; porosity false;
 CofR ({centre[0]} {centre[1]} 0);
 log false;
 executeControl timeStep; executeInterval 1;
 writeControl timeStep; writeInterval 1;
'''
 if typ=='forceCoeffs':s+=' liftDir (0 1 0); dragDir (1 0 0); pitchAxis (0 0 1); magUInf 0.1; lRef 0.1; Aref 0.001;\n'
 return s+'}\n'
def field_functions():return '''nativeVorticity { type vorticity; libs ("libfieldFunctionObjects.so"); field U; result vorticity; executeControl timeStep; executeInterval 1; writeControl timeStep; writeInterval 1; }
nativeBodyShear { type wallShearStress; libs ("libfieldFunctionObjects.so"); patches (cylinder); log false; executeControl timeStep; executeInterval 1; writeControl timeStep; writeInterval 1; }
'''
def functions(probe=False):
 s=''
 for body,centre in [('cylinder',(0,0))]:
  s+=force_function(('probe' if probe else '')+body[0].upper()+body[1:]+'Forces' if probe else body+'Forces','forces',body,centre)
  s+=force_function(('probe' if probe else '')+body[0].upper()+body[1:]+'Coefficients' if probe else body+'Coefficients','forceCoeffs',body,centre)
 return s+field_functions()
def dictionaries():
 d={'system/controlDict':'''solver incompressibleFluid;
startFrom startTime; startTime 0; stopAt endTime; endTime 0.1;
deltaT 0.02; adjustTimeStep false;
writeControl timeStep; writeInterval 1; purgeWrite 0;
writeFormat ascii; writePrecision 14; writeCompression off;
timeFormat general; timePrecision 10; runTimeModifiable false;
'''+ 'functions\n{\n'+functions()+'}\n','system/fvSchemes':'''ddtSchemes { default Euler; }
gradSchemes { default leastSquares; }
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
PIMPLE { nOuterCorrectors 30; nCorrectors 2; nNonOrthogonalCorrectors 0; momentumPredictor yes; consistent yes; }
''','system/decomposeParDict':'numberOfSubdomains 2;\nmethod simple;\nsimpleCoeffs { n (2 1 1); delta 0.001; }\n','constant/physicalProperties':'viscosityModel constant;\nnu 0.0005;\n','constant/momentumTransport':'simulationType laminar;\n','0/U':'''dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
 inlet { type fixedValue; value uniform (0.1 0 0); }
 outlet { type zeroGradient; }
 cylinder { type rotatingWallVelocity; origin (0 0 0); axis (0 0 1); omega -2 [rad/s]; }
 frontAndBackPlanes { type empty; }
}
''','0/p':'''dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
 inlet { type zeroGradient; }
 outlet { type fixedValue; value uniform 0; }
 cylinder { type zeroGradient; }
 frontAndBackPlanes { type empty; }
}
''','system/initialDiagnostics':'''nativeCentres { type writeCellCentres; libs ("libfieldFunctionObjects.so"); }
nativeVolumes { type writeCellVolumes; libs ("libfieldFunctionObjects.so"); }
nativeInitialFields { type writeObjects; libs ("libutilityFunctionObjects.so"); objects (U p phi); writeOption anyWrite; }
'''+field_functions(),'system/forcePreflight':'''nativeInitialFields { type writeObjects; libs ("libutilityFunctionObjects.so"); objects (U p phi); writeOption anyWrite; }
'''+functions(True)}
 return {n:header(Path(n).name,'volVectorField' if n=='0/U' else 'volScalarField' if n=='0/p' else 'dictionary')+v for n,v in d.items()}

def verify_inputs(case):
 case=checked_case(case);require(sha(Path(common.__file__))==COMMON_SHA,'Frozen common reader changed');m=json.loads((case/'provenance/prepared-input-sha256.json').read_text());require(all(sha(case/n)==h for n,h in m.items()),'Prepared source/control/physics input changed');return m

def prepare(mesh,case):
 import numpy as np
 case=checked_case(case);mesh=Path(mesh).resolve();require(not case.exists() and mesh==(run_root()/'mesh/magnus-cylinder.msh').resolve(),'Fresh dedicated case and selected shared mesh');require(sha(mesh)==MESH_SHA and sha(mesh.parent/'geometry.npz')==NPZ_SHA and sha(mesh.parent/'setup.json')==SETUP_SHA,'Actual accepted shared input hashes');require(sha(private()/'TARGETS.md')==TARGET_SHA and sha(private()/'EXECUTION-PLAN-v001.md')==PLAN_SHA and sha(Path(common.__file__))==COMMON_SHA,'Frozen target/plan/common source');setup=json.loads((mesh.parent/'setup.json').read_text());require(setup['cells']==8192 and setup['density_kg_m3']==1 and setup['kinematic_viscosity_m2_s']==.0005 and setup['inlet_velocity_m_s']==[.1,0] and setup['time_step_s']==.02 and setup['time_steps']==500 and setup['pilot_time_steps']==5,'Actual physical setup')
 case.mkdir(parents=True)
 for n,t in dictionaries().items():write_new(case/n,t)
 for n in ('U','p'):write_new(case/'provenance/prescribed-initial-fields'/n,dictionaries()['0/'+n])
 copied(mesh,case/'input-magnus-cylinder.msh')
 for p in (mesh.parent/'geometry.npz',mesh.parent/'setup.json',private()/'TARGETS.md',private()/'EXECUTION-PLAN-v001.md',Path(__file__),Path(common.__file__),private()/'run-openfoam-v001.sh',private()/'Launch-OpenFOAM-v001.ps1',private()/'OPENFOAM-PREREQUISITES-v001.md'):copied(p,case/'provenance'/p.name)
 with np.load(mesh.parent/'geometry.npz',allow_pickle=False) as f:src={n:np.asarray(f[n]).tolist() for n in ('nodes','quads','centres','areas','grid_shape','ix','iy','face_nodes','face_owner','face_neighbour','face_zone','face_area_vectors')}
 require(src['grid_shape']==[64,128] and len(src['quads'])==8192 and len(src['nodes'])==8320 and len(src['face_nodes'])==16512,'Accepted wrapped O-grid source shape');require(set(zip(src['ix'],src['iy']))=={(i,j) for j in range(64) for i in range(128)},'Complete wrapped angular/radial source slots');json_new(case/'provenance/source-geometry.json',src);json_new(case/'provenance/installed-source-list.json',list(INSTALLED_SOURCES));json_new(case/'provenance/official-primary-source-urls.json',{n:'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/'+n for n in INSTALLED_SOURCES})
 json_new(case/'declared-comparison.json',dict(phase='five_step_transient_pilot',initial_cell_U=[0,0,0],initial_cell_p=0,initial_state_kind='prescribed rest with active inlet and rotating wall; native initial phi retained, not solved-conservative',solved_gate_first_index=1,pilot_steps=5,dt_s=.02,pilot_final_s=.1,total_steps_prospective=500,additional_steps_prospective=495,final_s_prospective=10,outer_correctors=30,pressure_correctors=2,force_body_names=['cylinder'],cylinder_omega_rad_s=-2,force_reference_rho=1,force_reference_U=.1,force_reference_D=.1,force_reference_A_m2=.001,force_denominator_N=5e-6,force_coefficient_relative_absolute_target=1e-6,force_norm='abs(delta coefficient)<=1e-6*max(1,abs(native coefficient))',primary_restart_inputs=['U','p','phi','time','mesh','physics','sources'],derived_restart_outputs=['vorticity','wallShearStress','force function-object startup outputs'],target_sha256=TARGET_SHA,plan_sha256=PLAN_SHA,mesh_sha256=MESH_SHA,geometry_sha256=NPZ_SHA,setup_sha256=SETUP_SHA,common_source_sha256=COMMON_SHA,physical_validation=False,production_admitted=False))
 json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file() and p.relative_to(case).as_posix() not in ('0/U','0/p')});json_new(case/'run-status.json',dict(stage='prepared_only',solver_started=False,production_admitted=False));verify_inputs(case);return dict(case=str(case),prepared=True,solver_started=False,prepared_input_manifest_sha256=sha(case/'provenance/prepared-input-sha256.json'))

def validate_native(case):
 case=checked_case(case);verify_inputs(case);require(not (case/'solver.log').exists(),'Preflight before CFD');reg=json.loads((case/'provenance/native-registration-process-status.json').read_text())
 for label,ranks in [('serial',1),('parallel',2)]:
  r=reg[label];text=(case/r['log']).read_text();require(r['return_code']==0 and r['mpi_ranks']==ranks and r['log_sha256']==sha(case/r['log']) and re.search(r'^End\s*$',text,re.M) and 'FOAM FATAL' not in text,'Actual completed registration process/log');times=re.findall(r'^Time = ('+common.FLOAT+r')s?\s*$',text,re.M);require(times and all(float(t)==0 for t in times),'Actual native time0 only in preflight')
  if ranks==2:require(re.search(r'\bnProcs\s*:\s*2\b',text),'Actual parallel registration')
 text=(case/'mesh-check.log').read_text();require('Mesh OK.' in text and 'checkMesh -allGeometry -allTopology' in text,'Actual full checkMesh');nonorth=float(re.search(r'Mesh non-orthogonality Max:\s*('+common.FLOAT+')',text)[1]);skew=float(re.search(r'Max skewness\s*=\s*('+common.FLOAT+')',text)[1]);require(nonorth<=75 and skew<=4,'Native quality targets');poly=native_poly(case);pa=patch_table(poly);require(set(pa)==set(PATCH_COUNTS)|{SPEC.empty},'Exact native physical patch set')
 for n,q in pa.items():require(q['count']==(16384 if n==SPEC.empty else PATCH_COUNTS[n]) and q['kind']==('empty' if n==SPEC.empty else 'wall' if n=='cylinder' else 'patch'),'Actual native patch type/count '+n)
 g=native_geometry(poly,SPEC);C=field((case/'0/C').read_text(),8192,3,[0,1,0,0,0,0,0]);V=field((case/'0/Vc').read_text(),8192,1,[0,3,0,0,0,0,0]);src=json.loads((case/'provenance/source-geometry.json').read_text());mapping=source_mapping(g,src,C,V,SPEC,pa,ZONES);U=field((case/'0/U').read_text(),8192,3,[0,1,-1,0,0,0,0]);P=field((case/'0/p').read_text(),8192,1,[0,2,-2,0,0,0,0]);require(all(math.isfinite(v) and abs(v)<=1e-12 for row in U for v in row) and all(math.isfinite(v) and abs(v)<=1e-12 for v in P),'Actual prescribed initial cell rest/p0');ub=blocks((case/'0/U').read_text());pb=blocks((case/'0/p').read_text());require(kind(ub['inlet'])=='fixedValue' and all(math.dist(v,[.1,0,0])<=1e-12 for v in values(top(ub['inlet']),64,3,'value')) and kind(pb['inlet'])=='zeroGradient','Actual initial inlet types/value');require(kind(ub['outlet'])=='zeroGradient' and kind(pb['outlet'])=='fixedValue' and max(map(abs,values(top(pb['outlet']),64,1,'value')))<=1e-12,'Actual initial outlet p0')
 q=pa['cylinder'];ids=list(range(q['start'],q['start']+q['count']));require(kind(ub['cylinder'])=='rotatingWallVelocity' and kind(pb['cylinder'])=='zeroGradient','Actual rotating wall/open-flow pressure types');common.rotating_wall_configuration(ub['cylinder']);actual=values(top(ub['cylinder']),128,3,'value');wall_error=0.;wall_normal=0.
 for v,f in zip(actual,ids):
  sf=g['Sf'][f];mag=math.sqrt(sum(z*z for z in sf));normal=[z/mag for z in sf];x,y,z=g['Cf'][f];up=[2*y,-2*x,0.];expected=[a-b*sum(c*d for c,d in zip(normal,up)) for a,b in zip(up,normal)];wall_error=max(wall_error,max(abs(a-b) for a,b in zip(v,expected)));wall_normal=max(wall_normal,abs(sum(a*b for a,b in zip(v,normal))))
 require(wall_error<=1e-10 and wall_normal<=1e-12,'Actual constructor-evaluated projected Omega cross facet velocity')
 cells=[]
 for r in (0,1):
  n,b=counted(case/f'processor{r}/0/polyMesh/cellProcAddressing');a=list(map(int,b.split()));require(n==len(a),'Rank address count');cells.extend(a)
 require(sorted(cells)==list(range(8192)),'Actual MPI cell cover disjoint/complete')
 probe={}
 for name in ('Cylinder',):
  arrays=[]
  for typ,columns,filename in [('Forces',13,'forces.dat'),('Coefficients',6,'forceCoeffs.dat')]:
   p=case/f'postProcessing/probe{name}{typ}/0/{filename}';rows=[list(map(float,line.replace('(',' ').replace(')',' ').split())) for line in p.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')];require(len(rows)==1 and len(rows[0])==columns and rows[0][0]==0 and all(map(math.isfinite,rows[0])),'Actual time0 native force registration row');arrays.append(rows[0])
  f,k=arrays;cd=(f[1]+f[4])/5e-6;cl=(f[2]+f[5])/5e-6;require(abs(cd-k[2])<=1e-6*max(1,abs(k[2])) and abs(cl-k[3])<=1e-6*max(1,abs(k[3])),'Native registered raw force/coefficient references');probe[name]=dict(actual_native_time0_force_row=f,actual_native_time0_coefficient_row=k,registration_only_not_nonzero_validation=True)
 primary={};derived={}
 for folder in (case/'0',case/'processor0/0',case/'processor1/0'):
  for p in sorted(folder.rglob('*')):
   if p.is_file():
    rel=p.relative_to(case).as_posix()
    if p.name in ('vorticity','wallShearStress'):
     baseline='provenance/native-initial-derived-baseline/'+rel;copied(p,case/baseline);nc=8192 if folder==case/'0' else counted(folder/'polyMesh/cellProcAddressing')[0];derived[rel]=dict(sha256=sha(p),baseline_path=baseline,cells=nc,dimensions=[0,0,-1,0,0,0,0] if p.name=='vorticity' else [0,2,-2,0,0,0,0])
    else:primary[rel]=sha(p)
 require(all((case/f'processor{r}/0'/n).is_file() for r in (0,1) for n in ('U','p','phi','vorticity','wallShearStress')),'Actual initial rank fields');json_new(case/'provenance/native-mapping.json',mapping);json_new(case/'provenance/native-initial-primary-sha256.json',primary);json_new(case/'provenance/native-initial-derived-sha256.json',derived)
 receipt=dict(passed=True,solver_started=False,cells=8192,points=16640,faces=32896,internal_faces=16256,physical_faces=256,empty_faces=16384,maximum_nonorthogonality_degrees=nonorth,maximum_skewness=skew,min_positive_volume_m3=min(V),total_volume_m3=sum(V),full_source_cell_point_face_vertex_bijection=True,wrapped_one_hole_source_zones_proved=True,native_mesh_sha256={n:sha(poly/n) for n in ('points','faces','owner','neighbour','boundary')},actual_native_C_Vc_sha256={n:sha(case/'0'/n) for n in ('C','Vc')},initial_state_prescribed_not_solved_conservative=True,actual_initial_cell_rest=True,actual_rotating_wall_max_error_m_s=wall_error,actual_rotating_wall_max_normal_velocity_m_s=wall_normal,actual_native_force_registration=probe,actual_native_registration_process_proof=reg,primary_initial_manifest_sha256=sha(case/'provenance/native-initial-primary-sha256.json'),derived_initial_manifest_sha256=sha(case/'provenance/native-initial-derived-sha256.json'));json_new(case/'provenance/preflight-validation.json',receipt);return receipt

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,required=True);g=ap.add_mutually_exclusive_group(required=True);g.add_argument('--mesh',type=Path);g.add_argument('--validate-native',action='store_true');g.add_argument('--verify-inputs',action='store_true');a=ap.parse_args();result=prepare(a.mesh,a.case) if a.mesh else validate_native(a.case) if a.validate_native else dict(verified_files=len(verify_inputs(a.case)));print(json.dumps(result,indent=2))
if __name__=='__main__':main()