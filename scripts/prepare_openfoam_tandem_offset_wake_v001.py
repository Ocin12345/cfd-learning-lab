"""Prepare/check the frozen two-hole transient pilot; never invokes native programs or CFD."""
from pathlib import Path
import argparse,json,re,math
import tandem_offset_openfoam_common_v001 as common
from tandem_offset_openfoam_common_v001 import MeshSpec,require,sha,write_new,json_new,copied,header,field,blocks,top,kind,values,counted,native_poly,patch_table,native_geometry,source_mapping
COMMON_SHA='5e1129d55bda3cc681f94e6fffc4d982292d5b1b3a093f55a4fa74dd85b08d1c'
TARGET_SHA='426448b7e6909e75dbf08f065d93f90699b98f12be007891631d9ed9b10f71d7'
PLAN_SHA='79f894cdddbf65209b91b6cfd0121909eb69428e9dcf101ae267c0112ec5e055'
MESH_SHA='f5db9ca48d6e5bcc27d387febbc36fe1e2239f1deac62061362ce902cfe0583d';NPZ_SHA='099165a86f06dd6566c6d76860f5a0935fba4b956e16d6908f6f1b735270b178';SETUP_SHA='263f87d0068959e316215e85c6e1391627bfd5d6d87e68b5d4c6b13fa0324bb4'
RUN='2026-11-02-tandem-offset-wake-v001';PATCH_COUNTS=dict(inlet=60,outlet=60,top=160,bottom=160,body1=32,body2=32);ZONES=dict(inlet=4,outlet=5,top=6,bottom=7,body1=8,body2=9);SPEC=MeshSpec(9472,19446,38140,18692,PATCH_COUNTS)
INSTALLED_SOURCES=(
'applications/modules/incompressibleFluid/incompressibleFluid.C','applications/modules/incompressibleFluid/correctPressure.C',
'src/finiteVolume/fields/fvPatchFields/constraint/symmetry/symmetryFvPatchField.C',
'src/finiteVolume/cfdTools/general/solutionControl/pimpleControl/pimpleLoop/pimpleLoop.C','src/finiteVolume/cfdTools/general/solutionControl/pimpleControl/pimpleNoLoopControl/pimpleNoLoopControl.C','src/finiteVolume/cfdTools/general/solutionControl/pisoControl/pisoControl.C','src/finiteVolume/fvMatrices/fvMatrix/fvMatrixSolve.C',
'src/functionObjects/forces/forces/forces.C','src/functionObjects/forces/forcesBase/forcesBase.C','src/functionObjects/forces/forceCoeffs/forceCoeffs.C','src/functionObjects/field/wallShearStress/wallShearStress.C',
'src/OpenFOAM/db/functionObjects/functionObject/functionObject.C','src/OpenFOAM/db/functionObjects/functionObjectList/functionObjectList.C','src/functionObjects/field/fieldExpression/fieldExpression.C',
'src/functionObjects/utilities/writeObjects/writeObjects.C','src/functionObjects/field/writeCellCentres/writeCellCentres.C','src/functionObjects/field/writeCellVolumes/writeCellVolumes.C','src/functionObjects/field/vorticity/vorticity.C','src/finiteVolume/finiteVolume/fvc/fvcCurl.C','applications/utilities/postProcessing/foamPostProcess/foamPostProcess.C')
def root():return Path(__file__).resolve().parents[3]
def run_root():return root()/'native/runs'/RUN
def private():return root()/'local/tandem-offset-wake-v001'
def checked_case(p):
 p=Path(p).resolve();require(p==(run_root()/'openfoam-pilot-v001').resolve(),'Dedicated frozen tandem pilot only');return p

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
nativeBodyShear { type wallShearStress; libs ("libfieldFunctionObjects.so"); patches (body1 body2); log false; executeControl timeStep; executeInterval 1; writeControl timeStep; writeInterval 1; }
'''
def functions(probe=False):
 s=''
 for body,centre in [('body1',(.5,-.025)),('body2',(.9,.025))]:
  s+=force_function(('probe' if probe else '')+body[0].upper()+body[1:]+'Forces' if probe else body+'Forces','forces',body,centre)
  s+=force_function(('probe' if probe else '')+body[0].upper()+body[1:]+'Coefficients' if probe else body+'Coefficients','forceCoeffs',body,centre)
 return s+field_functions()
def dictionaries():
 d={'system/controlDict':'''solver incompressibleFluid;
startFrom startTime; startTime 0; stopAt endTime; endTime 0.25;
deltaT 0.05; adjustTimeStep false;
writeControl timeStep; writeInterval 1; purgeWrite 0;
writeFormat ascii; writePrecision 14; writeCompression off;
timeFormat general; timePrecision 10; runTimeModifiable false;
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
PIMPLE { nOuterCorrectors 30; nCorrectors 2; nNonOrthogonalCorrectors 0; momentumPredictor yes; consistent yes; }
''','system/decomposeParDict':'numberOfSubdomains 2;\nmethod simple;\nsimpleCoeffs { n (2 1 1); delta 0.001; }\n','constant/physicalProperties':'viscosityModel constant;\nnu 0.0001;\n','constant/momentumTransport':'simulationType laminar;\n','0/U':'''dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
 inlet { type fixedValue; value uniform (0.1 0 0); }
 outlet { type zeroGradient; }
 top { type symmetry; }
 bottom { type symmetry; }
 body1 { type noSlip; }
 body2 { type noSlip; }
 frontAndBackPlanes { type empty; }
}
''','0/p':'''dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
 inlet { type zeroGradient; }
 outlet { type fixedValue; value uniform 0; }
 top { type symmetry; }
 bottom { type symmetry; }
 body1 { type zeroGradient; }
 body2 { type zeroGradient; }
 frontAndBackPlanes { type empty; }
}
''','system/functions':functions(),'system/initialDiagnostics':'''nativeCentres { type writeCellCentres; libs ("libfieldFunctionObjects.so"); }
nativeVolumes { type writeCellVolumes; libs ("libfieldFunctionObjects.so"); }
nativeInitialFields { type writeObjects; libs ("libutilityFunctionObjects.so"); objects (U p phi); writeOption anyWrite; }
'''+field_functions(),'system/forcePreflight':'''nativeInitialFields { type writeObjects; libs ("libutilityFunctionObjects.so"); objects (U p phi); writeOption anyWrite; }
'''+functions(True)}
 return {n:header(Path(n).name,'volVectorField' if n=='0/U' else 'volScalarField' if n=='0/p' else 'dictionary')+v for n,v in d.items()}

def verify_inputs(case):
 case=checked_case(case);require(sha(Path(common.__file__))==COMMON_SHA,'Frozen common reader changed');m=json.loads((case/'provenance/prepared-input-sha256.json').read_text());require(all(sha(case/n)==h for n,h in m.items()),'Prepared source/control/physics input changed');return m

def prepare(mesh,case):
 import numpy as np
 case=checked_case(case);mesh=Path(mesh).resolve();require(not case.exists() and mesh==(run_root()/'mesh/tandem-offset-wake.msh').resolve(),'Fresh dedicated case and selected shared mesh');require(sha(mesh)==MESH_SHA and sha(mesh.parent/'geometry.npz')==NPZ_SHA and sha(mesh.parent/'setup.json')==SETUP_SHA,'Actual accepted shared input hashes');require(sha(private()/'TARGETS.md')==TARGET_SHA and sha(private()/'EXECUTION-PLAN-v001.md')==PLAN_SHA and sha(Path(common.__file__))==COMMON_SHA,'Frozen target/plan/common source');setup=json.loads((mesh.parent/'setup.json').read_text());require(setup['cells']==9472 and setup['density_kg_m3']==1 and setup['kinematic_viscosity_m2_s']==.0001 and setup['inlet_velocity_m_s']==[.1,0] and setup['time_step_s']==.05 and setup['time_steps']==800 and setup['pilot_time_steps']==5,'Actual physical setup')
 case.mkdir(parents=True)
 for n,t in dictionaries().items():write_new(case/n,t)
 for n in ('U','p'):write_new(case/'provenance/prescribed-initial-fields'/n,dictionaries()['0/'+n])
 copied(mesh,case/'input-tandem-offset-wake.msh')
 for p in (mesh.parent/'geometry.npz',mesh.parent/'setup.json',private()/'TARGETS.md',private()/'EXECUTION-PLAN-v001.md',Path(__file__),Path(common.__file__),private()/'run-openfoam-v001.sh',private()/'Launch-OpenFOAM-v001.ps1',private()/'OPENFOAM-PREREQUISITES-v001.md'):copied(p,case/'provenance'/p.name)
 with np.load(mesh.parent/'geometry.npz',allow_pickle=False) as f:src={n:np.asarray(f[n]).tolist() for n in ('nodes','quads','centres','areas','grid_shape','ix','iy','face_nodes','face_owner','face_neighbour','face_zone','face_area_vectors')}
 require(src['grid_shape']==[60,160] and len(src['quads'])==9472 and len(src['nodes'])==9723 and len(src['face_nodes'])==19196,'Accepted sparse source shape');slots=set(zip(src['ix'],src['iy']));expected={(i,j) for j in range(60) for i in range(160) if not(36<=i<44 and 24<=j<32) and not(68<=i<76 and 28<=j<36)};require(slots==expected,'Actual two-hole sparse mask');json_new(case/'provenance/source-geometry.json',src);json_new(case/'provenance/installed-source-list.json',list(INSTALLED_SOURCES));json_new(case/'provenance/official-primary-source-urls.json',{n:'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/'+n for n in INSTALLED_SOURCES})
 json_new(case/'declared-comparison.json',dict(phase='five_step_transient_pilot',initial_cell_U=[0,0,0],initial_cell_p=0,initial_state_kind='prescribed rest with active inlet; native initial phi retained, not solved-conservative',solved_gate_first_index=1,pilot_steps=5,dt_s=.05,pilot_final_s=.25,total_steps_prospective=800,additional_steps_prospective=795,final_s_prospective=40,outer_correctors=30,pressure_correctors=2,force_body_names=['body1','body2'],force_reference_rho=1,force_reference_U=.1,force_reference_D=.1,force_reference_A_m2=.001,force_denominator_N=5e-6,force_coefficient_relative_absolute_target=1e-6,force_norm='abs(delta coefficient)<=1e-6*max(1,abs(native coefficient))',primary_restart_inputs=['U','p','phi','time','mesh','physics','sources'],derived_restart_outputs=['vorticity','wallShearStress','force function-object startup outputs'],target_sha256=TARGET_SHA,plan_sha256=PLAN_SHA,mesh_sha256=MESH_SHA,geometry_sha256=NPZ_SHA,setup_sha256=SETUP_SHA,common_source_sha256=COMMON_SHA,physical_validation=False,production_admitted=False))
 json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file() and p.relative_to(case).as_posix() not in ('0/U','0/p')});json_new(case/'run-status.json',dict(stage='prepared_only',solver_started=False,production_admitted=False));verify_inputs(case);return dict(case=str(case),prepared=True,solver_started=False,prepared_input_manifest_sha256=sha(case/'provenance/prepared-input-sha256.json'))

def validate_native(case):
 case=checked_case(case);verify_inputs(case);require(not (case/'solver.log').exists(),'Preflight before CFD');reg=json.loads((case/'provenance/native-registration-process-status.json').read_text())
 for label,ranks in [('serial',1),('parallel',2)]:
  r=reg[label];text=(case/r['log']).read_text();require(r['return_code']==0 and r['mpi_ranks']==ranks and r['log_sha256']==sha(case/r['log']) and re.search(r'^End\s*$',text,re.M) and 'FOAM FATAL' not in text,'Actual completed registration process/log');times=re.findall(r'^Time = ('+common.FLOAT+r')s?\s*$',text,re.M);require(times and all(float(t)==0 for t in times),'Actual native time0 only in preflight')
  if ranks==2:require(re.search(r'\bnProcs\s*:\s*2\b',text),'Actual parallel registration')
 text=(case/'mesh-check.log').read_text();require('Mesh OK.' in text and 'checkMesh -allGeometry -allTopology' in text,'Actual full checkMesh');nonorth=float(re.search(r'Mesh non-orthogonality Max:\s*('+common.FLOAT+')',text)[1]);skew=float(re.search(r'Max skewness\s*=\s*('+common.FLOAT+')',text)[1]);require(nonorth<=75 and skew<=4,'Native quality targets');poly=native_poly(case);pa=patch_table(poly);require(set(pa)==set(PATCH_COUNTS)|{SPEC.empty},'Exact native physical patch set')
 for n,q in pa.items():require(q['count']==(18944 if n==SPEC.empty else PATCH_COUNTS[n]) and q['kind']==('empty' if n==SPEC.empty else 'symmetry' if n in ('top','bottom') else 'wall' if n.startswith('body') else 'patch'),'Actual native patch type/count '+n)
 g=native_geometry(poly,SPEC);C=field((case/'0/C').read_text(),9472,3,[0,1,0,0,0,0,0]);V=field((case/'0/Vc').read_text(),9472,1,[0,3,0,0,0,0,0]);src=json.loads((case/'provenance/source-geometry.json').read_text());mapping=source_mapping(g,src,C,V,SPEC,pa,ZONES);U=field((case/'0/U').read_text(),9472,3,[0,1,-1,0,0,0,0]);P=field((case/'0/p').read_text(),9472,1,[0,2,-2,0,0,0,0]);require(all(math.isfinite(v) and abs(v)<=1e-12 for row in U for v in row) and all(math.isfinite(v) and abs(v)<=1e-12 for v in P),'Actual prescribed initial cell rest/p0');ub=blocks((case/'0/U').read_text());pb=blocks((case/'0/p').read_text());require(kind(ub['inlet'])=='fixedValue' and all(math.dist(v,[.1,0,0])<=1e-12 for v in values(top(ub['inlet']),60,3,'value')) and kind(pb['inlet'])=='zeroGradient','Actual initial inlet types/value');require(kind(ub['outlet'])=='zeroGradient' and kind(pb['outlet'])=='fixedValue' and max(map(abs,values(top(pb['outlet']),60,1,'value')))<=1e-12,'Actual initial outlet p0')
 for n in ('body1','body2'):
  require(kind(ub[n])=='noSlip' and kind(pb[n])=='zeroGradient','Actual noSlip/zeroGradient bodies')
  if re.search(r'\bvalue\s+',top(ub[n])):require(max(abs(v) for row in values(top(ub[n]),32,3,'value') for v in row)<=1e-12,'Optional noSlip numerical value')
 for n in ('top','bottom'):require(kind(ub[n])==kind(pb[n])=='symmetry','Exact native symmetry fields')
 cells=[]
 for r in (0,1):
  n,b=counted(case/f'processor{r}/0/polyMesh/cellProcAddressing');a=list(map(int,b.split()));require(n==len(a),'Rank address count');cells.extend(a)
 require(sorted(cells)==list(range(9472)),'Actual MPI cell cover disjoint/complete')
 probe={}
 for name in ('Body1','Body2'):
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
     baseline='provenance/native-initial-derived-baseline/'+rel;copied(p,case/baseline);nc=9472 if folder==case/'0' else counted(folder/'polyMesh/cellProcAddressing')[0];derived[rel]=dict(sha256=sha(p),baseline_path=baseline,cells=nc,dimensions=[0,0,-1,0,0,0,0] if p.name=='vorticity' else [0,2,-2,0,0,0,0])
    else:primary[rel]=sha(p)
 require(all((case/f'processor{r}/0'/n).is_file() for r in (0,1) for n in ('U','p','phi','vorticity','wallShearStress')),'Actual initial rank fields');json_new(case/'provenance/native-mapping.json',mapping);json_new(case/'provenance/native-initial-primary-sha256.json',primary);json_new(case/'provenance/native-initial-derived-sha256.json',derived)
 receipt=dict(passed=True,solver_started=False,cells=9472,points=19446,faces=38140,internal_faces=18692,physical_faces=504,empty_faces=18944,maximum_nonorthogonality_degrees=nonorth,maximum_skewness=skew,min_positive_volume_m3=min(V),total_volume_m3=sum(V),full_source_cell_point_face_vertex_bijection=True,two_holes_source_zones_proved=True,native_mesh_sha256={n:sha(poly/n) for n in ('points','faces','owner','neighbour','boundary')},actual_native_C_Vc_sha256={n:sha(case/'0'/n) for n in ('C','Vc')},initial_state_prescribed_not_solved_conservative=True,actual_initial_cell_rest=True,actual_native_force_registration=probe,actual_native_registration_process_proof=reg,primary_initial_manifest_sha256=sha(case/'provenance/native-initial-primary-sha256.json'),derived_initial_manifest_sha256=sha(case/'provenance/native-initial-derived-sha256.json'));json_new(case/'provenance/preflight-validation.json',receipt);return receipt

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,required=True);g=ap.add_mutually_exclusive_group(required=True);g.add_argument('--mesh',type=Path);g.add_argument('--validate-native',action='store_true');g.add_argument('--verify-inputs',action='store_true');a=ap.parse_args();result=prepare(a.mesh,a.case) if a.mesh else validate_native(a.case) if a.validate_native else dict(verified_files=len(verify_inputs(a.case)));print(json.dumps(result,indent=2))
if __name__=='__main__':main()