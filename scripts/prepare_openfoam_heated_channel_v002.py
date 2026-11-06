"""Prepare the heated channel flow and separately admitted one-step frozen-flow thermal pilot; no solver launch."""
from pathlib import Path
import argparse,hashlib,json,math,re,shutil
RUN='2026-11-02-heated-channel-v001'
PATCHES=('inlet','outlet','top','bottom')
EMPTY='frontAndBackPlanes'
CELLS=2304
POINTS=4850
INTERNAL_FACES=4488
TOTAL_FACES=9336
TARGET_SHA='0a163061a9836c096d8b30b16c81410285e29bf3507f10d231619c43b64cab02'
FLOAT=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
INSTALLED_SOURCES=(
 'applications/modules/incompressibleFluid/incompressibleFluid.C',
 'applications/modules/functions/functions.C',
 'src/functionObjects/solvers/scalarTransport/scalarTransport.C',
 'src/functionObjects/solvers/scalarTransport/scalarTransport.H',
 'src/functionObjects/field/grad/grad.C',
 'src/functionObjects/field/grad/gradTemplates.C',
 'src/functionObjects/field/fieldExpression/fieldExpression.C',
 'src/finiteVolume/finiteVolume/gradSchemes/gaussGrad/gaussGrad.C',
 'src/functionObjects/utilities/writeObjects/writeObjects.C',
 'src/functionObjects/utilities/writeObjects/writeObjects.H',
 'src/OpenFOAM/db/regIOobject/regIOobjectWrite.C',
 'src/functionObjects/field/writeCellCentres/writeCellCentres.C',
 'src/functionObjects/field/writeCellVolumes/writeCellVolumes.C',
 'src/OpenFOAM/db/functionObjects/functionObjectList/functionObjectList.C',
 'applications/utilities/postProcessing/foamPostProcess/foamPostProcess.C',
 'applications/utilities/parallelProcessing/decomposePar/decomposePar.C',
 'src/OpenFOAM/db/Time/timeSelector.C',
 'tutorials/incompressibleFluid/pitzDailyScalarTransport/system/controlDict',
 'tutorials/incompressibleFluid/pitzDailyScalarTransport/system/functions',
 'etc/caseDicts/functions/solvers/scalarTransport.cfg',
)
def require(v,m):
 if not v:raise ValueError(m)
def lab_root():return Path(__file__).resolve().parents[3]
def run_root():return lab_root()/'native/runs'/RUN
def private_root():return lab_root()/'local/heated-channel-v001'
def checked_case(path,phase='flow'):
 p=Path(path).resolve();expected=run_root()/('openfoam-flow-v002' if phase=='flow' else 'openfoam-thermal-pilot-v002')
 require(p==expected.resolve(),'Only dedicated heated-channel '+phase+' case allowed');return p
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write_new(path,text):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('x',encoding='ascii',newline='\n') as f:f.write(text.rstrip()+'\n')
def json_new(path,value):write_new(path,json.dumps(value,indent=2,allow_nan=False))
def copied(src,dst):
 require(not dst.exists(),'Refuse overwrite '+str(dst));dst.parent.mkdir(parents=True,exist_ok=True);before=sha(src);shutil.copy2(src,dst);require(sha(src)==sha(dst)==before,'Copy hash differs')
def header(name,kind='dictionary'):return 'FoamFile\n{\n version 2.0;\n format ascii;\n class '+kind+';\n object '+name+';\n}\n\n'
def thermal_T():
 return header('T','volScalarField')+"""dimensions [0 0 0 1 0 0 0];
internalField uniform 300;
boundaryField
{
 inlet { type fixedValue; value uniform 300; }
 outlet { type zeroGradient; }
 top { type fixedValue; value uniform 350; }
 bottom { type fixedValue; value uniform 350; }
 frontAndBackPlanes { type empty; }
}
"""
def grad_function():return """temperatureGradient
{
 type grad;
 libs ("libfieldFunctionObjects.so");
 field T;
 result gradT;
 executeControl timeStep;
 executeInterval 1;
 writeControl timeStep;
 writeInterval 1;
}
"""
def frozen_function():return """frozenRuntimeFields
{
 type writeObjects;
 libs ("libutilityFunctionObjects.so");
 objects (U p phi);
 writeOption anyWrite;
 executeControl timeStep;
 executeInterval 1;
 writeControl timeStep;
 writeInterval 1;
}
"""
def dictionaries(phase='flow'):
 content={
 'system/controlDict':"""solver incompressibleFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 2000;
deltaT 1;
writeControl timeStep;
writeInterval 100;
purgeWrite 0;
writeFormat ascii;
writePrecision 14;
writeCompression off;
timeFormat general;
timePrecision 10;
runTimeModifiable false;
""",
 'system/fvSchemes':"""ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U); div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""",
 'system/fvSolution':"""solvers
{
 p { solver GAMG; tolerance 1e-10; relTol 0.01; smoother GaussSeidel; }
 U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01; }
}
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; }
relaxationFactors { equations { U 0.7; } }
""",
 'system/decomposeParDict':'numberOfSubdomains 2;\nmethod simple;\nsimpleCoeffs { n (2 1 1); delta 0.001; }\n',
 'constant/physicalProperties':'viscosityModel constant;\nnu 0.01;\n',
 'constant/momentumTransport':'simulationType laminar;\n',
 '0/U':"""dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
 inlet { type fixedValue; value uniform (1 0 0); }
 outlet { type zeroGradient; }
 top { type noSlip; }
 bottom { type noSlip; }
 frontAndBackPlanes { type empty; }
}
""",
 '0/p':"""dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
 inlet { type zeroGradient; }
 outlet { type fixedValue; value uniform 0; }
 top { type zeroGradient; }
 bottom { type zeroGradient; }
 frontAndBackPlanes { type empty; }
}
""",
 'system/initialDiagnostics':grad_function(),
 }
 if phase=='pilot':
  content.pop('0/U');content.pop('0/p');content.pop('system/initialDiagnostics')
  content['system/controlDict']="""solver functions;
subSolver incompressibleFluid;
subSolverTime 2000;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 0.005;
deltaT 0.005;
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
"""
  content['system/fvSchemes']="""ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,T) Gauss upwind; div(phi,U) bounded Gauss linearUpwind grad(U); div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear orthogonal; }
interpolationSchemes { default linear; }
snGradSchemes { default orthogonal; }
"""
  content['system/fvSolution']="""solvers
{
 T { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0; nCorrectors 0; }
 p { solver GAMG; tolerance 1e-10; relTol 0.01; smoother GaussSeidel; }
 U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01; }
}
PIMPLE { nOuterCorrectors 1; nCorrectors 0; nNonOrthogonalCorrectors 0; }
relaxationFactors { equations { T 1; } }
"""
  content['system/functions']="""transportTemperature
{
 type scalarTransport;
 libs ("libsolverFunctionObjects.so");
 field T;
 phi phi;
 diffusivity constant;
 D 0.01;
 executeControl timeStep;
 executeInterval 1;
 writeControl timeStep;
 writeInterval 1;
}
"""+grad_function()+frozen_function()
  content['system/baselineDiagnostics']=grad_function()+frozen_function()
 result={n:header(Path(n).name,'volVectorField' if n=='0/U' else 'volScalarField' if n=='0/p' else 'dictionary')+v for n,v in content.items()}
 result['0/T']=thermal_T()
 return result
def clean(text):
    return re.sub(r'/\*.*?\*/|//[^\n]*', '', text, flags=re.S)
def counted(path):
    text = clean(Path(path).read_text(encoding='ascii'))
    h = re.search(r'FoamFile\s*\{([^{}]*)\}', text)
    require(h is not None and re.search(r'\bformat\s+ascii\s*;',h[1]), 'Require ASCII native list')
    m = re.fullmatch(r'\s*(\d+)\s*\((.*)\)\s*', text[h.end():], re.S)
    require(m is not None, 'Invalid native counted list '+str(path))
    return int(m[1]), m[2]
def native_poly(case):
    candidates = [case/name for name in ('0/polyMesh','constant/polyMesh') if (case/name/'points').is_file()]
    require(len(candidates)==1, 'Exactly one initial native polyMesh required')
    return candidates[0]
def native_patches(poly):
    n,b=counted(poly/'boundary');rows={}
    for name,block in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',b):
        rows[name.strip('"')]=dict(kind=re.search(r'\btype\s+(\w+)\s*;',block)[1],count=int(re.search(r'\bnFaces\s+(\d+)\s*;',block)[1]),start=int(re.search(r'\bstartFace\s+(\d+)\s*;',block)[1]))
    require(n==5 and set(rows)==set(PATCHES)|{EMPTY},'Native square-wake patch set differs')
    require(all(rows[n]['kind']=='patch' and rows[n]['count']==24 for n in ('inlet','outlet')),'Native opening counts/types differ')
    require(all(rows[n]['kind']=='wall' and rows[n]['count']==96 for n in ('top','bottom')),'Native stationary wall counts/types differ')
    require(rows[EMPTY]['kind']=='empty' and rows[EMPTY]['count']==4608,'Native empty face count differs')
    ranges=sorted((v['start'],v['start']+v['count']) for v in rows.values())
    require(ranges[0][0]==INTERNAL_FACES and ranges[-1][1]==TOTAL_FACES and all(a[1]==b[0] for a,b in zip(ranges,ranges[1:])),'Native boundary ranges differ')
    return rows
def mesh_quality(log):
    require('checkMesh -allGeometry -allTopology' in log and 'Mesh OK.' in log and not re.search(r'Failed\s+\d+\s+mesh checks|FOAM FATAL',log), 'Full native checkMesh must pass')
    patterns = {'cells':r'\bcells:\s+(\d+)', 'minimum_volume_m3':r'Min volume\s*=\s*('+FLOAT+')',
        'maximum_nonorthogonality_degrees':r'Mesh non-orthogonality Max:\s*('+FLOAT+')', 'maximum_skewness':r'Max skewness\s*=\s*('+FLOAT+')'}
    result = {}
    for name, pattern in patterns.items():
        m=re.search(pattern,log); require(m is not None,'Missing quality '+name)
        result[name]=int(m[1]) if name=='cells' else float(m[1])
    require(result['cells']==CELLS and result['minimum_volume_m3']>0 and result['maximum_nonorthogonality_degrees']<=75 and result['maximum_skewness']<=4,'Mesh quality gate missed')
    return result
def native_geometry(poly):
    n, body = counted(poly/'points')
    points = [tuple(map(float,v.split())) for v in re.findall(r'\(([^()]*)\)',body)]
    require(n==POINTS==len(points) and all(len(p)==3 and all(map(math.isfinite,p)) for p in points),'Native point count/finite gate')
    bounds = [[min(p[j] for p in points),max(p[j] for p in points)] for j in range(3)]
    require(all(abs(a-b)<=1e-10 for row,target in zip(bounds,([0,1],[0,.1],[-.005,.005])) for a,b in zip(row,target)), 'Native heated channel/extrusion bounds differ')
    n,body=counted(poly/'faces')
    faces=[]
    for size,labels in re.findall(r'(\d+)\s*\(([^()]*)\)',body):
        face=list(map(int,labels.split()))
        require(len(face)==int(size)==4 and len(set(face))==4 and all(0<=i<len(points) for i in face),'Invalid native quadrilateral')
        faces.append(face)
    require(n==TOTAL_FACES==len(faces),'Native face count differs')
    def labels(name,expected):
        n,body=counted(poly/name); row=list(map(int,body.split()))
        require(n==expected==len(row) and all(0<=i<CELLS for i in row),'Native '+name+' count/range differs')
        return row
    owner=labels('owner',TOTAL_FACES); neighbour=labels('neighbour',INTERNAL_FACES)
    require(set(owner)==set(range(CELLS)) and all(owner[i]!=neighbour[i] for i in range(INTERNAL_FACES)), 'Native cell owner/neighbour coverage differs')
    volumes=[0.0]*CELLS
    for i,face in enumerate(faces):
        a=points[face[0]]; face_volume=0.0
        for j in range(1,len(face)-1):
            b,c=points[face[j]],points[face[j+1]]
            face_volume+=(a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0]))/6
        volumes[owner[i]]+=face_volume
        if i<len(neighbour):
            volumes[neighbour[i]]-=face_volume
    expected=(1/96)*(.1/24)*.01
    require(all(math.isfinite(v) and v>0 and abs(v-expected)<=expected*1e-8 for v in volumes),'Actual oriented native cell volumes differ from uniform extrusion')
    return dict(points=points,faces=faces,owner=owner,neighbour=neighbour,volumes=volumes,bounds=bounds)
def native_field_values(path,count,components,dimensions):
    text=clean(path.read_text(encoding='ascii'))
    require(re.search(r'\bformat\s+ascii\s*;',text),'ASCII native field required')
    d=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',text)
    require(d is not None and list(map(int,d[1].split()))==dimensions,'Native field dimensions differ')
    m=re.search(r'\binternalField\s+uniform\s+([^;]+);',text)
    if m:
        row=list(map(float,m[1].replace('(',' ').replace(')',' ').split()))
        require(len(row)==components and all(map(math.isfinite,row)),'Invalid uniform native field')
        return [row[:] for _ in range(count)] if components>1 else [row[0]]*count
    kind='vector' if components>1 else 'scalar'
    m=re.search(r'\binternalField\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
    require(m is not None and int(m[1])==count,'Native field count differs')
    data=list(map(float,m[2].replace('(',' ').replace(')',' ').split()))
    require(len(data)==count*components and all(map(math.isfinite,data)),'Invalid native field values')
    return [data[i:i+components] for i in range(0,len(data),components)] if components>1 else data

def verify_prepared(case):
    manifest=json.loads((case/'provenance/prepared-input-sha256.json').read_text())
    require(all(sha(case/name)==digest for name,digest in manifest.items()),'Prepared inputs changed')
    return manifest

def copied_provenance(case,extra=()):
    files=[Path(__file__),private_root()/'run-openfoam-v002.sh',private_root()/'Launch-OpenFOAM-v002.ps1',private_root()/'DESIGN-PROSPECTIVE.md',private_root()/'TARGETS.md',*extra]
    require(sha(private_root()/'TARGETS.md')==TARGET_SHA,'Frozen target document changed')
    for source in files:copied(source,case/'provenance'/source.name)

def manifest_new(case):
    json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file() and not(case.name=='openfoam-thermal-pilot-v002' and p.relative_to(case).as_posix() in ('0/U','0/p','0/phi'))})

def prepare(mesh,case):
    case=checked_case(case);mesh=Path(mesh).resolve()
    require(not case.exists(),'Fresh native flow case required')
    require(mesh==run_root()/'mesh/heated-channel.msh','Only selected heated-channel source mesh')
    raw=mesh.read_bytes();require(re.search(rb'\(2\s+2\)',raw),'Require2D FluentASCII')
    m=re.search(rb'^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)',raw,re.M);require(m is not None and int(m[1],16)==CELLS,'Require2304sourcecells')
    for name in ('fluid','interior',*PATCHES):require(re.search(rb'\(45\s+\([^)]*\s'+name.encode()+rb'\)\(\)\)',raw),'Missing named source zone '+name)
    setup=json.loads((mesh.parent/'setup.json').read_text(encoding='utf-8-sig'))
    require(setup['mesh_sha256']==sha(mesh) and setup['geometry_npz_sha256']==sha(mesh.parent/'geometry.npz'),'Source geometry hashes differ')
    expected=dict(cells=CELLS,nx=96,ny=24,length_m=1,height_m=.1,density_kg_m3=1,dynamic_viscosity_Pa_s=.01,kinematic_viscosity_m2_s=.01,
      specific_heat_J_kg_K=1000,thermal_conductivity_W_m_K=10,thermal_diffusivity_m2_s=.01,reference_depth_m=1,openfoam_extrusion_depth_m=.01,
      initial_pressure_Pa=0,outlet_gauge_pressure_Pa=0,initial_temperature_K=300,inlet_temperature_K=300,wall_temperature_K=350,time_step_s=.005,time_steps=600,final_time_s=3)
    require(all(setup[n]==v for n,v in expected.items()),'Shared material/geometry/time setup differs')
    require(setup['initial_velocity_m_s']==[0,0] and setup['inlet_velocity_m_s']==[1,0],'Shared initial/inlet U differs')
    require(shutil.disk_usage(case.anchor).free>=10*1024**3,'Require10GiB disk')
    case.mkdir(parents=True)
    for name,text in dictionaries().items():write_new(case/name,text)
    copied(mesh,case/'input-heated-channel.msh');copied_provenance(case,[mesh.parent/'setup.json',mesh.parent/'geometry.npz'])
    declaration=dict(model='Steady laminar channel followed by separately admitted passive temperature on actual frozen computed flow',physical_setup=setup,
      phase='steady_flow',initial_U_m_s=[0,0,0],initial_p=0,inlet_U_m_s=[1,0,0],depth_m=.01,
      solver='OpenFOAM Foundation14 incompressibleFluid',steady_iterations=2000,saved_iterations=list(range(100,2001,100)),selected_iterations=[500,1000,1500,1800,2000],
      max_printed_initial_p_Ux_Uy_residual=1e-6,global_relative_mass_imbalance_target=1e-5,local_divergence_scale_Dh_over_U_s=.2,local_divergence_max_target=1e-5,local_divergence_volume_weighted_RMS_target=1e-6,
      U_stability_iterations=[1800,2000],U_relative_L2_target=1e-3,mesh_source_sha256=sha(mesh),geometry_sha256=sha(mesh.parent/'geometry.npz'),target_document_sha256=TARGET_SHA,
      thermal_prospective_steps=600,thermal_dt_s=.005,thermal_production_admitted=False,physical_validation=False,mesh_independence=False,time_step_independence=False)
    json_new(case/'declared-comparison.json',declaration);manifest_new(case);json_new(case/'run-status.json',dict(stage='prepared_only',solver_started=False))
    return dict(case=str(case),cells=CELLS,solver_started=False,mesh_sha256=sha(mesh))

def validate_native(case):
    case=checked_case(case);require(not (case/'solver.log').exists() and not list(case.glob('processor*')),'Native preflight before solve/decomposition only');verify_prepared(case)
    poly=native_poly(case);patches=native_patches(poly);geometry=native_geometry(poly);quality=mesh_quality((case/'mesh-check.log').read_text())
    C=native_field_values(case/'0/C',CELLS,3,[0,1,0,0,0,0,0]);V=native_field_values(case/'0/Vc',CELLS,1,[0,3,0,0,0,0,0])
    require(all(v>0 and abs(v-w)<=1e-14 for v,w in zip(V,geometry['volumes'])),'Native Vc/oriented volume differs')
    keys=[]
    for x,y,z in C:
        ix=round(x*96-.5);iy=round(y*240-.5)
        require(0<=ix<96 and 0<=iy<24 and abs(x-(ix+.5)/96)<=1e-10 and abs(y-(iy+.5)/240)<=1e-10 and abs(z)<=1e-10,'Native C off source rectangle lattice')
        keys.append(iy*96+ix)
    require(sorted(keys)==list(range(CELLS)),'Native C lattice not bijective')
    native_field_values(case/'0/gradT',CELLS,3,[0,-1,0,1,0,0,0])
    receipt=dict(passed=True,solver_started=False,cells=CELLS,points=POINTS,faces=TOTAL_FACES,internal_faces=INTERNAL_FACES,actual_native_patches=patches,
      mesh_quality=quality,bounds_m=geometry['bounds'],actual_oriented_volume_min_m3=min(geometry['volumes']),actual_oriented_volume_max_m3=max(geometry['volumes']),total_volume_m3=sum(geometry['volumes']),
      native_C_grid_bijective=True,native_T_gradient_function_registered=True,scalarTransport_runtime_registration_deferred_to_separately_admitted_single_step=True,
      mesh_location=poly.relative_to(case).as_posix(),native_mesh_sha256={n:sha(poly/n) for n in ('points','faces','owner','neighbour','boundary')},native_initial_C_Vc_sha256={n:sha(case/'0'/n) for n in ('C','Vc')})
    json_new(case/'provenance/preflight-validation.json',receipt);return receipt

def prepare_pilot(case):
    case=checked_case(case,'pilot');flow=checked_case(run_root()/'openfoam-flow-v002');require(not case.exists(),'Fresh one-step pilot case required')
    status=json.loads((flow/'run-status.json').read_text());require(status.get('stage')=='finished_flow_execution','Native flow genuinely complete required')
    audit=private_root()/'openfoam-flow-checks-v001';checks=json.loads((audit/'checks.json').read_text());require(checks.get('all_checks_passed') is True and not checks.get('failed_checks'),'Flow audit must pass before thermal pilot')
    require(checks['preparation_source_sha256']==sha(Path(__file__)),'Flow audit preparation hash differs')
    verify_prepared(flow);pre=json.loads((flow/'provenance/preflight-validation.json').read_text());poly=native_poly(flow)
    require(all(sha(poly/n)==digest for n,digest in pre['native_mesh_sha256'].items()),'Frozen flow mesh changed')
    case.mkdir(parents=True)
    for name,text in dictionaries('pilot').items():write_new(case/name,text)
    for name in ('points','faces','owner','neighbour','boundary'):copied(poly/name,case/'0/polyMesh'/name)
    for name in ('U','p','phi','C','Vc'):copied(flow/'2000'/name,case/'2000'/name)
    # Actual baseline files are exact native flow copies, read/written by non-solver diagnostics.
    # Immutable provenance copies retain byte binding; native runtime0 outputs retain numeric binding.
    for name in ('U','p','phi'):
        copied(flow/'2000'/name,case/'0'/name)
        copied(flow/'2000'/name,case/'provenance/baseline-native-inputs'/name)
    copied_provenance(case,[flow/'provenance/setup.json',flow/'provenance/geometry.npz'])
    copied(audit/'checks.json',case/'provenance/flow-checks.json');copied(audit/'source-manifest.json',case/'provenance/flow-source-manifest.json')
    src=flow/'provenance/installed-sources'
    for p in sorted(src.rglob('*')):
        if p.is_file():copied(p,case/'provenance/installed-sources'/p.relative_to(src))
    copied(flow/'provenance/installed-source-sha256.txt',case/'provenance/installed-source-sha256.txt')
    references={}
    for name in ('U','p','phi','C','Vc'):
        require(sha(flow/'2000'/name)==sha(case/'2000'/name),'Exact frozen native copy required')
        references[name]=dict(flow_path=str(flow/'2000'/name),flow_sha256=sha(flow/'2000'/name),pilot_path='2000/'+name,pilot_sha256=sha(case/'2000'/name))
    json_new(case/'provenance/frozen-flow-references.json',dict(actual_flow_case=str(flow),flow_snapshot_steady_iteration=2000,flow_labels_are_not_thermal_time=True,fields=references))
    declaration=dict(phase='thermal_pilot_only',solver='OpenFOAM Foundation14 functions/subSolver incompressibleFluid',frozen_read_label=2000,thermal_start_time_s=0,thermal_end_time_s=.005,
      deltaT_s=.005,actual_steps_authorized=1,prospective_production_steps=600,prospective_production_end_time_s=3,production_admitted=False,
      rho_kg_m3=1,cp_J_kg_K=1000,k_W_m_K=10,alpha_m2_s=.01,depth_m=.01,initial_T_K=300,inlet_T_K=300,wall_T_K=350,
      thermal_advection='Gauss upwind',thermal_ddt='Euler',thermal_laplacian='Gauss linear orthogonal',thermal_nCorrectors=0,thermal_equation_relaxation=1,
      final_linear_T_residual_target=1e-8,T_bounds_tolerance_K=.01,frozen_field_max_abs_change=1e-12,global_Euler_relative_energy_target=1e-5,reconstructed_PDE_boundary_flux_target=.005,
      energy_reference_T_K=300,energy_normalization='per native depth; max(sum(abs(Qin_zone)),abs(storage),5000W/m)',native_diffusion='k*actual_boundary_gradT_dot_outward_normal*actual_face_area',
      upwind_advection='-rho*cp*actual_phi*(actual_upwind_T-300)',native_gradient_source='Installed gaussGrad correctBoundaryConditions explicitly imposes n.grad(T)=Tpatch.snGrad()',
      native_frozen_snapshot_source='writeObjects anyWrite exports registered U,p,phi without solving or regenerating phi',target_document_sha256=TARGET_SHA,physical_validation=False)
    json_new(case/'declared-comparison.json',declaration);manifest_new(case);json_new(case/'run-status.json',dict(stage='thermal_pilot_prepared',solver_started=False,production_admitted=False))
    return dict(case=str(case),solver_started=False,thermal_steps_authorized=1,production_admitted=False)

def verify_frozen(case):
    case=checked_case(case,'pilot');r=json.loads((case/'provenance/frozen-flow-references.json').read_text())
    require(all(sha(Path(v['flow_path']))==v['flow_sha256'] and sha(case/v['pilot_path'])==v['pilot_sha256'] for v in r['fields'].values()),'Actual frozen flow sources/copies changed')
    require(all(sha(case/'provenance/baseline-native-inputs'/n)==r['fields'][n]['flow_sha256'] for n in ('U','p','phi')),'Immutable baseline native flow input copies changed')
    for n,count,components,dims in [('U',CELLS,3,[0,1,-1,0,0,0,0]),('p',CELLS,1,[0,2,-2,0,0,0,0]),('phi',INTERNAL_FACES,1,[0,3,-1,0,0,0,0])]:
        baseline=native_field_values(case/'0'/n,count,components,dims);reference=native_field_values(case/'2000'/n,count,components,dims)
        left=[v for row in baseline for v in row] if components>1 else baseline;right=[v for row in reference for v in row] if components>1 else reference
        require(max(abs(a-b) for a,b in zip(left,right))<=1e-12,'Runtime baseline internal '+n+' differs from actual frozen reference')
    return dict(frozen_source_bytes_unchanged=True,actual_native_phi_preserved=True,immutable_baseline_input_bytes_match_flow=True,runtime_baseline_internal_U_p_phi_match_frozen=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--case',required=True,type=Path)
    mode=parser.add_mutually_exclusive_group(required=True);mode.add_argument('--mesh',type=Path);mode.add_argument('--validate-native',action='store_true');mode.add_argument('--prepare-thermal-pilot',action='store_true');mode.add_argument('--verify-frozen',action='store_true')
    a=parser.parse_args()
    result=prepare(a.mesh,a.case) if a.mesh else validate_native(a.case) if a.validate_native else prepare_pilot(a.case) if a.prepare_thermal_pilot else verify_frozen(a.case)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
