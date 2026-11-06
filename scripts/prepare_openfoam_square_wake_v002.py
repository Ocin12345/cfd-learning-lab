"""Prepare/check the dedicated steady laminar square-obstacle wake; starts no solver."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import shutil

RUN = '2026-11-01-square-wake-v001'
PATCHES = ('inlet','outlet','top','bottom','body')
CELLS=15900
INTERNAL_FACES=31500
TOTAL_FACES=63900
POINTS=32400
EMPTY = 'frontAndBackPlanes'
FLOAT = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
INSTALLED_SOURCES=(
 'applications/modules/incompressibleFluid/incompressibleFluid.C',
 'src/functionObjects/forces/forces/forces.C',
 'src/functionObjects/forces/forcesBase/forcesBase.C',
 'src/functionObjects/forces/forceCoeffs/forceCoeffs.C',
 'src/functionObjects/field/wallShearStress/wallShearStress.C',
 'src/OpenFOAM/db/functionObjects/functionObjectList/functionObjectList.C',
 'applications/utilities/postProcessing/foamPostProcess/foamPostProcess.C',
 'src/functionObjects/field/writeCellCentres/writeCellCentres.C',
 'src/functionObjects/field/writeCellVolumes/writeCellVolumes.C',
 'src/finiteVolume/fields/fvPatchFields/constraint/symmetryPlane/symmetryPlaneFvPatchField.C',
 'tutorials/incompressibleFluid/pitzDailySteady/system/fvSchemes',
 'tutorials/incompressibleFluid/pitzDailySteady/system/fvSolution',
)
SOURCES={n:'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/'+n for n in INSTALLED_SOURCES}

def require(value, message):
    if not value:
        raise ValueError(message)
def lab_root():
    return Path(__file__).resolve().parents[3]
def checked_case(path):
    path = Path(path).resolve()
    require(path == (lab_root()/'native/runs'/RUN/'openfoam-v002').resolve(), 'Only dedicated square-wake-v001/openfoam-v002 allowed')
    return path
def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()
def write_new(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='ascii', newline='\n') as stream:
        stream.write(text.rstrip()+'\n')
def json_new(path, value):
    write_new(path, json.dumps(value, indent=2, allow_nan=False))
def copied(source, target):
    require(not target.exists(), 'Refuse overwrite '+str(target))
    target.parent.mkdir(parents=True, exist_ok=True)
    before = sha(source)
    shutil.copy2(source, target)
    require(before == sha(source) == sha(target), 'Copy hash differs')
def header(name, kind='dictionary'):
    return 'FoamFile\n{\n version 2.0;\n format ascii;\n class '+kind+';\n object '+name+';\n}\n\n'
def force_function(name,kind,interval):
    text=f'''{name}
{{
 type {kind};
 libs ("libforces.so");
 patches (body);
 rho rhoInf;
 rhoInf 1;
 p p;
 U U;
 pRef 0;
 porosity false;
 CofR (0.5 0 0);
 log false;
 executeControl timeStep;
 executeInterval {interval};
 writeControl timeStep;
 writeInterval {interval};
'''
    if kind=='forceCoeffs':
        text+=' liftDir (0 1 0);\n dragDir (1 0 0);\n pitchAxis (0 0 1);\n magUInf 0.1;\n lRef 0.1;\n Aref 0.001;\n'
    return text+'}\n'
def shear_function(name,interval):
    return f'''{name}
{{
 type wallShearStress;
 libs ("libfieldFunctionObjects.so");
 patches (body);
 log false;
 executeControl timeStep;
 executeInterval {interval};
 writeControl timeStep;
 writeInterval {interval};
}}
'''
def dictionaries():
    main_functions=force_function('bodyForces','forces',1)+force_function('bodyCoefficients','forceCoeffs',1)+shear_function('bodyShear',100)
    probe_functions=force_function('probeForces','forces',1)+force_function('probeCoefficients','forceCoeffs',1)+shear_function('probeShear',1)
    content={
      'system/controlDict':'''solver incompressibleFluid;
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
timePrecision 8;
runTimeModifiable false;
functions
{
'''+main_functions+'}\n',
      'system/forcePreflight':probe_functions,
      'system/fvSchemes':'''ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
 default none;
 div(phi,U) bounded Gauss linearUpwind grad(U);
 div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''',
      'system/fvSolution':'''solvers
{
 p { solver GAMG; tolerance 1e-10; relTol 0.01; smoother GaussSeidel; }
 U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01; }
}
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; }
relaxationFactors { equations { U 0.7; } }
''',
      'system/decomposeParDict':'''numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
''',
      'constant/physicalProperties':'viscosityModel constant;\nnu 0.0005;\n',
      'constant/momentumTransport':'simulationType laminar;\n',
      '0/U':'''dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
 inlet { type fixedValue; value uniform (0.1 0 0); }
 outlet { type zeroGradient; }
 top { type symmetryPlane; }
 bottom { type symmetryPlane; }
 body { type noSlip; }
 frontAndBackPlanes { type empty; }
}
''',
      '0/p':'''dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
 inlet { type zeroGradient; }
 outlet { type fixedValue; value uniform 0; }
 top { type symmetryPlane; }
 bottom { type symmetryPlane; }
 body { type zeroGradient; }
 frontAndBackPlanes { type empty; }
}
''',
    }
    return {name:header(Path(name).name,'volVectorField' if name=='0/U' else 'volScalarField' if name=='0/p' else 'dictionary')+body for name,body in content.items()}
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
    require(n==6 and set(rows)==set(PATCHES)|{EMPTY},'Native square-wake patch set differs')
    require(all(rows[n]['kind']=='patch' and rows[n]['count']==80 for n in ('inlet','outlet')),'Native opening counts/types differ')
    require(all(rows[n]['kind']=='symmetryPlane' and rows[n]['count']==200 for n in ('top','bottom')),'Native slip/symmetry plane counts/types differ')
    require(rows['body']['kind']=='wall' and rows['body']['count']==40,'Native body wall must be40faces')
    require(rows[EMPTY]['kind']=='empty' and rows[EMPTY]['count']==31800,'Native empty face count differs')
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
    require(all(abs(a-b)<=1e-10 for row,target in zip(bounds,([0,2],[-.4,.4],[-.005,.005])) for a,b in zip(row,target)), 'Native square-obstacle wake/extrusion bounds differ')
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
    expected=.01*.01*.01
    require(all(math.isfinite(v) and v>0 and abs(v-expected)<=expected*1e-8 for v in volumes),'Actual oriented native cell volumes differ from uniform extrusion')
    return dict(points=points,faces=faces,owner=owner,neighbour=neighbour,volumes=volumes,bounds=bounds)
def prepare(mesh,case):
    mesh,case=Path(mesh).resolve(),checked_case(case)
    require(not case.exists(),'Fresh native square-wake case required')
    require(mesh==case.parent/'mesh/square-wake.msh','Require selected shared square-wake mesh')
    raw=mesh.read_bytes();require(re.search(rb'\(2\s+2\)',raw),'Require2D FluentASCII')
    m=re.search(rb'^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)',raw,re.M)
    require(m is not None and int(m[1],16)==CELLS,'Require15900sharedquads')
    for n in ('fluid','interior',*PATCHES):require(re.search(rb'\(45\s+\([^)]*\s'+n.encode()+rb'\)\(\)\)',raw),'Missing named source zone '+n)
    setup=json.loads((mesh.parent/'setup.json').read_text(encoding='utf-8-sig'))
    require(setup['cells']==CELLS and setup['mesh_sha256']==sha(mesh),'Shared mesh count/hash differs')
    for n,v in dict(nx=200,ny=80,dx_m=.01,dy_m=.01,density_kg_m3=1,dynamic_viscosity_Pa_s=.0005,kinematic_viscosity_m2_s=.0005,
      Re_body=20,body_side_m=.1,reference_length_m=.1,reference_velocity_m_s=.1,reference_area_m2=.1,reference_depth_m=1,openfoam_extrusion_depth_m=.01,operating_pressure_Pa=0,outlet_gauge_pressure_Pa=0).items():
        require(setup[n]==v,'Shared physical setup differs '+n)
    require(setup['initial_velocity_m_s']==[0,0] and setup['initial_pressure_Pa']==0 and setup['inlet_velocity_m_s']==[.1,0] and setup['geometry_npz_sha256']==sha(mesh.parent/'geometry.npz'),'Shared inlet/geometry hash differs')
    require(shutil.disk_usage(case.anchor).free>=10*1024**3,'Require10GiBdisk')
    case.mkdir(parents=True)
    for n,text in dictionaries().items():write_new(case/n,text)
    copied(mesh,case/'input-square-wake.msh')
    for src in (Path(__file__),lab_root()/'local/square-wake-v001/run-openfoam-v002.sh',lab_root()/'local/square-wake-v001/Launch-OpenFOAM-v002.ps1',mesh.parent/'setup.json',mesh.parent/'geometry.npz',lab_root()/'local/square-wake-v001/DESIGN.md',lab_root()/'local/square-wake-v001/FORCE_TARGETS.md',lab_root()/'local/square-wake-v001/mesh-independent-review-v001.json'):
        copied(src,case/'provenance'/src.name)
    declaration=dict(model='2D steady laminar confined square-obstacle wake',intended_cells=CELLS,nx=200,ny=80,dx_m=.01,dy_m=.01,body_bounds_m=[[.45,.55],[-.05,.05]],
      domain_bounds_m=[[0,2],[-.4,.4]],extrusion_depth_m=.01,density_kg_m3=1,dynamic_viscosity_Pa_s=.0005,kinematic_viscosity_m2_s=.0005,Re_body=20,
      inlet_U_m_s=[.1,0,0],initial_U_m_s=[0,0,0],initial_p=0,outlet='zeroGradient U/fixedValue p0',body='noSlip U/zeroGradient p',top_bottom='symmetryPlane U/p; zero normal velocity',
      solver='OpenFOAM Foundation14 incompressibleFluid',algorithm='steady SIMPLE consistent',U_relaxation=.7,nonorthogonal_correctors=0,advection='bounded Gauss linearUpwind grad(U)',
      max_iterations=2000,saved_iterations=list(range(100,2001,100)),selected_iterations=[500,1000,1500,1800,2000],mpi_processes=2,
      printed_p_Ux_Uy_residual_target=1e-6,local_divergence_scale_D_over_U_s=1,local_divergence_max_target=1e-5,local_divergence_RMS_target=1e-6,
      relative_all_boundary_mass_imbalance_target=1e-4,stability_iterations=[1800,2000],U_relative_L2_target=1e-3,
      force_reference_rho_kg_m3=1,force_reference_U_m_s=.1,force_reference_D_m=.1,force_reference_area_m2=.001,force_denominator_N=5e-6,
      force_moment_reference_m=[.5,0,0],force_patches=['body'],force_pressure_reference_Pa=0,force_sample_interval=1,
      force_coefficient_absolute_crosscheck_target=1e-8,body_field_force_coefficient_equivalent_crosscheck_target=1e-8,
      last200_force_window=[1801,2000],last200_Cd_peak_to_peak_over_abs_mean_target=1e-3,last200_abs_mean_Cl_target=1e-3,
      zero_mean_Cd_policy='Undefined relative range yields false numerical gate with null metric/reason, never parser failure',
      force_field_reconstruction='PressureF=sum(rho*p_boundary*Sf); zeroGradient p uses actual ownerp; viscousF=-rho*sum(nativewallShearStress*facearea)',
      pressure_comparison='Raw gauge Pa=rho*p, common outlet0',mesh_source_sha256=sha(mesh),geometry_sha256=sha(mesh.parent/'geometry.npz'),
      force_target_document_sha256=sha(lab_root()/'local/square-wake-v001/FORCE_TARGETS.md'),official_sources=SOURCES,physical_validation=False,mesh_independence_established=False,unbounded_freestream_established=False,blockage_ratio=.125,domain_height_over_D=8,
      iteration_indices_are_not_physical_time=True,transient_shedding_established=False)
    json_new(case/'declared-comparison.json',declaration)
    json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file()})
    json_new(case/'run-status.json',dict(stage='prepared_only',solver_started=False))
    return dict(case=str(case),cells=CELLS,mesh_sha256=sha(mesh),solver_started=False)
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
def validate_native(case):
    case=checked_case(case);require(not (case/'solver.log').exists() and not list(case.glob('processor*')),'Preflight must precede solver/decomposition')
    manifest=json.loads((case/'provenance/prepared-input-sha256.json').read_text());require(all(sha(case/n)==v for n,v in manifest.items()),'Prepared source changed')
    poly=native_poly(case);patches=native_patches(poly);geometry=native_geometry(poly);quality=mesh_quality((case/'mesh-check.log').read_text())
    C=native_field_values(case/'0/C',CELLS,3,[0,1,0,0,0,0,0]);V=native_field_values(case/'0/Vc',CELLS,1,[0,3,0,0,0,0,0])
    require(all(v>0 and abs(v-w)<=1e-14 for v,w in zip(V,geometry['volumes'])),'Native Vc differs from oriented volume')
    keys=[]
    for x,y,z in C:
        ix=round(x/.01-.5);iy=round((y+.4)/.01-.5)
        require(0<=ix<200 and 0<=iy<80 and not(45<=ix<55 and 35<=iy<45) and abs(x-(ix+.5)*.01)<=1e-10 and abs(y-(-.4+(iy+.5)*.01))<=1e-10 and abs(z)<=1e-10,'Native C off source hole grid')
        keys.append(iy*200+ix)
    expected=[j*200+i for j in range(80) for i in range(200) if not(45<=i<55 and 35<=j<45)]
    require(sorted(keys)==expected,'Actual native C hole occupancy differs')
    def numeric_rows(path,columns):
        text=path.read_text();require(not re.search(r'FOAM FATAL|(?<![A-Za-z0-9_])[-+]?(?:nan|inf(?:inity)?)(?![A-Za-z0-9_])',text,re.I),'Invalid force probe')
        rows=[list(map(float,re.findall(FLOAT,l))) for l in text.splitlines() if l.strip() and not l.lstrip().startswith('#')]
        require(len(rows)==1 and len(rows[0])==columns and rows[0][0]==0 and all(map(math.isfinite,rows[0])),'Require one actual iteration0 force probe row')
        return rows[0]
    forces=numeric_rows(case/'postProcessing/probeForces/0/forces.dat',13)
    coeff=numeric_rows(case/'postProcessing/probeCoefficients/0/forceCoeffs.dat',6)
    Cd=(forces[1]+forces[4])/5e-6;Cl=(forces[2]+forces[5])/5e-6
    require(abs(Cd-coeff[2])<=1e-8 and abs(Cl-coeff[3])<=1e-8,'Actual native registered force/coefficient probe differs')
    require((case/'0/wallShearStress').is_file(),'Actual wallShearStress force reconstruction probe missing')
    receipt=dict(passed=True,solver_started=False,actual_native_patches=patches,actual_mesh_quality=quality,native_bounds_m=geometry['bounds'],points=POINTS,faces=TOTAL_FACES,
      internal_faces=INTERNAL_FACES,cells=CELLS,actual_oriented_volume_min_m3=min(geometry['volumes']),actual_oriented_volume_max_m3=max(geometry['volumes']),
      actual_oriented_total_volume_m3=sum(geometry['volumes']),actual_native_C_hole_grid_bijective=True,force_and_coefficient_registration_at_iteration0_passed=True,
      force_probe_Cd=Cd,force_probe_Cl=Cl,prepared_inputs_unchanged=True,mesh_location=poly.relative_to(case).as_posix(),
      native_mesh_sha256={n:sha(poly/n) for n in ('points','faces','owner','neighbour','boundary')},native_initial_C_Vc_sha256={n:sha(case/'0'/n) for n in ('C','Vc')})
    json_new(case/'provenance/preflight-validation.json',receipt);return receipt
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--case',required=True,type=Path)
    group=parser.add_mutually_exclusive_group(required=True);group.add_argument('--mesh',type=Path);group.add_argument('--validate-native',action='store_true')
    args=parser.parse_args()
    print(json.dumps(prepare(args.mesh,args.case) if args.mesh else validate_native(args.case),indent=2))
if __name__=='__main__':
    main()

