"""Prepare/check the dedicated laminar planar 90-degree channel bend; starts no solver."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import shutil

RUN = '2026-11-01-elbow-v001'
PATCHES = ('inlet', 'outlet', 'walls')
CELLS=3200
INTERNAL_FACES=6220
TOTAL_FACES=12980
POINTS=6762
EMPTY = 'frontAndBackPlanes'
FLOAT = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
SOURCES = {
 'velocity': 'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/tutorials/incompressibleFluid/pitzDailySteady/0/U',
 'pressure': 'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/tutorials/incompressibleFluid/pitzDailySteady/0/p',
 'steady_schemes': 'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/tutorials/incompressibleFluid/pitzDailySteady/system/fvSchemes',
 'steady_SIMPLE': 'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/tutorials/incompressibleFluid/pitzDailySteady/system/fvSolution',
 'noSlip': 'https://cpp.openfoam.org/v14/classFoam_1_1noSlipFvPatchVectorField.html',
}

def require(value, message):
    if not value:
        raise ValueError(message)
def lab_root():
    return Path(__file__).resolve().parents[3]
def checked_case(path):
    path = Path(path).resolve()
    require(path == (lab_root()/'native/runs'/RUN/'openfoam').resolve(), 'Only dedicated elbow-v001/openfoam allowed')
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
def dictionaries():
    content = {
        'system/controlDict': '''solver incompressibleFluid;
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
''',
        'system/fvSchemes': '''ddtSchemes { default steadyState; }
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
        'system/fvSolution': '''solvers
{
 p { solver GAMG; tolerance 1e-10; relTol 0.01; smoother GaussSeidel; }
 U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01; }
}
SIMPLE
{
 nNonOrthogonalCorrectors 0;
 consistent yes;
}
relaxationFactors { equations { U 0.7; } }
''',
        'system/decomposeParDict': '''numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
''',
        'constant/physicalProperties': 'viscosityModel constant;\nnu 0.001;\n',
        'constant/momentumTransport': 'simulationType laminar;\n',
        '0/U': '''dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
 inlet { type fixedValue; value uniform (0.1 0 0); }
 outlet { type zeroGradient; }
 walls { type noSlip; }
 frontAndBackPlanes { type empty; }
}
''',
        '0/p': '''dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
 inlet { type zeroGradient; }
 outlet { type fixedValue; value uniform 0; }
 walls { type zeroGradient; }
 frontAndBackPlanes { type empty; }
}
''',
    }
    return {name: header(Path(name).name, 'volVectorField' if name=='0/U' else 'volScalarField' if name=='0/p' else 'dictionary')+body for name,body in content.items()}
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
    count,body=counted(poly/'boundary');rows={}
    for name,block in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',body):
        rows[name.strip('"')]=dict(kind=re.search(r'\btype\s+(\w+)\s*;',block)[1],
            count=int(re.search(r'\bnFaces\s+(\d+)\s*;',block)[1]),
            start=int(re.search(r'\bstartFace\s+(\d+)\s*;',block)[1]))
    require(count==4 and set(rows)==set(PATCHES)|{EMPTY},'Unexpected native elbow patches')
    require(rows['walls']['kind']=='wall' and rows['walls']['count']==320,'Require320 native wall faces')
    require(all(rows[n]['count']==20 and rows[n]['kind']=='patch' for n in ('inlet','outlet')),'Each opening must be20 native patch faces')
    require(rows[EMPTY]['kind']=='empty' and rows[EMPTY]['count']==6400,'Require6400 empty plane faces')
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
    require(all(abs(a-b)<=1e-10 for row,target in zip(bounds,([-.6,.2],[0,1],[-.005,.005])) for a,b in zip(row,target)), 'Native planar channel bend/extrusion bounds differ')
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
    require(not case.exists(),'Fresh native elbow case required; refuse overwrite')
    require(mesh==case.parent/'mesh/elbow.msh','Require selected shared elbow.msh')
    raw=mesh.read_bytes()
    require(re.search(rb'\(2\s+2\)',raw),'Require2D Fluent ASCII')
    m=re.search(rb'^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)',raw,re.M)
    require(m is not None and int(m[1],16)==CELLS,'Require3200 shared quads')
    for name in ('fluid','interior',*PATCHES):
        require(re.search(rb'\(45\s+\([^)]*\s'+name.encode()+rb'\)\(\)\)',raw),'Missing shared zone '+name)
    setup=json.loads((mesh.parent/'setup.json').read_text(encoding='utf-8-sig'))
    require(setup['cells']==CELLS and setup['mesh_sha256']==sha(mesh),'Shared hash/count differs')
    for name,value in dict(dx_m=.01,dy_m=.01,density_kg_m3=1,dynamic_viscosity_Pa_s=.001,kinematic_viscosity_m2_s=.001,Re_inlet=20,reference_length_m=.2,reference_velocity_m_s=.1,openfoam_extrusion_depth_m=.01,operating_pressure_Pa=0,outlet_gauge_pressure_Pa=0,initial_pressure_Pa=0).items():
        require(setup[name]==value,'Shared physical setup differs '+name)
    require(setup['inlet_velocity_m_s']==[.1,0] and setup['initial_velocity_m_s']==[0,0],'Shared inlet/initial vectors differ')
    require(setup['geometry_npz_sha256']==sha(mesh.parent/'geometry.npz'),'Shared geometry hash differs')
    require(shutil.disk_usage(case.anchor).free>=10*1024**3,'Require10GiB output disk')
    case.mkdir(parents=True)
    for name,text in dictionaries().items():
        write_new(case/name,text)
    copied(mesh,case/'input-elbow.msh')
    for source in (Path(__file__),lab_root()/'local/elbow-v001/run-openfoam.sh',lab_root()/'local/elbow-v001/Launch-OpenFOAM.ps1',
            mesh.parent/'setup.json',mesh.parent/'geometry.npz',lab_root()/'local/elbow-v001/DESIGN.md'):
        copied(source,case/'provenance'/source.name)
    declaration=dict(model='2D steady laminar planar 90-degree channel bend',intended_cells=CELLS,dx_m=.01,dy_m=.01,
        density_kg_m3=1,dynamic_viscosity_Pa_s=.001,kinematic_viscosity_m2_s=.001,Re_inlet=20,
        reference_length_m=.2,reference_velocity_m_s=.1,extrusion_depth_m=.01,x_bounds_m=[-.6,.2],y_bounds_m=[0,1],z_bounds_m=[-.005,.005],
        geometry='Union horizontal x[-.6,.2] y[0,.2] and upper vertical x[0,.2] y[.2,1]; sparse3200-cell occupancy in100x80 background',
        inlet_U_m_s=[.1,0,0],outlet='zeroGradient U,fixedValue kinematic p0',walls='stationary noSlip',
        initial_U_m_s=[0,0,0],initial_p=0,pressure='kinematic p; gauge Pa=rho*p; outlet0; zeroGradient p at inlet and walls',
        solver='OpenFOAM Foundation14 incompressibleFluid',algorithm='steady SIMPLE consistent',
        U_equation_relaxation=.7,nonorthogonal_correctors=0,gradient='Gauss linear',advection='bounded Gauss linearUpwind grad(U)',
        mpi_processes=2,max_iterations=2000,saved_iterations=list(range(100,2001,100)),selected_iterations=[500,1000,1500,1800,2000],stability_iterations=[1800,2000],
        printed_p_Ux_Uy_residual_target=1e-6,relative_L2_U_target=1e-3,
        native_local_divergence_definition='sum signed outward native face phi / actual oriented native cell volume * H/Uref; H/Uref=2s',
        native_local_divergence_max_target=1e-5,native_local_divergence_RMS_target=1e-6,relative_all_boundary_mass_imbalance_target=1e-4,
        mass_normalization='abs(sum signed outward phi on inlet,outlet,walls) / actual sum magnitude inward native phi at inlet',
        expected_per_depth_volume_flow_m2_s=dict(inlet=.02,outlet=.02),
        mesh_source_path=str(mesh),mesh_source_sha256=sha(mesh),geometry_sha256=sha(mesh.parent/'geometry.npz'),official_sources=SOURCES,
        physical_validation=False,mesh_independence_established=False,residual_norms_equivalent_to_fluent=False,scalar_mixing_simulated=False,
        three_dimensional_secondary_flow_established=False)
    json_new(case/'declared-comparison.json',declaration)
    json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file()})
    json_new(case/'run-status.json',dict(stage='prepared_only',solver_started=False))
    return dict(case=str(case),cells=CELLS,solver_started=False,mesh_sha256=sha(mesh))
def validate_native(case):
    case=checked_case(case)
    require(not (case/'solver.log').exists() and not list(case.glob('processor*')),'Native preflight precedes decomposition/solve')
    manifest=json.loads((case/'provenance/prepared-input-sha256.json').read_text())
    require(all(sha(case/name)==digest for name,digest in manifest.items()),'Prepared inputs changed')
    poly=native_poly(case); patches=native_patches(poly)
    quality=mesh_quality((case/'mesh-check.log').read_text(encoding='utf-8'))
    geometry=native_geometry(poly)
    receipt=dict(passed=True,solver_started=False,actual_native_patches=patches,actual_mesh_quality=quality,native_bounds_m=geometry['bounds'],
        points=len(geometry['points']),faces=len(geometry['faces']),internal_faces=len(geometry['neighbour']),cells=len(geometry['volumes']),
        actual_oriented_volume_min_m3=min(geometry['volumes']),actual_oriented_volume_max_m3=max(geometry['volumes']),
        actual_oriented_total_volume_m3=sum(geometry['volumes']),extrusion_depth_m=.01,prepared_inputs_unchanged=True,
        mesh_location=poly.relative_to(case).as_posix(),native_mesh_sha256={name:sha(poly/name) for name in ('points','faces','owner','neighbour','boundary')})
    json_new(case/'provenance/preflight-validation.json',receipt)
    return receipt
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--case',required=True,type=Path)
    group=parser.add_mutually_exclusive_group(required=True);group.add_argument('--mesh',type=Path);group.add_argument('--validate-native',action='store_true')
    args=parser.parse_args()
    print(json.dumps(prepare(args.mesh,args.case) if args.mesh else validate_native(args.case),indent=2))
if __name__=='__main__':
    main()

