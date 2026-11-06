"""Prepare/check the dedicated transient solid heat-conduction slab; starts no solver."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import shutil

RUN = '2026-11-01-thermal-slab-v001'
PATCHES = ('hot','cold','insulated_bottom','insulated_top')
CELLS=2000
INTERNAL_FACES=3880
TOTAL_FACES=8120
POINTS=4242
EMPTY = 'frontAndBackPlanes'
FLOAT = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
REFERENCE_SHA='550c050d4294c08d0fb3905d06a2eeda7c82bb44397072a3dd261db1a6d9846a'
INSTALLED_SOURCES=(
 'applications/solvers/foamRun/foamRun.C',
 'applications/modules/solid/solid.C',
 'src/thermophysicalModels/solidThermo/constSolidThermo/constSolidThermo.C',
 'src/thermophysicalModels/solidThermo/solidThermo/solidThermo.C',
 'src/ThermophysicalTransportModels/solid/isotropic/isotropic.C',
 'src/ThermophysicalTransportModels/solid/solidThermophysicalTransportModel/solidThermophysicalTransportModel.C',
 'src/functionObjects/field/wallHeatFlux/wallHeatFlux.C',
 'src/functionObjects/field/wallHeatFlux/wallHeatFlux.H',
 'src/functionObjects/field/writeCellCentres/writeCellCentres.C',
 'src/functionObjects/field/writeCellVolumes/writeCellVolumes.C',
 'src/finiteVolume/finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.C',
 'tutorials/multiRegion/CHT/coolingCylinder2D/constant/solid/physicalProperties',
 'tutorials/multiRegion/CHT/coolingCylinder2D/system/solid/fvSolution',
 'tutorials/multiRegion/CHT/coolingCylinder2D/system/solid/fvSchemes',
)
SOURCES={name:'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/'+name for name in INSTALLED_SOURCES}

def require(value, message):
    if not value:
        raise ValueError(message)
def lab_root():
    return Path(__file__).resolve().parents[3]
def checked_case(path):
    path = Path(path).resolve()
    require(path == (lab_root()/'native/runs'/RUN/'openfoam-v002').resolve(), 'Only dedicated thermal-slab-v001/openfoam-v002 allowed')
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
    content={
      'system/controlDict':'''solver solid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 200;
deltaT 2;
adjustTimeStep false;
writeControl timeStep;
writeInterval 1;
purgeWrite 0;
writeFormat ascii;
writePrecision 14;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable false;
functions
{
 wallHeatRates
 {
  type wallHeatFlux;
  libs ("libfieldFunctionObjects.so");
  patches (hot cold insulated_bottom insulated_top);
  log false;
  executeControl timeStep;
  executeInterval 1;
  writeControl timeStep;
  writeInterval 1;
 }
}
''',
      'system/fvSchemes':'''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear uncorrected; }
interpolationSchemes { default linear; }
snGradSchemes { default uncorrected; }
''',
      'system/fvSolution':'''solvers
{
 e { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }
 eFinal { $e; relTol 0; }
}
PIMPLE { nOuterCorrectors 1; nNonOrthogonalCorrectors 0; }
relaxationFactors { equations { e 1; } }
''',
      'system/decomposeParDict':'''numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
''',
      'constant/physicalProperties':'''thermoType constSolidThermo;
rho { type uniform; value 1000; }
Cv { type uniform; value 1000; }
kappa { type uniform; value 10; }
''',
      '0/T':'''dimensions [0 0 0 1 0 0 0];
internalField uniform 300;
boundaryField
{
 hot { type fixedValue; value uniform 400; }
 cold { type fixedValue; value uniform 300; }
 insulated_bottom { type zeroGradient; }
 insulated_top { type zeroGradient; }
 frontAndBackPlanes { type empty; }
}
''',
    }
    return {name:header(Path(name).name,'volScalarField' if name=='0/T' else 'dictionary')+body for name,body in content.items()}
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
            count=int(re.search(r'\bnFaces\s+(\d+)\s*;',block)[1]),start=int(re.search(r'\bstartFace\s+(\d+)\s*;',block)[1]))
    require(count==5 and set(rows)==set(PATCHES)|{EMPTY},'Unexpected native slab patch names')
    require(all(rows[n]['kind']=='wall' and rows[n]['count']==(20 if n in ('hot','cold') else 100) for n in PATCHES),'All physical thermal patches must be wall type with exact counts')
    require(rows[EMPTY]['kind']=='empty' and rows[EMPTY]['count']==4000,'Require4000 empty plane faces')
    ranges=sorted((v['start'],v['start']+v['count']) for v in rows.values())
    require(ranges[0][0]==INTERNAL_FACES and ranges[-1][1]==TOTAL_FACES and all(a[1]==b[0] for a,b in zip(ranges,ranges[1:])),'Native boundary face ranges differ')
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
    require(all(abs(a-b)<=1e-10 for row,target in zip(bounds,([0,.1],[0,.02],[-.005,.005])) for a,b in zip(row,target)), 'Native slab/extrusion bounds differ')
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
    expected=.001*.001*.01
    require(all(math.isfinite(v) and v>0 and abs(v-expected)<=expected*1e-8 for v in volumes),'Actual oriented native cell volumes differ from uniform extrusion')
    return dict(points=points,faces=faces,owner=owner,neighbour=neighbour,volumes=volumes,bounds=bounds)
def prepare(mesh,case):
    mesh,case=Path(mesh).resolve(),checked_case(case)
    require(not case.exists(),'Fresh native thermal slab case required; refuse overwrite')
    require(mesh==case.parent/'mesh/slab.msh','Require selected shared slab.msh')
    raw=mesh.read_bytes();require(re.search(rb'\(2\s+2\)',raw),'Require2D Fluent ASCII')
    m=re.search(rb'^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)',raw,re.M)
    require(m is not None and int(m[1],16)==CELLS,'Require2000 shared quads')
    for name in ('slab','interior',*PATCHES):
        require(re.search(rb'\(45\s+\([^)]*\s'+name.encode()+rb'\)\(\)\)',raw),'Missing shared zone '+name)
    setup=json.loads((mesh.parent/'setup.json').read_text(encoding='utf-8-sig'))
    require(setup['cells']==CELLS and setup['mesh_sha256']==sha(mesh),'Shared mesh hash/count differs')
    properties=dict(nx=100,ny=20,length_m=.1,height_m=.02,dx_m=.001,dy_m=.001,density_kg_m3=1000,specific_heat_J_kg_K=1000,
      thermal_conductivity_W_m_K=10,thermal_diffusivity_m2_s=1e-5,initial_temperature_K=300,hot_temperature_K=400,cold_temperature_K=300,
      time_step_s=2,time_steps=100,final_time_s=200,openfoam_extrusion_depth_m=.01)
    for name,value in properties.items():require(setup[name]==value,'Shared thermal setup differs '+name)
    require(setup['cell_zone_name']=='slab','Shared slab zone name differs')
    require(setup['geometry_npz_sha256']==sha(mesh.parent/'geometry.npz'),'Shared geometry hash differs')
    reference=lab_root()/'repo/cfd-learning-lab/scripts/thermal_slab_reference_v001.py'
    require(sha(reference)==REFERENCE_SHA,'Frozen analytical reference changed')
    require(shutil.disk_usage(case.anchor).free>=10*1024**3,'Require10GiB output disk')
    case.mkdir(parents=True)
    for name,text in dictionaries().items():write_new(case/name,text)
    copied(mesh,case/'input-slab.msh')
    for source in (Path(__file__),reference,lab_root()/'local/thermal-slab-v001/run-openfoam-v002.sh',lab_root()/'local/thermal-slab-v001/Launch-OpenFOAM-v002.ps1',
        mesh.parent/'setup.json',mesh.parent/'geometry.npz',lab_root()/'local/thermal-slab-v001/DESIGN.md'):
        copied(source,case/'provenance'/source.name)
    declaration=dict(model='Genuine transient 2D constant-property solid heat conduction',intended_cells=CELLS,nx=100,ny=20,
      length_m=.1,height_m=.02,extrusion_depth_m=.01,rho_kg_m3=1000,Cv_J_kg_K=1000,Cp_J_kg_K=1000,kappa_W_m_K=10,alpha_m2_s=1e-5,
      initial_T_K=300,hot_T_K=400,cold_T_K=300,insulated_patches=['insulated_bottom','insulated_top'],temperature_only=True,
      solver='OpenFOAM Foundation14 solid',material_model='constSolidThermo; Cp() returns Cv(); e=Cv*T internally, e NO_WRITE',
      time_scheme='Euler',time_step_s=2,time_steps=100,final_time_s=200,saved_times_s=list(range(0,201,2)),GIF_times_s=list(range(0,201,4)),
      mpi_processes=2,PIMPLE_outer_correctors=1,nonorthogonal_correctors=0,linear_energy_solver='PCG DIC; tolerance1e-12; relTol0',energy_relaxation=1,
      final_linear_energy_residual_target=1e-10,temperature_limits_K=[300,400],temperature_limit_tolerance_K=1e-6,y_spread_target_K=1e-6,
      native_heat_rate_convention='wallHeatFlux positive INTO solid; native field W/m2; native table patch QW integrates actualfacearea once',
      energy_reference_T_K=300,energy_definition='sum actualcells rho*Cv*V*(T-300)',Qref_W=2,
      Euler_energy_relative_definition='abs((En-Eprev)/dt-sumQin)/max(sum(abs(Qin_each_wall)),abs((En-Eprev)/dt),Qref)',
      Euler_energy_relative_target=1e-5,analytic_final_max_error_target_K=.5,analytic_final_RMS_error_target_K=.25,
      analytic_final_only_gate=True,analytic_series_terms=[80,160],analytic_series_difference_target_K=1e-9,
      reference_sha256=REFERENCE_SHA,mesh_source_path=str(mesh),mesh_source_sha256=sha(mesh),geometry_sha256=sha(mesh.parent/'geometry.npz'),
      official_sources=SOURCES,physical_validation=False,mesh_independence_established=False,time_step_independence_established=False,
      limitation='First-order dt2s coarsely resolves initial incompatible hot-boundary jump; all-time analytic errors retained without all-time accuracy gate')
    json_new(case/'declared-comparison.json',declaration)
    json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file()})
    json_new(case/'run-status.json',dict(stage='prepared_only',solver_started=False))
    return dict(case=str(case),cells=CELLS,solver_started=False,mesh_sha256=sha(mesh),reference_sha256=REFERENCE_SHA)
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
    case=checked_case(case)
    require(not (case/'solver.log').exists() and not list(case.glob('processor*')),'Native preflight must precede decomposition/solve')
    manifest=json.loads((case/'provenance/prepared-input-sha256.json').read_text())
    require(all(sha(case/name)==digest for name,digest in manifest.items()),'Prepared inputs changed')
    poly=native_poly(case);patches=native_patches(poly);quality=mesh_quality((case/'mesh-check.log').read_text(encoding='utf-8'))
    geometry=native_geometry(poly)
    centres=native_field_values(case/'0/C',CELLS,3,[0,1,0,0,0,0,0]);volumes=native_field_values(case/'0/Vc',CELLS,1,[0,3,0,0,0,0,0])
    require(all(v>0 and abs(v-w)<=1e-16 for v,w in zip(volumes,geometry['volumes'])),'Native emitted Vc differs from oriented volume')
    keys=[]
    for x,y,z in centres:
        ix=round(x/.001-.5);iy=round(y/.001-.5)
        require(0<=ix<100 and 0<=iy<20 and abs(x-(ix+.5)*.001)<=1e-10 and abs(y-(iy+.5)*.001)<=1e-10 and abs(z)<=1e-10,'Native emitted C off full slab grid')
        keys.append(iy*100+ix)
    require(sorted(keys)==list(range(CELLS)),'Native emitted C grid bijection failed')
    receipt=dict(passed=True,solver_started=False,actual_native_patches=patches,actual_mesh_quality=quality,native_bounds_m=geometry['bounds'],
      points=POINTS,faces=TOTAL_FACES,internal_faces=INTERNAL_FACES,cells=CELLS,actual_oriented_volume_min_m3=min(geometry['volumes']),
      actual_oriented_volume_max_m3=max(geometry['volumes']),actual_oriented_total_volume_m3=sum(geometry['volumes']),
      native_emitted_Vc_min_m3=min(volumes),native_emitted_Vc_max_m3=max(volumes),native_emitted_C_grid_bijective=True,solid_runtime_postprocess_constructed=True,
      extrusion_depth_m=.01,prepared_inputs_unchanged=True,mesh_location=poly.relative_to(case).as_posix(),
      native_mesh_sha256={name:sha(poly/name) for name in ('points','faces','owner','neighbour','boundary')},
      native_initial_C_Vc_sha256={name:sha(case/'0'/name) for name in ('C','Vc')})
    json_new(case/'provenance/preflight-validation.json',receipt)
    return receipt
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--case',required=True,type=Path)
    group=parser.add_mutually_exclusive_group(required=True);group.add_argument('--mesh',type=Path);group.add_argument('--validate-native',action='store_true')
    args=parser.parse_args()
    print(json.dumps(prepare(args.mesh,args.case) if args.mesh else validate_native(args.case),indent=2))
if __name__=='__main__':
    main()

