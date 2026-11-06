"""Prepare/check the fresh ideal-gas inviscid 2D nozzle. Starts no solver."""
from pathlib import Path
import argparse, hashlib, json, math, re, shutil
RUN='2026-11-01-nozzle-v001'
F=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
ROOT=Path(__file__).resolve().parents[3]
SOURCES={
'module':'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/applications/modules/shockFluid/shockFluid.C',
'LTS':'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/applications/modules/shockFluid/setRDeltaT.C',
'schemes':'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/tutorials/shockFluid/forwardStep/system/fvSchemes',
'gas':'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/tutorials/shockFluid/forwardStep/constant/physicalProperties',
'total_pressure':'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/src/finiteVolume/fields/fvPatchFields/derived/totalPressure/totalPressureFvPatchScalarField.H',
'pressure_velocity':'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/src/finiteVolume/fields/fvPatchFields/derived/pressureInletOutletVelocity/pressureInletOutletVelocityFvPatchVectorField.H',
}
def require(v,m):
 if not v: raise ValueError(m)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(1048576),b''): h.update(b)
 return h.hexdigest()
def write(p,s):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x',encoding='ascii',newline='\n') as f:f.write(s.rstrip()+'\n')
def jwrite(p,v):write(p,json.dumps(v,indent=2,allow_nan=False))
def copy(a,b):
 before=sha(a); b.parent.mkdir(parents=True,exist_ok=True)
 require(not b.exists(),'Refuse overwrite '+str(b)); shutil.copy2(a,b)
 require(sha(a)==before==sha(b),'Copy hash differs '+str(b))
def header(name,kind='dictionary'):
 return 'FoamFile\n{\n version 2.0;\n format ascii;\n class '+kind+';\n object '+name+';\n}\n\n'
def dictionaries():
 control='''solver shockFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 6000;
deltaT 1;
writeControl timeStep;
writeInterval 500;
purgeWrite 0;
writeFormat ascii;
writePrecision 14;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable false;
adjustTimeStep false;
functions
{
'''
 for patch in ('inlet','outlet','upper_wall','lower_wall'):
  control+=''' flux_%s
 {
  type surfaceFieldValue;
  libs ("libfieldFunctionObjects.so");
  patch %s;
  operation sum;
  fields (phi);
  writeFields false;
  log false;
  writeControl timeStep;
  writeInterval 1;
 }
'''%(patch,patch)
 control+='}\n'
 content={
'system/controlDict':control,
'system/fvSchemes':'''fluxScheme Kurganov;
ddtSchemes { default localEuler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes
{
 default linear;
 reconstruct(rho) vanLeer;
 reconstruct(U) vanLeerV;
 reconstruct(T) vanLeer;
}
snGradSchemes { default corrected; }
''',
'system/fvSolution':'''solvers
{
 "(rho|rhoU|rhoE).*" { solver diagonal; }
 "U.*" { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-9; relTol 0.01; }
 "e.*" { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-10; relTol 0; }
}
PIMPLE { maxCo 0.4; rDeltaTSmoothingCoeff 0.02; }
''',
'system/decomposeParDict':'''numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
''',
'constant/physicalProperties':'''thermoType
{
 type hePsiThermo;
 mixture pureMixture;
 transport const;
 thermo hConst;
 equationOfState perfectGas;
 specie specie;
 energy sensibleInternalEnergy;
}
mixture
{
 specie { molWeight %.16g; }
 thermodynamics { Cp 1004.675; hf 0; }
 transport { mu 0; Pr 1; }
}
'''%(8314.46261815324/287.05),
'0/U':'''dimensions [0 1 -1 0 0 0 0];
internalField uniform (50 0 0);
boundaryField
{
 inlet { type pressureInletOutletVelocity; value uniform (50 0 0); }
 outlet { type zeroGradient; }
 upper_wall { type slip; }
 lower_wall { type slip; }
 "(frontAndBack.*|defaultFaces)" { type empty; }
}
''',
'0/p':'''dimensions [1 -1 -2 0 0 0 0];
internalField uniform 140000;
boundaryField
{
 inlet { type totalPressure; rho none; psi psi; gamma 1.4; p0 uniform 200000; value uniform 190000; }
 outlet { type fixedValue; value uniform 140000; }
 upper_wall { type zeroGradient; }
 lower_wall { type zeroGradient; }
 "(frontAndBack.*|defaultFaces)" { type empty; }
}
''',
'0/T':'''dimensions [0 0 0 1 0 0 0];
internalField uniform 300;
boundaryField
{
 inlet { type totalTemperature; psi psi; gamma 1.4; T0 uniform 300; value uniform 296.8; }
 outlet { type zeroGradient; }
 upper_wall { type zeroGradient; }
 lower_wall { type zeroGradient; }
 "(frontAndBack.*|defaultFaces)" { type empty; }
}
'''}
 return {n:header(Path(n).name,'volVectorField' if n=='0/U' else 'volScalarField' if n in ('0/p','0/T') else 'dictionary')+b for n,b in content.items()}
def checked_case(case):
 case=Path(case).resolve()
 require(case==(ROOT/'native/runs'/RUN/'openfoam').resolve(),'Only exact dedicated nozzle case allowed')
 return case
def clean(s):return re.sub(r'/\*.*?\*/|//[^\n]*','',s,flags=re.S)
def counted(path):
 s=clean(path.read_text(encoding='ascii')); h=re.search(r'FoamFile\s*\{([^{}]*)\}',s)
 require(h is not None and re.search(r'\bformat\s+ascii\s*;',h[1]),'Require ASCII native list')
 m=re.fullmatch(r'\s*(\d+)\s*\((.*)\)\s*',s[h.end():],re.S)
 require(m is not None,'Invalid counted list '+str(path));return int(m[1]),m[2]
def poly(case):
 candidates=[case/p for p in ('0/polyMesh','constant/polyMesh') if (case/p/'points').is_file()]
 require(len(candidates)==1,'Require exactly one native initial polyMesh');return candidates[0]
def patches(path):
 n,b=counted(path/'boundary'); rows={}
 for name,body in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',b):
  rows[name.strip('"')]=dict(kind=re.search(r'\btype\s+(\w+)\s*;',body)[1],count=int(re.search(r'\bnFaces\s+(\d+)\s*;',body)[1]),start=int(re.search(r'\bstartFace\s+(\d+)\s*;',body)[1]))
 require(n==len(rows) and set(rows)=={'inlet','outlet','upper_wall','lower_wall','frontAndBackPlanes'},'Unexpected converted patch set')
 require(rows['upper_wall']['kind']=='wall' and rows['lower_wall']['kind']=='wall' and rows['frontAndBackPlanes']['kind']=='empty','Require slipwall geometrical wall and empty extrudedplanes')
 require([rows[k]['count'] for k in ('inlet','outlet','upper_wall','lower_wall','frontAndBackPlanes')]==[40,40,240,240,19200],'Converted patch facecounts differ')
 return rows
def mesh_quality(log):
 require('checkMesh -allGeometry -allTopology' in log and 'Mesh OK.' in log and not re.search(r'Failed\s+\d+\s+mesh checks|FOAM FATAL',log),'Native full checkMesh must pass')
 patterns={'cells':r'\bcells:\s+(\d+)','minimum_volume_m3':r'Min volume\s*=\s*('+F+')','maximum_nonorthogonality_degrees':r'Mesh non-orthogonality Max:\s*('+F+')','maximum_skewness':r'Max skewness\s*=\s*('+F+')'}
 q={}
 for k,p in patterns.items():
  m=re.search(p,log); require(m is not None,'Missing mesh quality '+k);q[k]=int(m[1]) if k=='cells' else float(m[1])
 require(q['cells']==9600 and q['minimum_volume_m3']>0 and q['maximum_nonorthogonality_degrees']<=75 and q['maximum_skewness']<=4,'Mesh numerical quality gate missed');return q
def prepare(mesh,case):
 mesh=Path(mesh).resolve();case=checked_case(case)
 require(not case.exists(),'Refuse existing native case')
 require(mesh==case.parent/'mesh/nozzle.msh','Require selected source mesh/nozzle.msh')
 raw=mesh.read_bytes();require(re.search(rb'\(2\s+2\)',raw),'Require 2D FluentASCII')
 m=re.search(rb'^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)',raw,re.M)
 require(m is not None and int(m[1],16)==9600,'Require exact9600 sharedquads')
 for name in ('fluid','interior','upper_wall','lower_wall','inlet','outlet'):
  require(re.search(rb'\(45\s+\([^)]*\s'+name.encode()+rb'\)\(\)\)',raw),'Missing named zone '+name)
 require(shutil.disk_usage(case.anchor).free>=10*1024**3,'Require10GiBdisk')
 setup=json.loads((mesh.parent/'setup.json').read_text(encoding='utf-8-sig'))
 require(setup['cells']==9600 and setup['mesh_sha256']==sha(mesh),'Mesh setup/hash differs')
 case.mkdir(parents=True)
 for n,s in dictionaries().items():write(case/n,s)
 copy(mesh,case/'input-nozzle.msh')
 for src in (Path(__file__),ROOT/'local/nozzle-v001/run-openfoam.sh',ROOT/'local/nozzle-v001/Launch-OpenFOAM.ps1',mesh.parent/'setup.json'):
  copy(src,case/'provenance'/src.name)
 target=dict(model='2D ideal-gas inviscid Euler converging-diverging nozzle',full_height_inlet_m=.2,full_height_throat_m=.08,full_height_exit_m=.16,length_m=1,throat_x_m=.4,extrusion_depth_m=.01,nx=240,ny=40,intended_cells=9600,R_J_kg_K=287.05,gamma=1.4,Cp_J_kg_K=1004.675,molecular_weight_kg_kmol=8314.46261815324/287.05,mu_Pa_s=0,inlet_total_p_Pa=200000,inlet_total_T_K=300,outlet_static_p_Pa=140000,solver='OpenFOAM Foundation14 shockFluid',flux_scheme='Kurganov',ddt_scheme='localEuler',maxCo=.4,max_pseudo_iterations=6000,saved_pseudo_iterations=list(range(500,6001,500)),last_three_saved_states=[5000,5500,6000],mpi_processes=2,initial_U_m_s=[50,0,0],initial_p_Pa=140000,initial_T_K=300,mass_flux_interval_iterations=1,relative_mass_imbalance_target=1e-4,last_three_relative_L2_U_p_T_target=1e-3,last_three_shock_location_range_over_length_target=.01,diagonal_residuals_are_not_convergence_evidence=True,iteration_indices_are_not_physical_time=True,steady_field_tracers_are_not_transient_CFD=True,pressure_note='Absolute staticp in Pa; phi is kg/s through0.01m extruded depth; divide by0.01 forkg/s/m',mesh_source_path=str(mesh),mesh_source_sha256=sha(mesh),physical_validation=False,mesh_independence_established=False,official_sources=SOURCES)
 jwrite(case/'declared-comparison.json',target)
 jwrite(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file()})
 jwrite(case/'run-status.json',dict(stage='prepared_only',solver_started=False))
 return dict(case=str(case),cells=9600,solver_started=False,mesh_sha256=sha(mesh))
def validate(case):
 case=checked_case(case);require(not (case/'solver.log').exists() and not list(case.glob('processor*')),'Preflight must precede any solver/decomposition')
 manifest=json.loads((case/'provenance/prepared-input-sha256.json').read_text())
 require(all(sha(case/n)==v for n,v in manifest.items()),'Prepared input source changed')
 p=poly(case);rows=patches(p);q=mesh_quality((case/'mesh-check.log').read_text())
 n,b=counted(p/'points');pts=[tuple(map(float,v.split())) for v in re.findall(r'\(([^()]*)\)',b)]
 require(n==len(pts) and all(len(v)==3 and all(map(math.isfinite,v)) for v in pts),'Invalid nativepoints')
 bounds=[[min(v[j] for v in pts),max(v[j] for v in pts)] for j in range(3)]
 require(all(abs(a-b)<1e-8 for row,t in zip(bounds,([0,1],[-.1,.1],[-.005,.005])) for a,b in zip(row,t)),'Unexpected native bounds')
 receipt=dict(passed=True,solver_started=False,mesh_location=p.relative_to(case).as_posix(),actual_native_patches=rows,actual_mesh_quality=q,native_bounds_m=bounds,extrusion_depth_m=.01,prepared_inputs_unchanged=True,native_mesh_sha256={n:sha(p/n) for n in ('points','faces','owner','neighbour','boundary')})
 jwrite(case/'provenance/preflight-validation.json',receipt);return receipt
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',required=True,type=Path);g=ap.add_mutually_exclusive_group(required=True);g.add_argument('--mesh',type=Path);g.add_argument('--validate-native',action='store_true');a=ap.parse_args();print(json.dumps(prepare(a.mesh,a.case) if a.mesh else validate(a.case),indent=2))
if __name__=='__main__':main()
