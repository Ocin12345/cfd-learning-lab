"""SOURCE DRAFT: separately authorized five-step frozen-flow tracer qualification; no production mode."""
from pathlib import Path
import argparse, hashlib, json, math, os, re, shutil, sys
ROOT=Path('E:/CFD-Learning-Lab' if os.name=='nt' else '/mnt/e/CFD-Learning-Lab')
PRIVATE=ROOT/'local/prospective-tracer-v001'
SCRIPTS=ROOT/'repo/cfd-learning-lab/scripts'
RUN=ROOT/'native/runs/2026-11-02-tracer-tjunction-v001'
ANCESTOR=ROOT/'native/runs/2026-11-01-tjunction-v001/openfoam'
QUAL_SHA='dccbb17e813f2552298332dc3fd7ff713e30df71a70508316a30831995ccd312'
HELPERS={'prepare_openfoam_tjunction_v001.py':'16ded13dcc9d8e7e315f5aaefe24d84af2963e910efda6b48c0d355fe04965d4','prepare_openfoam_heated_channel_pilot_v003.py':'807c30de1c7c7400ab38d2b4592896db774be4c98ec2b853afd2f9fff2ed90b0'}
SOURCES=('prepare_openfoam_tracer_qualification_v004.py','run-openfoam-tracer-qualification-v004.sh','Launch-OpenFOAM-TracerQualification-v004.ps1')
KEYS=('preparation_sha256','wrapper_sha256','launcher_sha256')
CASES={'Control':'openfoam-control-v003','Pulse':'openfoam-pulse-pilot-v003'}
DT=.025

def require(value,message):
    if not value: raise ValueError(message)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()
def readj(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def write_new(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='ascii',newline='\n') as f: f.write(text.rstrip()+'\n')
def json_new(path,value): write_new(path,json.dumps(value,indent=2,allow_nan=False))
def copy_new(src,dst):
    before=sha(src); require(not dst.exists(),'Refuse overwrite '+str(dst)); dst.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(src,dst); require(sha(src)==sha(dst)==before,'Source/copy hash changed '+str(src))
def project_path(raw):
    text=str(raw).replace('\\','/')
    for prefix in ('E:/CFD-Learning-Lab','/mnt/e/CFD-Learning-Lab'):
        if text.startswith(prefix+'/'):
            path=(ROOT/text[len(prefix)+1:]).resolve(); require(ROOT.resolve() in path.parents,'Path escapes project'); return path
    raise ValueError('Expected exact E project prefix '+text)
for name,digest in HELPERS.items(): require(sha(SCRIPTS/name)==digest,'Pinned native parser changed '+name)
sys.path.insert(0,str(SCRIPTS))
import prepare_openfoam_tjunction_v001 as old
import prepare_openfoam_heated_channel_pilot_v003 as scalar_helper

LINEAR_SOURCE_SHA256={'src/OpenFOAM/matrices/lduMatrix/lduMatrix/lduMatrixSolver.C': '95eda050d2718403b71ad92df5ed17440bff88928c8e2996a44b3da4a1824e03', 'src/OpenFOAM/matrices/lduMatrix/solvers/PBiCGStab/PBiCGStab.C': '49b7698e2ce6bd4d4022776548d4af1b3c62c736abe36bd46a92f15f70a49d7b', 'src/OpenFOAM/matrices/lduMatrix/solvers/smoothSolver/smoothSolver.C': '4abaa82c47b08b4c01b787fab8c97559c706be5ec82393a11a6ece091e26454d', 'src/OpenFOAM/matrices/lduMatrix/preconditioners/DILUPreconditioner/DILUPreconditioner.C': '7bda547d004e84e5cbf7e48f459c3f718e15427761254ebbd18d5bece15facfb'}

INSTALLED=(
 'applications/modules/functions/functions.C',
 'src/functionObjects/solvers/scalarTransport/scalarTransport.C',
 'src/functionObjects/solvers/scalarTransport/scalarTransport.H',
 'src/OpenFOAM/primitives/functions/Function1/Table/Table.C',
 'src/OpenFOAM/interpolations/interpolationWeights/stepInterpolationWeights/stepInterpolationWeights.C',
 'src/OpenFOAM/containers/Lists/ListOps/ListOpsTemplates.C',
 'src/finiteVolume/fields/fvPatchFields/derived/uniformFixedValue/uniformFixedValueFvPatchField.C',
 'src/finiteVolume/fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.C',
 'src/finiteVolume/interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H',
 'src/finiteVolume/finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.C',
 'src/finiteVolume/finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.C',
 'src/functionObjects/utilities/writeObjects/writeObjects.C',
 'src/functionObjects/field/writeCellVolumes/writeCellVolumes.C',
 'src/functionObjects/field/writeCellCentres/writeCellCentres.C',
 'src/finiteVolume/fvMatrices/fvMatrix/fvMatrixSolve.C',
)+tuple(LINEAR_SOURCE_SHA256)
def field(path,count,components,dimensions):
    # Original pinned parser with only Foundation14 dimensionless [] normalization.
    # Native files are never rewritten. [] is not accepted for nonzero dimensions.
    text=scalar_helper.clean(Path(path).read_text(encoding='ascii'))
    require(re.search(r'\bformat\s+ascii\s*;',text),'ASCII native field required')
    d=re.search(r'\bdimensions\s*\[([^]]*)\]\s*;',text)
    require(d is not None,'Native field dimensions missing')
    native_dimensions=list(map(int,d[1].split())) if d[1].strip() else [0]*7
    require(len(native_dimensions)==7 and native_dimensions==dimensions,'Native field dimensions differ')
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
def flat(values):
    return [v for row in values for v in row] if values and isinstance(values[0],(list,tuple)) else list(values)
def difference(a,b):
    a,b=flat(a),flat(b); require(len(a)==len(b),'Native vector length differs')
    return max((abs(x-y) for x,y in zip(a,b)),default=0.)
def patch_value(text,count):
    m=re.search(r'\bvalue\s+uniform\s+([^;]+);',text)
    if m:
        value=float(m[1]); require(math.isfinite(value),'Nonfinite boundary scalar'); return [value]*count
    m=re.search(r'\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
    require(m is not None and int(m[1])==count,'Native scalar patch count differs')
    data=list(map(float,m[2].split())); require(len(data)==count and all(map(math.isfinite,data)),'Invalid scalar patch values'); return data

def authorization(path,phase):
    path=Path(path).resolve(); require(PRIVATE.resolve() in path.parents,'Authorization outside candidate private folder')
    a=readj(path)
    require(a.get('approved') is True and a.get('approved_by')=='root' and a.get('action')=='run_tracer_openfoam_qualification_v004','Exact root authorization required')
    require(a.get('allowed_phase')==phase and a.get('actual_steps_authorized')==5 and a.get('final_time_s')==.125 and a.get('production_admitted') is False,'Only exact separate five-step phase allowed')
    require(a.get('fluent_sessions_closed') is True and a.get('native_flow_qualification_sha256')==QUAL_SHA,'Closed Fluent and exact qualified ancestry required')
    require(a.get('constant_field_tolerance')==1e-6 and a.get('scalar_bounds_tolerance')==1e-8 and a.get('final_matrix_residual_target')==1e-8,'Predeclared qualification tolerances required')
    for key,name in zip(KEYS,SOURCES): require(a.get(key)==sha(PRIVATE/name),'Reviewed source differs '+name)
    for name,digest in HELPERS.items(): require(a.get('helper_sha256',{}).get(name)==digest,'Authorization must bind native parser '+name)
    require(set(a.get('installed_source_sha256',{}))==set(INSTALLED) and all(re.fullmatch('[0-9a-f]{64}',v) for v in a['installed_source_sha256'].values()),'Exact reviewed current installed-source hashes required')
    require(all(a['installed_source_sha256'].get(name)==digest for name,digest in LINEAR_SOURCE_SHA256.items()),'Actual captured linear solver source pins required')
    if phase=='Pulse': require(a.get('pulse_endpoint_source_fixture_passed') is True and a.get('pulse_continuum_edge_s')==2.0125 and a.get('pulse_step80_value')==1 and a.get('pulse_step81_value')==0,'Pulse endpoint source fixture required; actual80/81 remain future')
    return path,a

def verify_sources():
    qpath=PRIVATE/'native-flow-qualification-v001.json'; require(sha(qpath)==QUAL_SHA,'Qualified ancestry changed'); q=readj(qpath)
    require(q['source_bytes_unchanged_during_qualification'] is True and q['old_inputs_modified'] is False,'Original qualification did not preserve inputs')
    rows=[]
    for row in q['source_manifest']:
        require(row['sha256_before']==row['sha256_after'],'Original source changed during qualification')
        text=str(row['path']).replace('\\','/')
        # D: Fluent headers are unrelated to this OpenFOAM phase and inaccessible in its WSL prefix.
        if os.name!='nt' and not text.startswith('E:/CFD-Learning-Lab/'): continue
        p=Path(row['path']) if os.name=='nt' else project_path(row['path'])
        require(sha(p)==row['sha256_after'],'Qualified source now changed '+str(p)); rows.append(dict(path=str(p),sha256=sha(p)))
    checks=readj(ROOT/'local/tjunction-v001/openfoam-checks-v001/checks.json')
    require(checks['all_declared_checks_passed'] is True and len(checks['check_results'])==17 and all(v is True for v in checks['check_results'].values()),'Accepted 17-check T-junction flow required')
    return rows

def header(name,kind='dictionary'): return old.header(name,kind)
def frozen_function():
    return '''frozenRuntimeFields
{
 type writeObjects;
 libs ("libutilityFunctionObjects.so");
 objects (U p phi);
 writeOption anyWrite;
 executeControl timeStep; executeInterval 1;
 writeControl timeStep; writeInterval 1;
}
'''
def dictionaries(phase):
    require(phase in CASES,'Unknown qualification phase')
    c1=1 if phase=='Control' else 0
    top='type fixedValue; value uniform 1;' if phase=='Control' else '''type uniformFixedValue;
 uniformValue { type table; values ((0 1) (2.0125 0) (20 0)); interpolationScheme step; outOfBounds clamp; }
 value uniform 1;'''
    c=header('tracer','volScalarField')+f'''dimensions [0 0 0 0 0 0 0];
internalField uniform {c1};
boundaryField
{{
 inlet_left {{ type fixedValue; value uniform {c1}; }}
 inlet_top {{ {top} }}
 outlet {{ type zeroGradient; }}
 walls {{ type zeroGradient; }}
 frontAndBackPlanes {{ type empty; }}
}}
'''
    base=old.dictionaries()
    data={'system/decomposeParDict':base['system/decomposeParDict'],'0/tracer':c,
      'system/controlDict':header('controlDict')+'''solver functions;
subSolver incompressibleFluid;
subSolverTime 2000;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 0.125;
deltaT 0.025;
adjustTimeStep false;
writeControl timeStep;
writeInterval 1;
purgeWrite 0;
writeFormat ascii;
writePrecision 17;
writeCompression off;
timeFormat general;
timePrecision 10;
runTimeModifiable false;
''',
      'system/fvSchemes':header('fvSchemes')+'''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,tracer) Gauss upwind; div(phi,U) bounded Gauss linearUpwind grad(U); div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear orthogonal; }
interpolationSchemes { default linear; }
snGradSchemes { default orthogonal; }
''',
      'system/fvSolution':header('fvSolution')+'''solvers
{
 tracer { solver PBiCGStab; preconditioner DILU; tolerance 1e-12; relTol 0; maxIter 1000; nCorrectors 0; }
 tracerFinal { $tracer; }
 p { solver GAMG; tolerance 1e-10; relTol 0.01; smoother GaussSeidel; }
 U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01; }
}
PIMPLE { nOuterCorrectors 1; nCorrectors 0; nNonOrthogonalCorrectors 0; }
relaxationFactors { equations { tracer 1; } }
''',
      'system/functions':header('functions')+'''transportTracer
{
 type scalarTransport;
 libs ("libsolverFunctionObjects.so");
 field tracer; phi phi;
 diffusivity constant; D 0.001;
 executeControl timeStep; executeInterval 1;
 writeControl timeStep; writeInterval 1;
}
'''+frozen_function(),
      'system/baselineDiagnostics':header('baselineDiagnostics')+'''nativeCellCentres { type writeCellCentres; libs ("libfieldFunctionObjects.so"); }
nativeCellVolumes { type writeCellVolumes; libs ("libfieldFunctionObjects.so"); }
'''+frozen_function()}
    return data

def checked_case(phase): return RUN/CASES[phase]
def prepare(phase,authpath):
    authpath,a=authorization(authpath,phase); watched=verify_sources(); case=checked_case(phase)
    require(not case.exists(),'Fresh versioned case required'); require(shutil.disk_usage(ROOT).free>=10*1024**3,'Need 10GiB E disk')
    poly=old.native_poly(ANCESTOR); old.native_patches(poly); old.native_geometry(poly)
    for name,text in dictionaries(phase).items(): write_new(case/name,text)
    for name in ('points','faces','owner','neighbour','boundary'): copy_new(poly/name,case/'0/polyMesh'/name)
    for name in ('U','p','phi'):
        for folder in ('0','2000','provenance/baseline-native-inputs'): copy_new(ANCESTOR/'2000'/name,case/folder/name)
    copy_new(ANCESTOR/'2000/C',case/'provenance/ancestor-C')
    for name in ('physicalProperties','momentumTransport'): copy_new(ANCESTOR/'constant'/name,case/'constant'/name)
    for name in SOURCES: copy_new(PRIVATE/name,case/'provenance'/name)
    for name in HELPERS: copy_new(SCRIPTS/name,case/'provenance'/name)
    copy_new(authpath,case/'provenance/root-run-authorization.json'); copy_new(PRIVATE/'native-flow-qualification-v001.json',case/'provenance/native-flow-qualification-v001.json')
    copy_new(case/'0/tracer',case/'provenance/initial-tracer')
    json_new(case/'provenance/original-source-snapshot.json',dict(sources=watched,qualification_sha256=QUAL_SHA))
    immutable={p.relative_to(case).as_posix():sha(p) for p in case.rglob('*') if p.is_file() and not p.relative_to(case).as_posix().startswith('0/')}
    json_new(case/'provenance/immutable-input-sha256.json',immutable)
    json_new(case/'run-status.json',dict(stage='prepared_only',phase=phase,actual_steps_authorized=5,production_admitted=False,solver_started=False))
    return dict(case=str(case),phase=phase,prepared_only=True,solver_started=False,qualification_sha256=QUAL_SHA)

def blocks(raw):
    text=old.clean(raw); m=re.search(r'\bboundaryField\s*\{',text); require(m is not None,'Missing boundaryField'); pos=m.end(); rows={}
    while True:
        while pos<len(text) and text[pos].isspace(): pos+=1
        if text[pos]=='}': break
        m=re.match(r'([A-Za-z_][A-Za-z_0-9]*)\s*\{',text[pos:]); require(m is not None,'Invalid native boundary name'); name=m[1]; start=pos+m.end(); end=start; depth=1
        while depth and end<len(text):
            depth+=(text[end]=='{')-(text[end]=='}'); end+=1
        require(depth==0 and name not in rows,'Invalid/duplicate boundary block'); rows[name]=text[start:end-1]; pos=end
    return rows

def all_phi(path,poly):
    nfaces,_=old.counted(poly/'faces'); ninternal,_=old.counted(poly/'neighbour'); n,body=old.counted(poly/'boundary')
    patches={name:dict(kind=re.search(r'\btype\s+(\w+)\s*;',b)[1],count=int(re.search(r'\bnFaces\s+(\d+)\s*;',b)[1]),start=int(re.search(r'\bstartFace\s+(\d+)\s*;',b)[1])) for name,b in re.findall(r'([A-Za-z_][A-Za-z_0-9]*)\s*\{([^{}]*)\}',body)}
    require(len(patches)==n,'Rank/native boundary count differs'); values=[0.]*nfaces; covered=[False]*nfaces; covered[:ninternal]=[True]*ninternal
    values[:ninternal]=field(path,ninternal,1,[0,3,-1,0,0,0,0]); bs=blocks(path.read_text()); require(set(bs)==set(patches),'phi boundary coverage differs')
    for name,p in patches.items():
        lo,hi=p['start'],p['start']+p['count']; require(0<=lo<=hi<=nfaces and not any(covered[lo:hi]),'Invalid/overlapping native patch'); covered[lo:hi]=[True]*(hi-lo)
        if p['kind']!='empty': values[lo:hi]=patch_value(bs[name],p['count'])
    require(all(covered),'Missing native face coverage'); return values

def verify(phase,authpath,completed=False):
    _,auth=authorization(authpath,phase); verify_sources(); case=checked_case(phase)
    installed=readj(case/'provenance/installed-source-sha256.json'); require(installed==auth['installed_source_sha256'],'Captured installed-source selection changed')
    for name,digest in installed.items(): require(sha(case/'provenance/installed-sources'/name)==digest,'Captured installed source changed '+name)
    for name,digest in readj(case/'provenance/immutable-input-sha256.json').items(): require(sha(case/name)==digest,'Immutable prepared input changed '+name)
    for name,text in dictionaries(phase).items():
        if name!='0/tracer': require((case/name).read_text()==text,'Reviewed dictionary changed '+name)
    poly=old.native_poly(case); original_poly=old.native_poly(ANCESTOR)
    for name in ('points','faces','owner','neighbour','boundary'): require(sha(poly/name)==sha(original_poly/name),'Frozen native mesh copy changed '+name)
    patches=old.native_patches(poly); geom=old.native_geometry(poly); volumes=geom['volumes']
    old.mesh_quality((case/'mesh-check.log').read_text())
    require(difference(field(case/'0/Vc',4400,1,[0,3,0,0,0,0,0]),volumes)<=1e-14,'Actual native Vc differs from oriented mesh')
    require(difference(field(case/'0/C',4400,3,[0,1,0,0,0,0,0]),field(case/'provenance/ancestor-C',4400,3,[0,1,0,0,0,0,0]))<=1e-12,'Native C/source ordering differs')
    dims={'U':(3,[0,1,-1,0,0,0,0]),'p':(1,[0,2,-2,0,0,0,0])}; baseline={n:field(case/'provenance/baseline-native-inputs'/n,4400,k,d) for n,(k,d) in dims.items()}
    q0=all_phi(case/'provenance/baseline-native-inputs/phi',poly); outlet=patches['outlet']; require(all(v>=0 for v in q0[outlet['start']:outlet['start']+outlet['count']]),'Frozen outlet backflow is unsupported in qualification')
    indices=range(6) if completed else (0,); history=[]; scalar_states={}
    for i in indices:
        label=f'{i*DT:.10g}'; folder=case/label
        for n,(k,d) in dims.items(): require(difference(field(folder/n,4400,k,d),baseline[n])<=1e-12,'Frozen '+n+' changed at'+label)
        require(difference(all_phi(folder/'phi',poly),q0)<=1e-12,'Frozen native phi changed at'+label)
        raw=(folder/'tracer').read_text(); require(re.search(r'\bclass\s+volScalarField\s*;',raw) and re.search(r'\bobject\s+tracer\s*;',raw),'Expected native tracer scalar object'); c=field(folder/'tracer',4400,1,[0,0,0,0,0,0,0]); scalar_states[i]=c; bs=blocks(raw); require(set(bs)==set(patches),'Scalar patch coverage differs')
        for name in ('outlet','walls'): require(re.search(r'\btype\s+zeroGradient\s*;',bs[name]),'Scalar zeroGradient changed '+name)
        require(re.search(r'\btype\s+empty\s*;',bs[old.EMPTY]),'Scalar empty planes changed')
        require(re.search(r'\btype\s+fixedValue\s*;',bs['inlet_left']),'Left scalar inlet type changed')
        require(re.search(r'\btype\s+'+('fixedValue' if phase=='Control' else 'uniformFixedValue')+r'\s*;',bs['inlet_top']),'Top scalar inlet type changed')
        for name,expected in [('inlet_left',1 if phase=='Control' else 0),('inlet_top',1)]: require(max(abs(v-expected) for v in patch_value(bs[name],20))<=1e-12,'Actual scalar inlet value changed '+name)
        if i==0: require(max(abs(v-(1 if phase=='Control' else 0)) for v in c)<=1e-12,'Initial scalar differs')
        elif phase=='Control': require(max(abs(v-1) for v in c)<=1e-6,'Constant-c1 field tolerance missed')
        else: require(min(c)>=-1e-8 and max(c)<=1+1e-8 and max(map(abs,c))>0,'Pulse scalar bounds/no-change failure')
        inventory=sum(v*x for v,x in zip(volumes,c))/.01
        if phase=='Control': require(abs(inventory/.44-1)<=1e-6,'Constant-c1 inventory tolerance missed')
        history.append(dict(index=i,time_s=i*DT,min_c=min(c),max_c=max(c),inventory_kg_per_m=inventory))
    if completed:
        log=(case/'solver.log').read_text(); require(re.search(r'^End\s*$',log,re.M) is not None and not re.search(r'FOAM FATAL|Floating point exception|Segmentation fault',log),'Native solver did not close cleanly')
        times=[float(v) for v in re.findall(r'^Time = ('+old.FLOAT+r')s?\s*$',log,re.M)]; require(len(times)==5 and max(abs(t-i*DT) for i,t in enumerate(times,1))<=1e-12,'Exactly five native steps required')
        residuals=re.findall(r'Solving for tracer, Initial residual = ('+old.FLOAT+r'), Final residual = ('+old.FLOAT+r'), No Iterations (\d+)',log)
        require(len(residuals)==5 and all(float(v[1])<=1e-8 for v in residuals),'Five actual scalar solves/final residual gate required'); require(not re.search(r'Solving for (?:Ux|Uy|Uz|p),',log),'Unexpected frozen-flow solve')
        rank_cells=[]
        for rank in (0,1):
            rp=case/('processor'+str(rank)); rpoly=old.native_poly(rp); nc,body=old.counted(rpoly/'cellProcAddressing'); cellmap=list(map(int,body.split())); require(len(cellmap)==nc and len(set(cellmap))==nc and all(0<=v<4400 for v in cellmap),'Invalid rank cell addressing'); rank_cells+=cellmap
            nf,body=old.counted(rpoly/'faceProcAddressing'); facemap=list(map(int,body.split())); require(len(facemap)==nf and all(v!=0 for v in facemap),'Invalid signed face addressing')
            globalfaces=[abs(v)-1 for v in facemap]; require(all(v<len(q0) for v in globalfaces),'Rank face IDs outside global mesh')
            require(nf==old.counted(rpoly/'faces')[0],'Rank face addressing length differs from mesh')
            npoints,pointbody=old.counted(rpoly/'pointProcAddressing'); pointmap=list(map(int,pointbody.split())); require(len(pointmap)==npoints and len(set(pointmap))==npoints and all(0<=v<9282 for v in pointmap),'Invalid rank point addressing')
            native_npoints,pointbody=old.counted(rpoly/'points'); coords=[list(map(float,v.split())) for v in re.findall(r'\(([^()]*)\)',pointbody)]; require(native_npoints==npoints and difference(coords,[geom['points'][g] for g in pointmap])<=1e-12,'Rank/source native points differ')
            _,facebody=old.counted(rpoly/'faces'); rankfaces=[list(map(int,v.split())) for size,v in re.findall(r'(\d+)\s*\(([^()]*)\)',facebody)]; require(len(rankfaces)==nf,'Rank native faces differ')
            for localface,globalface in zip(rankfaces,globalfaces): require(len(localface)==4 and {pointmap[v] for v in localface}==set(geom['faces'][globalface]),'Rank face/source vertex mapping differs')
            expected_phi=[q0[g]*(1 if signed>0 else -1) for g,signed in zip(globalfaces,facemap)]
            for i in range(6):
                label=f'{i*DT:.10g}'; folder=rp/label
                for n,(k,d) in dims.items(): require(difference(field(folder/n,nc,k,d),[baseline[n][g] for g in cellmap])<=1e-12,'Actual rank frozen '+n+' changed')
                require(difference(all_phi(folder/'phi',rpoly),expected_phi)<=1e-12,'Actual rank frozen all-face phi changed')
                require(difference(field(folder/'tracer',nc,1,[0,0,0,0,0,0,0]),[scalar_states[i][g] for g in cellmap])<=1e-12,'Rank/reconstructed c differs')
                if i:
                    text=(folder/'uniform/time').read_text(); value=re.search(r'\bvalue\s+('+old.FLOAT+r')\s*;',text); index=re.search(r'\bindex\s+(\d+)\s*;',text); delta=re.search(r'\bdeltaT\s+('+old.FLOAT+r')\s*;',text)
                    require(value and index and delta and abs(float(value[1])-i*DT)<=1e-12 and int(index[1])==i and abs(float(delta[1])-DT)<=1e-12,'Rank native time/index/dt differs')
        require(sorted(rank_cells)==list(range(4400)),'Two ranks do not biject to global cells')
    verify_sources()
    return dict(phase=phase,completed_five_steps=completed,states_verified=len(history),frozen_U_p_all_phi_unchanged=True,native_Vc_geometry_agree=True,history=history,full_signed_scalar_budget_analysis_pending=True,production_admitted=False,physical_validation=False)

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--phase',choices=tuple(CASES),required=True); p.add_argument('--authorization',type=Path,required=True)
    mode=p.add_mutually_exclusive_group(required=True); mode.add_argument('--prepare',action='store_true'); mode.add_argument('--verify-baseline',action='store_true'); mode.add_argument('--verify-completed',action='store_true'); mode.add_argument('--verify-sources',action='store_true')
    a=p.parse_args(); authorization(a.authorization,a.phase)
    result=prepare(a.phase,a.authorization) if a.prepare else dict(verified_source_files=len(verify_sources())) if a.verify_sources else verify(a.phase,a.authorization,a.verify_completed)
    print(json.dumps(result,indent=2,allow_nan=False))
if __name__=='__main__': main()