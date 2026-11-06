"""SOURCE ONLY: independent streaming actual tracer CSV/GIF/native-input readback.
No plotting import, image regeneration, CFD, WSL, GUI or mutation of render/native inputs.
"""
from pathlib import Path
import argparse,csv,ctypes,datetime,gzip,hashlib,importlib.util,json,re,sys
import h5py,numpy as np
from PIL import Image
ROOT=Path('E:/CFD-Learning-Lab');PRIVATE=ROOT/'local/prospective-tracer-v001'
GEOMETRY=ROOT/'native/runs/2026-11-01-tjunction-v001/mesh/geometry.npz';GEOMETRY_SHA='671236b7bf7658fcbae5330b4cd6ce66c8e2dee9ff06753ebfc45c0a4cfb2949'
ADAPTER=PRIVATE/'scalar_native_geometry_v003.py';ADAPTER_SHA='06f6b4089116d8e7cdab8b4050c30d2939d3876381ebbaa878c7908e6f2d1903'
FV=PRIVATE/'scalar_fv_audit_v002.py';FV_SHA='75ceec60add2f5c53fb82585d105dcc5cd7c9d8d19e1ea8ac3f2dde14aa455cd'
Q=PRIVATE/'prepare_openfoam_tracer_qualification_v004.py';Q_SHA='278e6948a244e82887d8eb8dc2da90b502a14a7165de2238c3ba1e01b2890888'
HELPERS={'prepare_openfoam_tjunction_v001.py':'16ded13dcc9d8e7e315f5aaefe24d84af2963e910efda6b48c0d355fe04965d4','prepare_openfoam_heated_channel_pilot_v003.py':'807c30de1c7c7400ab38d2b4592896db774be4c98ec2b853afd2f9fff2ed90b0'}
N=4400;STATES=801;DT=.025;NAMES=('fluent','openfoam');STEPS=list(range(0,801,8));PATCHES=('inlet_left','inlet_top','outlet','walls');PREFIX='tracer-tjunction'
FIELD_HEADER=['solver','step','time_s','source_cell_index_zero_based','native_cell_index_zero_based','x_m','y_m','area_m2','c']
HISTORY_HEADER=['solver','step','time_s','inventory_kg_per_m','scalar_min','scalar_max','negative_cell_count','above_one_cell_count','top_inlet_c','outlet_outflow_weighted_c','global_Euler_defect_relative','own_terminal_scalar_residual']+[p+'_'+t+'_inward_kg_m_s' for p in PATCHES for t in ('advective','diffusive','total')]
FLOAT=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'

def require(v,m):
    if not v:raise ValueError(m)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def project_path(raw):
    p=Path(raw);require(p.is_absolute(),'Absolute project path required');p=p.resolve();require(ROOT.resolve() in p.parents,'Project path required');return p
class Evidence:
    def __init__(self):self.rows={};self.cache={}
    def pin(self,ref):
        p=project_path(ref['path']);h=ref['sha256'];require(re.fullmatch('[0-9a-f]{64}',h) is not None,'Explicit lowercase SHA256 required')
        if p not in self.rows:require(p.is_file() and sha(p)==h,'Missing/changed selected input: '+str(p));self.rows[p]=h
        require(self.rows[p]==h,'Conflicting input hashes');return p
    def touch(self,p):return self.pin(dict(path=str(p),sha256=sha(p)))
    def json(self,ref):
        p=self.pin(ref)
        if p not in self.cache:self.cache[p]=json.loads(p.read_text(encoding='utf-8-sig'))
        return self.cache[p]
    def value(self,ref):
        require(isinstance(ref.get('key_path'),list),'Explicit key_path required');v=self.json(ref)
        for k in ref['key_path']:v=v[k]
        return v
    def module(self,name,p,h):
        p=self.pin(dict(path=str(p),sha256=h));require(name not in sys.modules,'Fresh verified module namespace');s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
    def finish(self):
        rows=[]
        for p,h in self.rows.items():
            actual=sha(p);require(actual==h,'Input changed during verification: '+str(p));rows.append(dict(path=str(p),bytes=p.stat().st_size,sha256_before=h,sha256_after=actual))
        return rows

def memory_guard():
    class M(ctypes.Structure):_fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(k,ctypes.c_ulonglong) for k in ('total','available','pageTotal','pageAvailable','virtualTotal','virtualAvailable','extended')]
    m=M();m.length=ctypes.sizeof(m);require(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)) and m.available>=2.5*1024**3,'Fresh read-only2.5GiB RAM floor');return int(m.available)
def load_geometry(e):
    e.pin(dict(path=str(FV),sha256=FV_SHA));module=e.module('independent_tracer_visual_geometry',ADAPTER,ADAPTER_SHA);native,proof=module.qualified_geometry()
    for row in proof['consumed_inputs']:e.pin(row)
    p=e.pin(dict(path=str(GEOMETRY),sha256=GEOMETRY_SHA))
    with np.load(p,allow_pickle=False) as z:g={k:z[k].copy() for k in ('nodes','quads','centres','areas')}
    require(g['quads'].shape==(N,4) and g['centres'].shape==(N,2),'Exact actual canonical geometry shape');return g,native

def load_bundle(name,selected,g,native,e):
    require(selected['full_actual_audit_accepted_by_root'] is True,'Root-accepted full scalar export required')
    for ref in selected['accepted_evidence']:e.pin(ref)
    p=e.pin(selected['fields'])
    with np.load(p,allow_pickle=False) as z:b={k:z[k].copy() for k in ('step','time_s','c','native_cell_to_source_cell','source_native_cell_index','source_centres_m','source_areas_m2')}
    require(b['c'].shape==(STATES,N) and b['c'].dtype==np.float64 and np.isfinite(b['c']).all(),'Actual finite source-order scalar history required')
    require(np.array_equal(b['step'],np.arange(STATES)) and np.issubdtype(b['step'].dtype,np.integer) and b['time_s'].shape==(STATES,) and np.max(abs(b['time_s']-b['step']*DT))<=1e-12,'Actual step/time arrays differ')
    mapping=np.asarray(native[name]['native_cell_to_source_cell']);inverse=np.argsort(mapping);stencil=native[name]['stencil']
    require(np.array_equal(b['native_cell_to_source_cell'],mapping) and np.array_equal(b['source_native_cell_index'],inverse),'Actual source/native mapping differs')
    require(np.array_equal(b['source_centres_m'],stencil['cell_centres'][inverse]) and np.array_equal(b['source_areas_m2'],(stencil['volumes']/stencil['depth_m'])[inverse]),'Exact audited native-derived geometry metadata differs')
    coord=float(np.max(abs(b['source_centres_m']-g['centres'])));area=float(np.max(abs(b['source_areas_m2']-g['areas'])));require(coord<=1e-12 and area<=1e-14,'Existing native/canonical geometry correspondence exceeded')
    rows=e.json(selected['manifest'])['states'];require(len(rows)==STATES,'Actual801-state adapter required')
    for i,row in enumerate(rows):require(row['step']==i and row['time_s']==b['time_s'][i],'Adapter state/time differs')
    # Independent native mass-weighted integral. This preserves the specified native reduction order.
    b['inventory']=np.sum(b['c'][:,mapping]*(native[name]['rho']*stencil['volumes'])[None,:],axis=1)/stencil['depth_m']
    b.update(rows=rows,native_geometry_error_m=coord,native_geometry_area_error_m2=area,name=name)
    return b
def verify_cell_csv(path,bundles,g,full):
    dtype=np.dtype([('solver','U8'),('step','i8'),('time','f8'),('cell','i8'),('native','i8'),('x','f8'),('y','f8'),('area','f8'),('c','f8')]);count=0
    opener=gzip.open if full else open
    with opener(path,'rt',encoding='utf-8',newline='') as stream:
        require(next(csv.reader([stream.readline()]))==FIELD_HEADER,'Actual scalar CSV header differs')
        for name in NAMES:
            b=bundles[name]
            for step in (range(STATES) if full else (800,)):
                rows=np.loadtxt(stream,delimiter=',',dtype=dtype,max_rows=N,comments=None,encoding='utf-8')
                require(rows.shape==(N,) and np.all(rows['solver']==name) and np.all(rows['step']==step) and np.all(rows['time']==b['time_s'][step]),'CSV actual solver/state/time differs')
                require(np.array_equal(rows['cell'],np.arange(N)) and np.array_equal(rows['native'],b['source_native_cell_index']),'CSV source/native cell order differs')
                for key,expected in [('x',g['centres'][:,0]),('y',g['centres'][:,1]),('area',g['areas']),('c',b['c'][step])]:require(np.array_equal(rows[key],expected),'CSV canonical geometry or actual scalar differs: '+name+'/'+str(step)+'/'+key)
                count+=N
                if full and step%100==0:memory_guard();print(json.dumps(dict(streaming_actual_CSV_solver=name,step=step,rows_checked=count)),flush=True)
        require(not stream.read(),'Unexpected extra scalar CSV content')
    require(count==(7048800 if full else 8800),'Complete actual CSV row total differs');return count

def verify_history(path,bundles,e):
    count=0;max_inventory_budget_error=0.
    with path.open(encoding='utf-8',newline='') as stream:
        reader=csv.DictReader(stream);require(reader.fieldnames==HISTORY_HEADER,'Actual history CSV header differs')
        for name in NAMES:
            b=bundles[name]
            for step in range(STATES):
                row=next(reader);require(None not in row and row['solver']==name and int(row['step'])==step,'History actual row/order differs')
                expected=dict(time_s=b['time_s'][step],inventory_kg_per_m=b['inventory'][step],scalar_min=b['c'][step].min(),scalar_max=b['c'][step].max(),negative_cell_count=int((b['c'][step]<0).sum()),above_one_cell_count=int((b['c'][step]>1).sum()))
                for key,value in expected.items():require(float(row[key])==value,'History native-field statistic differs: '+name+'/'+str(step)+'/'+key)
                state=b['rows'][step]
                if step==0:
                    require(state.get('budget') is None and state.get('residual') is None and all(row[key]=='' for key in HISTORY_HEADER[8:]),'Initial unavailable diagnostics were invented or changed')
                else:
                    budget=e.value(state['budget']);residual=e.value(state['residual']);require('RECONSTRUCTED' in budget['classification'],'Actual reconstructed budget identity missing')
                    require(float(row['top_inlet_c'])==state['prescribed_top_c'],'History actual inlet sample differs')
                    require(float(row['global_Euler_defect_relative'])==budget['global_Euler_defect_relative'] and float(row['own_terminal_scalar_residual'])==residual,'History actual budget/residual differs')
                    value=budget['reconstructed_outlet_outflow_weighted_mean_c'];require(row['outlet_outflow_weighted_c']=='' if value is None else float(row['outlet_outflow_weighted_c'])==value,'History outlet value/unavailability differs')
                    inventory_error=abs(float(row['inventory_kg_per_m'])-budget['inventory_kg_per_m']);max_inventory_budget_error=max(max_inventory_budget_error,inventory_error)
                    require(np.isclose(float(row['inventory_kg_per_m']),budget['inventory_kg_per_m'],rtol=1e-13,atol=1e-15),'Existing inventory/budget roundoff correspondence exceeded')
                    for patch in PATCHES:
                        for term in ('advective','diffusive','total'):
                            key=term+'_inward_kg_m_s';require(float(row[patch+'_'+key])==budget['reconstructed_boundary'][patch][key],'History actual boundary component differs')
                count+=1
        require(next(reader,None) is None,'Unexpected extra history CSV rows')
    require(count==1602,'Actual history CSV row count differs');return dict(rows=count,all_field_statistics_exact=True,all_available_budget_residual_outlet_boundary_values_exact=True,initial_unavailable_preserved=True,inventory_exact_independent_native_field_reduction=True,inventory_budget_existing_roundoff_rtol=1e-13,inventory_budget_existing_roundoff_atol=1e-15,maximum_inventory_budget_roundoff_difference_kg_per_m=max_inventory_budget_error)

def clock(raw,keys,scheme):
    values={}
    for key in keys:
        pattern=r'\('+re.escape(key)+r'\s+('+FLOAT+r')\)' if scheme else r'\b'+re.escape(key)+r'\s+('+FLOAT+r')\s*;'
        found=re.findall(pattern,raw);require(len(found)==1,'Unambiguous actual native clock required');values[key]=float(found[0])
    return values
def native_cell_scalar(h,path):
    root=h[path];datasets=[root] if isinstance(root,h5py.Dataset) else list(root.values());seen=np.zeros(N,bool);c=np.empty(N)
    for dataset in datasets:
        require(isinstance(dataset,h5py.Dataset),'Actual native scalar chunk expected');lo=int(np.asarray(dataset.attrs['minId']).ravel()[0])-1;hi=int(np.asarray(dataset.attrs['maxId']).ravel()[0]);values=np.asarray(dataset[()]).ravel()
        require(0<=lo<hi<=N and len(values)==hi-lo and not seen[lo:hi].any(),'Actual native scalar coverage differs');seen[lo:hi]=True;c[lo:hi]=values
    require(seen.all() and np.isfinite(c).all(),'Complete actual finite native scalar required');return c

def verify_native_frames(bundles,render_checks,e):
    for name,pin in HELPERS.items():e.pin(dict(path=str(ROOT/'repo/cfd-learning-lab/scripts'/name),sha256=pin))
    parser=e.module('independent_tracer_visual_OF_parser',Q,Q_SHA);require(parser.HELPERS==HELPERS,'Qualified scalar parser helper binding differs')
    for attr,name in [('old','prepare_openfoam_tjunction_v001.py'),('scalar_helper','prepare_openfoam_heated_channel_pilot_v003.py')]:require(Path(getattr(parser,attr).__file__).resolve()==(ROOT/'repo/cfd-learning-lab/scripts'/name).resolve(),'Loaded native scalar helper identity differs')
    result={}
    for name in NAMES:
        b=bundles[name];proof=render_checks['raw_selected_native_readback'][name];require(len(proof)==101,'Actual selected native proof count differs');verified=[]
        for frame,step in enumerate(STEPS):
            state=b['rows'][step];native=state['native'];declared=proof[frame]
            require(declared['step']==step and declared['time_s']==b['time_s'][step] and declared['raw_scalar_exact'] is True,'Selected native proof step/time differs')
            if name=='fluent':
                case=e.pin(native['case']);data=e.pin(native['data']);require(project_path(declared['case'])==case and project_path(declared['data'])==data and declared['scalar_dataset']==native['scalar_dataset'],'Selected Fluent proof paths/dataset differ')
                with h5py.File(data,'r') as h:current=native_cell_scalar(h,native['scalar_dataset']);actual=clock(h['settings/Data Variables'][0].decode(),('flow-time','time-step','physical-time-step'),True)
                require(actual==declared['native_clock'] and abs(actual['flow-time']-step*DT)<=1e-12 and actual['time-step']==step and abs(actual['physical-time-step']-DT)<=1e-12,'Actual selected Fluent native clock differs')
            else:
                field=e.pin(native['field']);require(project_path(declared['field'])==field,'Selected OF proof path differs');raw=field.read_text(encoding='ascii');require(re.search(r'\bobject\s+tracer\s*;',raw),'Actual OF tracer identity differs');current=np.asarray(parser.field(field,N,1,[0]*7),dtype=float)
                for ref in native.get('lineage',[]):e.pin(ref)
                if step:
                    path=e.pin(native['time_metadata']);actual=clock(path.read_text(encoding='ascii'),('value','index','deltaT'),False)
                    require(project_path(declared['time_metadata'])==path and actual==declared['native_clock'] and abs(actual['value']-step*DT)<=1e-12 and actual['index']==step and abs(actual['deltaT']-DT)<=1e-12,'Actual selected OF native clock differs')
                else:require(declared['initial_clock_verified_in_selected_audit'] is True and 'native_clock' not in declared,'Initial OF clock proof was invented')
            require(np.array_equal(current[b['source_native_cell_index']],b['c'][step]),'Actual selected native scalar differs from accepted canonical NPZ')
            verified.append(dict(step=step,time_s=float(b['time_s'][step]),native_scalar_exact=True))
        result[name]=dict(selected_states=verified,raw_native_cell_values_checked=101*N,mesh_identity_scope='Authenticated old geometry plus selected accepted full native audit; no duplicate all801 mesh audit')
    return result
def verify_gif(render,checks):
    a=checks['animation'];require(a['frames']==101 and a['actual_native_steps']==STEPS and a['actual_times_s']==[i*DT for i in STEPS] and a['playback_frame_delay_ms']==60,'Actual GIF state/time metadata differs')
    require(a['fixed_common_c_range']==[0,1] and a['input_arrays_clipped'] is False and a['temporal_interpolation'] is False and a['particles'] is False and a['equal_physical_xy_aspect'] is True and a['exact_sparse_source_quads'] is True,'Actual GIF physical/display scope differs')
    p=render/(PREFIX+'-actual-concentration.gif');require(sha(p)==a['sha256'] and p.stat().st_size==a['bytes'],'Actual GIF hash/size differs');full=[];maps={n:[] for n in NAMES};contacts={}
    for name,crop in a['map_crop_pixel_bounds'].items():require(name in NAMES and len(crop)==4 and 0<=crop[0]<crop[2]<=1280 and 0<=crop[1]<crop[3]<=640,'Actual map crop bounds invalid')
    with Image.open(p) as gif:
        require(gif.n_frames==101 and gif.size==(1280,640) and gif.info.get('loop')==0,'Decoded actual GIF frames/size/loop differs')
        for i in range(101):
            gif.seek(i);require(gif.info.get('duration')==60,'Decoded actual GIF delay differs');rgb=gif.convert('RGB');full.append(hashlib.sha256(rgb.tobytes()).hexdigest())
            for name in NAMES:maps[name].append(hashlib.sha256(rgb.crop(a['map_crop_pixel_bounds'][name]).tobytes()).hexdigest())
            if i in (0,50,100):contacts[i]=rgb.copy()
    require(full==a['decoded_frame_hashes'] and len(set(full))==a['distinct_decoded_RGB_frames'],'Every decoded actual frame hash differs from render proof')
    require({n:len(set(v)) for n,v in maps.items()}==a['map_region_unique_frames'],'Decoded map-region uniqueness proof differs')
    contact=render/(PREFIX+'-decoded-contact-sheet.png')
    with Image.open(contact) as image:
        require(image.size==(1280,1998),'Actual decoded contact sheet size differs');image=image.convert('RGB')
        for slot,i in enumerate((0,50,100)):
            crop=image.crop((0,slot*666+26,1280,slot*666+666));require(crop.tobytes()==contacts[i].tobytes(),'Actual contact-sheet frame crop differs from decoded GIF')
    return dict(decoded_frames=101,decoded_delay_ms=60,pixel_dimensions=[1280,640],all_decoded_frame_hashes_exact=True,decoded_contact_frame_indices=[0,50,100],decoded_contact_native_steps=[0,400,800],decoded_contact_frame_pixels_exact=True,map_region_unique_frames={n:len(set(v)) for n,v in maps.items()},repeated_frames_allowed=True,image_to_scalar_scope='Raw selected scalar/time inputs independently verified; recorded rendered frame hashes and contact pixels checked. Academic image framing/color interpretation still receives root visual inspection.')

def manifests(render,e):
    input_path=render/'input-manifest.json';output_path=render/'output-manifest.json';inputs=e.json(dict(path=str(input_path),sha256=sha(input_path)));outputs=e.json(dict(path=str(output_path),sha256=sha(output_path)))
    require(inputs['inputs_unchanged'] is True and outputs['manifest_self_excluded'] is True,'Actual render manifest completion flags missing')
    for row in inputs['files']:
        require(row['sha256_before']==row['sha256_after'],'Render input changed during render');p=e.pin(dict(path=row['path'],sha256=row['sha256_after']));require(p.stat().st_size==row['bytes'],'Render input size differs')
    expected={PREFIX+'-final-concentration.png',PREFIX+'-final-concentration.svg',PREFIX+'-actual-concentration.gif',PREFIX+'-decoded-contact-sheet.png',PREFIX+'-selected-concentration.png',PREFIX+'-selected-concentration.svg',PREFIX+'-transport-history.png',PREFIX+'-transport-history.svg',PREFIX+'-numerical-history.png',PREFIX+'-numerical-history.svg','actual-native-scalar-cell-fields.csv.gz','record-final-scalar-fields.csv','record-scalar-step-history.csv','input-manifest.json','checks.json'}
    declared=set()
    for row in outputs['files']:
        p=project_path(row['path']);require(p.parent==render and row['relative_path']==p.name and p.name not in declared,'Actual output manifest duplicate/path escape');declared.add(p.name);e.pin(dict(path=str(p),sha256=row['sha256']));require(p.stat().st_size==row['bytes'],'Actual rendered output size differs')
    require(declared==expected and {p.name for p in render.iterdir()}==expected|{'output-manifest.json'},'Complete actual render output coverage differs')
    source=PRIVATE/'plot_tracer_tjunction_native_v004.py';source_sha='60c1bb85858beedbd6f68c08818bc620a942339ed4640a76d061bf5cabc906d2'
    require(any(project_path(row['path'])==source and row['sha256_after']==source_sha for row in inputs['files']),'Actual render did not bind reviewed plotterv004 source')
    for p in render.glob('*.png'):
        with Image.open(p) as image:image.verify()
    return dict(render_input_count=len(inputs['files']),render_output_count=len(outputs['files']),complete_output_manifest_coverage=True,all_declared_output_hashes_and_sizes_match=True,PNG_decode_integrity=True)

def write_json(path,value):
    with path.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--selection',type=Path,required=True);parser.add_argument('--selection-sha256',required=True);parser.add_argument('--render-dir',type=Path,required=True);parser.add_argument('--source-sha256',required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    require(sha(__file__)==args.source_sha256,'Reviewed independent verifier source differs');available=memory_guard();render=project_path(args.render_dir);require(render.parent==ROOT/'exports/figures' and render.is_dir(),'Actual versioned E render directory required');out=project_path(args.out);require(out.parent==PRIVATE and re.fullmatch(r'visual-export-verification-v\d{3}',out.name) and not out.exists(),'Fresh private versioned verification output required')
    e=Evidence();e.pin(dict(path=str(Path(__file__).resolve()),sha256=args.source_sha256));selection=e.json(dict(path=str(args.selection),sha256=args.selection_sha256));require(selection['root_authorizes_offline_render'] is True and selection['actual_paired_physical_steps']==800 and selection['dt_s']==DT,'Actual reviewed root full render selection required')
    for ref in selection['context_inputs']:e.pin(ref)
    manifest_proof=manifests(render,e);checks=e.json(dict(path=str(render/'checks.json'),sha256=sha(render/'checks.json')));require(checks['input_hashes_stable'] is True and checks['input_arrays_modified'] is False and checks['initial_solved_budget_claim'] is False,'Actual render field/provenance flags differ')
    g,native=load_geometry(e);bundles={name:load_bundle(name,selection['solvers'][name],g,native,e) for name in NAMES}
    full=verify_cell_csv(render/'actual-native-scalar-cell-fields.csv.gz',bundles,g,True);final=verify_cell_csv(render/'record-final-scalar-fields.csv',bundles,g,False);history=verify_history(render/'record-scalar-step-history.csv',bundles,e)
    for name,b in bundles.items():
        expected=dict(minimum=b['c'].min(axis=1),maximum=b['c'].max(axis=1),negative=(b['c']<0).sum(axis=1),above_one=(b['c']>1).sum(axis=1))
        for key,value in expected.items():require(np.array_equal(np.asarray(checks['scalar_extrema_counts'][name][key]),value),'Render all-state raw concentration extrema/counts differ')
    raw=verify_native_frames(bundles,checks,e);animation=verify_gif(render,checks);inputs=e.finish();after_ram=memory_guard()
    result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),kind='Independent actual paired tracer CSV/GIF/native-frame verification',all_declared_readback_checks_passed=True,source_sha256=args.source_sha256,full_scalar_CSV_rows=full,final_scalar_CSV_rows=final,all_CSV_c_values_and_order_exact=True,CSV_geometry_scope='Canonical source polygon centres and areas as actually plotted; not a literal native-coordinate dump',canonical_CSV_geometry_exact=True,native_NPZ_metadata_exact_authenticated_stencil=True,native_canonical_geometry_errors={n:dict(coordinates_m=b['native_geometry_error_m'],areas_m2=b['native_geometry_area_error_m2']) for n,b in bundles.items()},history=history,selected_raw_native=raw,animation=animation,manifests=manifest_proof,all_declared_input_output_hashes_stable=True,initial_diagnostics_unavailable_preserved=True,render_or_native_inputs_modified=False,CFD_or_WSL_or_GUI_started=False,plotter_or_builder_main_called=False,actual_image_regeneration_performed=False,root_visual_layout_review_pending=True,physical_validation=False,mesh_independence=False,time_step_independence=False,complete_washout_established=False,available_ram_before_bytes=available,available_ram_after_bytes=after_ram)
    out.mkdir();write_json(out/'verification.json',result);write_json(out/'input-manifest.json',dict(all_inputs_unchanged=True,files=inputs));outputs=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.iterdir())];write_json(out/'output-manifest.json',dict(manifest_self_excluded=True,files=outputs));print(json.dumps(dict(out=str(out),verification_sha256=sha(out/'verification.json'),full_rows=full,history_rows=history['rows'],final_rows=final,decoded_frames=animation['decoded_frames'])),flush=True)
if __name__=='__main__':main()
