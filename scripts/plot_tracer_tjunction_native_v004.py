"""Offline paired actual tracer exports. No CFD, WSL, GUI or publication.
Root supplies accepted files and explicit budget/residual JSON key paths.
"""
from pathlib import Path
import argparse,csv,ctypes,datetime,gzip,hashlib,importlib.util,json,os,re,sys
os.environ['TEMP']=os.environ['TMP']='E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-tracer-v001'
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image,ImageDraw
ROOT=Path('E:/CFD-Learning-Lab'); PRIVATE=ROOT/'local/prospective-tracer-v001'
GEOMETRY=ROOT/'native/runs/2026-11-01-tjunction-v001/mesh/geometry.npz'
GEOMETRY_SHA='671236b7bf7658fcbae5330b4cd6ce66c8e2dee9ff06753ebfc45c0a4cfb2949'
ADAPTER=PRIVATE/'scalar_native_geometry_v003.py'
ADAPTER_SHA='06f6b4089116d8e7cdab8b4050c30d2939d3876381ebbaa878c7908e6f2d1903'
FV=PRIVATE/'scalar_fv_audit_v002.py';FV_SHA='75ceec60add2f5c53fb82585d105dcc5cd7c9d8d19e1ea8ac3f2dde14aa455cd'
OF_SCALAR_PARSER=PRIVATE/'prepare_openfoam_tracer_qualification_v004.py'
OF_SCALAR_PARSER_SHA='278e6948a244e82887d8eb8dc2da90b502a14a7165de2238c3ba1e01b2890888'
OF_PARSER_HELPERS={'prepare_openfoam_tjunction_v001.py':'16ded13dcc9d8e7e315f5aaefe24d84af2963e910efda6b48c0d355fe04965d4','prepare_openfoam_heated_channel_pilot_v003.py':'807c30de1c7c7400ab38d2b4592896db774be4c98ec2b853afd2f9fff2ed90b0'}
N,STATES,DT=4400,801,.025
FRAMES=np.arange(0,801,8,dtype=int)
NAMES={'fluent':'ANSYS Fluent','openfoam':'OpenFOAM'}
COLORS={'fluent':'#008498','openfoam':'#d65d28'}
PATCHES=('inlet_left','inlet_top','outlet','walls')
FLOAT=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?';PREFIX='tracer-tjunction'
def require(v,m):
    if not v:raise ValueError(m)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def within(path,root):
    require(Path(path).is_absolute(),'Absolute project path required')
    path,root=Path(path).resolve(),Path(root).resolve()
    require(path==root or root in path.parents,'Path outside project: '+str(path));return path
class MemoryStatus(ctypes.Structure):
    _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(k,ctypes.c_ulonglong) for k in ('total','available','total_page','available_page','total_virtual','available_virtual','available_extended')]
def resources(minimum_gib):
    import shutil
    m=MemoryStatus();m.length=ctypes.sizeof(m)
    require(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)),'RAM query failed')
    require(m.available>=minimum_gib*1024**3,'Fresh render RAM floor not met')
    free=shutil.disk_usage(ROOT).free;require(free>=10*1024**3,'E: free-space floor not met')
    return dict(available_ram_bytes=m.available,E_free_bytes=free,floor_ram_gib=minimum_gib)
class Inputs:
    def __init__(self):self.hashes={};self.cache={}
    def pin(self,ref):
        require(isinstance(ref,dict) and set(ref)>={'path','sha256'},'Explicit file/SHA reference required')
        p=within(ref['path'],ROOT);require(p.is_file(),'Missing input: '+str(p));h=ref['sha256'].lower()
        require(re.fullmatch('[0-9a-f]{64}',h) is not None,'Valid SHA256 required')
        if p not in self.hashes:require(sha(p)==h,'Input SHA mismatch: '+str(p));self.hashes[p]=h
        require(self.hashes[p]==h,'Conflicting input hashes');return p
    def json(self,ref):
        p=self.pin(ref)
        if p not in self.cache:self.cache[p]=json.loads(p.read_text(encoding='utf-8-sig'))
        return self.cache[p]
    def value(self,ref):
        require(isinstance(ref.get('key_path'),list),'Explicit JSON key_path required');v=self.json(ref)
        for k in ref['key_path']:
            require(isinstance(k,(str,int)) and not isinstance(k,bool),'Literal JSON key/index required');v=v[k]
        return v
    def finish(self):
        rows=[]
        for p,h in self.hashes.items():
            actual=sha(p);require(actual==h,'Input changed: '+str(p));rows.append(dict(path=str(p),bytes=p.stat().st_size,sha256_before=h,sha256_after=actual))
        return rows
def load_geometry(inputs):
    inputs.pin(dict(path=str(FV),sha256=FV_SHA));p=inputs.pin(dict(path=str(ADAPTER),sha256=ADAPTER_SHA))
    s=importlib.util.spec_from_file_location('tracer_visual_qualified_geometry',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    native,proof=m.qualified_geometry()
    for row in proof['consumed_inputs']:inputs.pin(row)
    p=inputs.pin(dict(path=str(GEOMETRY),sha256=GEOMETRY_SHA))
    with np.load(p,allow_pickle=False) as z:g={k:z[k].copy() for k in ('nodes','quads','centres','areas')}
    require(g['nodes'].shape==(4641,2) and g['quads'].shape==(N,4),'Canonical T geometry shape')
    g['polygons']=g['nodes'][g['quads']];require(abs(float(g['areas'].sum())-.44)<1e-12,'Canonical domain area')
    return g,native,proof
def validate_arrays(step,time,c):
    require(step.shape==(STATES,) and np.issubdtype(step.dtype,np.integer) and np.array_equal(step,np.arange(STATES)),'Complete integer0..800 steps required')
    require(time.shape==(STATES,) and np.isfinite(time).all() and np.max(abs(time-step*DT))<=1e-12,'Actual0..20s saved time history required')
    require(c.shape==(STATES,N) and c.dtype==np.float64 and np.isfinite(c).all(),'Finite double source-order scalar history required')
    return dict(minimum=c.min(axis=1),maximum=c.max(axis=1),negative=(c<0).sum(axis=1),above_one=(c>1).sum(axis=1))
def read_budget(inputs,row):
    b=inputs.value(row['budget']);require('RECONSTRUCTED' in b['classification'],'Explicit reconstructed budget classification')
    require(set(b['reconstructed_boundary'])==set(PATCHES),'Four actual budget patches required');return b
def verify_export_geometry(centres,areas,g,native,inverse):
    # Exports retain actual native metadata; canonical polygons have qualified roundoff differences.
    stencil=native['stencil'];expected_centres=stencil['cell_centres'][inverse]
    expected_areas=(stencil['volumes']/stencil['depth_m'])[inverse]
    require(np.array_equal(centres,expected_centres) and np.array_equal(areas,expected_areas),'Export metadata differs from authenticated native stencil')
    coordinate_error=float(np.max(abs(centres-g['centres'])));area_error=float(np.max(abs(areas-g['areas'])))
    require(coordinate_error<=1e-12 and area_error<=1e-14,'Existing qualified native/canonical geometry tolerances exceeded')
    return dict(export_metadata_exact_native=True,maximum_native_canonical_coordinate_error_m=coordinate_error,maximum_native_canonical_area_error_m2=area_error,existing_coordinate_tolerance_m=1e-12,existing_area_tolerance_m2=1e-14,metadata_modified=False)
def load_of_scalar_parser(inputs):
    # Authenticate the qualifier and all of its nonstdlib imported helpers before import.
    paths={name:inputs.pin(dict(path=str(ROOT/'repo/cfd-learning-lab/scripts'/name),sha256=pin)) for name,pin in OF_PARSER_HELPERS.items()}
    p=inputs.pin(dict(path=str(OF_SCALAR_PARSER),sha256=OF_SCALAR_PARSER_SHA));name='tracer_visual_qualified_scalar_parser_v004'
    if name in sys.modules:
        module=sys.modules[name];require(Path(module.__file__).resolve()==p,'Pinned scalar parser namespace collision')
    else:
        spec=importlib.util.spec_from_file_location(name,p);module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    require(module.HELPERS==OF_PARSER_HELPERS,'Scalar parser helper declarations differ')
    for attribute,filename in (('old','prepare_openfoam_tjunction_v001.py'),('scalar_helper','prepare_openfoam_heated_channel_pilot_v003.py')):
        require(Path(getattr(module,attribute).__file__).resolve()==paths[filename],'Loaded scalar parser helper identity differs')
    return module
def load_history(name,selection,g,native,inputs):
    require(selection['full_actual_audit_accepted_by_root'] is True,'Root full-audit acceptance required')
    require(selection['accepted_evidence'],'Literal accepted audit evidence required')
    for ref in selection['accepted_evidence']:inputs.pin(ref)
    p=inputs.pin(selection['fields'])
    with np.load(p,allow_pickle=False) as z:
        step,time,c=(z[k].copy() for k in ('step','time_s','c'));mapping=z['native_cell_to_source_cell'].copy();inverse=z['source_native_cell_index'].copy()
        centres,areas=(z[k].copy() for k in ('source_centres_m','source_areas_m2'))
    stats=validate_arrays(step,time,c);actual=np.asarray(native[name]['native_cell_to_source_cell'])
    require(np.array_equal(mapping,actual) and np.array_equal(np.sort(mapping),np.arange(N)),'Native/source bijection differs')
    require(np.array_equal(inverse,np.argsort(mapping)),'Inverse source/native mapping differs')
    geometry_metadata=verify_export_geometry(centres,areas,g,native[name],inverse)
    manifest=inputs.json(selection['manifest']);rows=manifest['states'];require(len(rows)==STATES,'Literal paths for all801 states required')
    budgets,residuals,schedule=[],[],[]
    for index,row in enumerate(rows):
        require(row['step']==index and abs(row['time_s']-time[index])<=1e-12,'Manifest step/time differs')
        require(row['native_field_name']==('SV_UDS_0' if name=='fluent' else 'tracer'),'Observed native scalar identity differs')
        if name=='fluent':inputs.pin(row['native']['case']);inputs.pin(row['native']['data'])
        else:
            inputs.pin(row['native']['field'])
            for ref in row['native'].get('lineage',[]):inputs.pin(ref)
            if index:inputs.pin(row['native']['time_metadata'])
        if index:
            require(row['prescribed_top_c']==(1. if index<=80 else 0.) and row['prescribed_left_c']==0.,'Actual sampled inlet schedule differs')
            schedule.append(float(row['prescribed_top_c']));budgets.append(read_budget(inputs,row));r=float(inputs.value(row['residual']))
            require(np.isfinite(r) and r>=0,'Finite own terminal residual required');residuals.append(r)
        else:require(row.get('budget') is None and row.get('residual') is None,'Initial has no solved budget/residual')
    rho=native[name]['rho'];stencil=native[name]['stencil'];native_c=c[:,mapping]
    inventory=np.sum(native_c*(rho*stencil['volumes'])[None,:],axis=1)/stencil['depth_m']
    observed=np.array([b['inventory_kg_per_m'] for b in budgets])
    require(np.allclose(observed,inventory[1:],rtol=1e-13,atol=1e-15),'Budget inventory differs from actual c')
    boundary=np.array([[[b['reconstructed_boundary'][patch][term+'_inward_kg_m_s'] for term in ('advective','diffusive','total')] for patch in PATCHES] for b in budgets])
    outlet=np.array([np.nan if b['reconstructed_outlet_outflow_weighted_mean_c'] is None else b['reconstructed_outlet_outflow_weighted_mean_c'] for b in budgets])
    defect=np.array([b['global_Euler_defect_relative'] for b in budgets]);require(np.isfinite(boundary).all() and np.isfinite(defect).all(),'Finite reconstructed histories')
    return dict(name=name,c=c,step=step,time=time,stats=stats,mapping=mapping,inverse=inverse,geometry_metadata=geometry_metadata,rows=rows,budgets=budgets,inventory=inventory,boundary=boundary,outlet=outlet,defect=defect,residual=np.array(residuals),residual_label=selection['residual_label'],schedule=np.array(schedule),accepted_evidence=selection['accepted_evidence'])
def hdf_scalar(data,dataset):
    root=data[dataset];chunks=[root] if isinstance(root,h5py.Dataset) else list(root.values());out=np.empty(N);seen=np.zeros(N,bool)
    for d in chunks:
        require(isinstance(d,h5py.Dataset),'Explicit scalar chunks required');lo=int(np.asarray(d.attrs['minId']).ravel()[0])-1;hi=int(np.asarray(d.attrs['maxId']).ravel()[0]);v=np.asarray(d[()]).ravel()
        require(0<=lo<hi<=N and len(v)==hi-lo and not seen[lo:hi].any(),'Scalar coverage differs');out[lo:hi]=v;seen[lo:hi]=True
    require(seen.all() and np.isfinite(out).all(),'Complete finite native scalar required');return out
def clock_values(raw,keys,syntax):
    out={}
    for key in keys:
        pattern=r'\('+re.escape(key)+r'\s+('+FLOAT+r')\)' if syntax=='scheme' else r'\b'+re.escape(key)+r'\s+('+FLOAT+r')\s*;'
        found=re.findall(pattern,raw);require(len(found)==1,'Unambiguous actual clock: '+key);out[key]=float(found[0])
    return out
def verify_native_frames(b,inputs):
    proof=[];parser=load_of_scalar_parser(inputs) if b['name']=='openfoam' else None
    for index in FRAMES:
        row=b['rows'][int(index)];native=row['native']
        if b['name']=='fluent':
            case=inputs.pin(native['case']);p=inputs.pin(native['data'])
            with h5py.File(p,'r') as h:current=hdf_scalar(h,native['scalar_dataset']);raw=h['settings/Data Variables'][0].decode()
            t=clock_values(raw,('flow-time','time-step','physical-time-step'),'scheme')
            require(abs(t['flow-time']-index*DT)<=1e-12 and t['time-step']==index and abs(t['physical-time-step']-DT)<=1e-12,'Actual Fluent selected clock differs')
            detail=dict(case=str(case),data=str(p),scalar_dataset=native['scalar_dataset'],native_clock=t)
        else:
            p=inputs.pin(native['field']);raw=p.read_text(encoding='utf-8-sig');require(re.search(r'\bobject\s+tracer\s*;',raw),'Actual OF tracer object required')
            current=np.asarray(parser.field(p,N,1,[0,0,0,0,0,0,0]),dtype=float);detail=dict(field=str(p),initial_clock_verified_in_selected_audit=bool(index==0))
            if index:
                tp=inputs.pin(native['time_metadata']);t=clock_values(tp.read_text(encoding='utf-8-sig'),('value','index','deltaT'),'foam')
                require(abs(t['value']-index*DT)<=1e-12 and t['index']==index and abs(t['deltaT']-DT)<=1e-12,'Actual OF selected clock differs');detail.update(time_metadata=str(tp),native_clock=t)
        require(np.array_equal(current[b['inverse']],b['c'][index]),'Raw selected c differs from accepted canonical export')
        proof.append(dict(step=int(index),time_s=float(b['time'][index]),raw_scalar_exact=True,**detail))
    return proof
def walls():return [[[-.4,-.1],[1.2,-.1]],[[-.4,.1],[-.1,.1]],[[-.1,.1],[-.1,.7]],[[.1,.1],[.1,.7]],[[.1,.1],[1.2,.1]]]
def cmap():
    c=matplotlib.colormaps['viridis'].copy();c.set_under('#e34878');c.set_over('#f07b27');return c
def decorate(ax,title):
    ax.add_collection(LineCollection(walls(),colors='#294653',linewidths=1.1,zorder=5))
    ax.plot([-.4,-.4],[-.1,.1],color='#008498',lw=1.5,ls='--');ax.plot([-.1,.1],[.7,.7],color='#d65d28',lw=1.5,ls='--');ax.plot([1.2,1.2],[-.1,.1],color='#6d8793',lw=1,ls=':')
    ax.set(xlim=(-.4,1.2),ylim=(-.1,.7),xlabel='x (m)',ylabel='y (m)');ax.set_aspect('equal');ax.set_title(title,loc='left',fontsize=10,pad=9,fontweight='bold');ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8);ax.set_facecolor('white')
def text_bounds(fig):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();bad=[]
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text().strip():
            box=t.get_window_extent(renderer)
            if box.x0< -1 or box.y0< -1 or box.x1>fig.bbox.width+1 or box.y1>fig.bbox.height+1:bad.append(t.get_text())
    require(not bad,'Text outside canvas: '+repr(bad))
def canvas(g,bundles):
    fig=plt.figure(figsize=(12.8,6.4),dpi=100);fig.text(.07,.95,'A dye pulse through a T-junction',fontsize=21,fontweight='bold')
    fig.text(.07,.895,'Frozen laminar flow | Re = 20 at each inlet | D = 0.001 m2/s | dt = 0.025 s',fontsize=10,color='#425e6d');clock=fig.text(.07,.838,'',fontsize=9.5,color='#304b59');artists={};axes={}
    for left,name in ((.07,'fluent'),(.56,'openfoam')):
        ax=fig.add_axes([left,.30,.40,.40]);col=PolyCollection(g['polygons'],array=bundles[name]['c'][0],cmap=cmap(),norm=Normalize(0,1,clip=False),edgecolors='none',antialiased=False);ax.add_collection(col);decorate(ax,NAMES[name]+' | actual native concentration');artists[name]=col;axes[name]=ax
    cax=fig.add_axes([.285,.18,.43,.022]);fig.colorbar(ScalarMappable(norm=Normalize(0,1,clip=False),cmap=cmap()),cax=cax,orientation='horizontal',extend='both',label='Dimensionless concentration c | fixed common scale')
    fig.text(.07,.085,'Actual scalar transport on saved native fields; no animated particles or interpolation.',fontsize=8.5,color='#425e6d')
    fig.text(.07,.048,'2D | first-order Euler/upwind | 20 s early washout | mesh/time sensitivity and physical validation not established',fontsize=8,color='#425e6d');return fig,artists,axes,clock
def palette():
    colors=np.rint(matplotlib.colormaps['viridis'](np.linspace(0,1,128))[:,:3]*255).astype('uint8').tolist();colors += [[v,v,v] for v in range(0,256,4)]
    colors += [[255,255,255],[41,70,83],[66,94,109],[48,75,89],[0,132,152],[214,93,40],[227,72,120],[240,123,39]];colors += [[255,255,255]]*(256-len(colors));p=Image.new('P',(1,1));p.putpalette(np.array(colors,dtype='uint8').ravel().tolist());return p
def render_animation(g,bundles,out):
    fig,artists,axes,clock=canvas(g,bundles);frames=[];crops={};pal=palette()
    for j,step in enumerate(FRAMES):
        if j%10==0:resources(2.5)
        for name,b in bundles.items():artists[name].set_array(b['c'][step])
        phase='saved initial c = 0 field' if step==0 else 'top inlet c = '+str(int(bundles['fluent']['schedule'][step-1]))
        counts=[f'{bundles[n]["stats"]["negative"][step]}/{bundles[n]["stats"]["above_one"][step]}' for n in bundles]
        clock.set_text(f'Actual physical time {step*DT:.2f} s | native step {step} | {phase} | below 0 / above 1 cells F, OF: '+', '.join(counts))
        if j in (0,len(FRAMES)-1):text_bounds(fig)
        fig.canvas.draw();im=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy());frames.append(im.quantize(palette=pal,dither=Image.Dither.NONE))
        for name,ax in axes.items():
            box=ax.get_window_extent();crops[name]=[round(box.x0)+2,round(im.height-box.y1)+2,round(box.x1)-2,round(im.height-box.y0)-2]
    fig.savefig(out/(PREFIX+'-final-concentration.png'),dpi=170);fig.savefig(out/(PREFIX+'-final-concentration.svg'));p=out/(PREFIX+'-actual-concentration.gif');frames[0].save(p,save_all=True,append_images=frames[1:],duration=60,loop=0,optimize=False,disposal=2);plt.close(fig)
    require(p.stat().st_size<=15*1024**2,'GIF exceeds asset size target; preserve attempt');hashes=[];regions={n:[] for n in bundles};delays=[]
    with Image.open(p) as gif:
        require(gif.n_frames==101 and gif.size==(1280,640),'Decoded GIF count/dimensions differ');sheet=Image.new('RGB',(1280,666*3),'#e9eff2');draw=ImageDraw.Draw(sheet)
        for j in range(gif.n_frames):
            gif.seek(j);im=gif.convert('RGB');delays.append(gif.info.get('duration'));hashes.append(hashlib.sha256(im.tobytes()).hexdigest())
            for name,crop in crops.items():regions[name].append(hashlib.sha256(im.crop(crop).tobytes()).hexdigest())
            if j in (0,50,100):
                top=(0,50,100).index(j)*666;draw.text((12,top+6),f'Decoded actual step {FRAMES[j]} | physical time {FRAMES[j]*DT:.2f} s',fill='#203c4b');sheet.paste(im,(0,top+26))
        sheet.save(out/(PREFIX+'-decoded-contact-sheet.png'))
    require(delays==[60]*101,'Decoded GIF delays differ')
    return dict(frames=101,actual_native_steps=FRAMES.tolist(),actual_times_s=(FRAMES*DT).tolist(),playback_frame_delay_ms=60,actual_solver_scalar_evolution=True,temporal_interpolation=False,particles=False,fixed_common_c_range=[0,1],display_under_over_colors=True,input_arrays_clipped=False,equal_physical_xy_aspect=True,exact_sparse_source_quads=True,pixel_dimensions=[1280,640],decoded_frame_hashes=hashes,distinct_decoded_RGB_frames=len(set(hashes)),map_region_unique_frames={n:len(set(v)) for n,v in regions.items()},repeated_map_frames_retained=True,map_crop_pixel_bounds=crops,sha256=sha(p),bytes=p.stat().st_size)
def save_figure(fig,out,stem):
    text_bounds(fig);fig.savefig(out/(stem+'.png'),dpi=170);fig.savefig(out/(stem+'.svg'));plt.close(fig)
def render_snapshots(g,bundles,out):
    fig,axes=plt.subplots(4,2,figsize=(12.8,12.6));fig.subplots_adjust(left=.065,right=.95,top=.90,bottom=.11,hspace=.50,wspace=.22)
    fig.text(.065,.96,'Injection, transit, arrival and early washout',fontsize=19,fontweight='bold');fig.text(.065,.928,'Actual native concentration | common 0..1 scale | frozen laminar velocity | dt = 0.025 s',fontsize=9.5,color='#425e6d')
    for row,step in enumerate((80,240,480,800)):
        for col,(name,b) in enumerate(bundles.items()):
            ax=axes[row,col];ax.add_collection(PolyCollection(g['polygons'],array=b['c'][step],cmap=cmap(),norm=Normalize(0,1,clip=False),edgecolors='none',antialiased=False));decorate(ax,f'{NAMES[name]} | actual {step*DT:.2f} s | step {step}')
    cax=fig.add_axes([.285,.055,.43,.012]);fig.colorbar(ScalarMappable(norm=Normalize(0,1,clip=False),cmap=cmap()),cax=cax,orientation='horizontal',extend='both',label='Dimensionless c | values retained without clipping')
    fig.text(.065,.018,'2D, first-order Euler/upwind; 20 s is a finite window. Complete washout and physical validation are not established.',fontsize=8,color='#425e6d');save_figure(fig,out,PREFIX+'-selected-concentration')
def plot_style(axes):
    for ax in np.asarray(axes).ravel():ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.18);ax.tick_params(labelsize=8);ax.legend(fontsize=8,frameon=False)
def render_histories(bundles,out):
    fig,axes=plt.subplots(2,2,figsize=(12.8,8.6));fig.subplots_adjust(left=.085,right=.96,top=.85,bottom=.13,hspace=.38,wspace=.28)
    fig.text(.085,.95,'Pulse transport and conservation',fontsize=20,fontweight='bold');fig.text(.085,.905,'Saved scalar fields and RECONSTRUCTED Euler/upwind budgets; positive boundary transport is inward.',fontsize=9,color='#425e6d')
    axes[0,0].plot(bundles['fluent']['time'][1:],bundles['fluent']['schedule'],ls='none',marker='o',ms=2,color='#354e5a',label='Top inlet actual samples');axes[0,0].set(xlabel='Actual solved time (s)',ylabel='Prescribed c',title='The injected pulse')
    for name,b in bundles.items():
        color=COLORS[name];label=NAMES[name];t=b['time'];axes[0,1].plot(t,b['inventory'],color=color,label=label);axes[1,0].plot(t[1:],b['outlet'],color=color,label=label)
        integrated=np.r_[0.,np.cumsum(b['boundary'][:,:,2].sum(axis=1)*np.diff(t))]
        axes[1,1].plot(t,b['inventory']-b['inventory'][0],color=color,label=label+' inventory change');axes[1,1].plot(t,integrated,color=color,lw=1.1,ls='--',label=label+' integrated boundary input')
    axes[0,1].set(xlabel='Actual physical time (s)',ylabel='Inventory per unit depth (kg/m)',title='Stored scalar under rho = 1 kg/m3');axes[1,0].set(xlabel='Actual solved time (s)',ylabel='Outlet outflow-weighted c',title='Arrival and early washout');axes[1,1].set(xlabel='Actual physical time (s)',ylabel='Scalar inventory change (kg/m)',title='Euler accumulation includes inlet diffusion');plot_style(axes)
    fig.text(.085,.045,'Inventory is mass-weighted dimensionless concentration, not calibrated dye mass. Dashed curves integrate actual solved boundary rates times dt.',fontsize=8,color='#425e6d');save_figure(fig,out,PREFIX+'-transport-history')
    fig,axes=plt.subplots(2,2,figsize=(12.8,8.6));fig.subplots_adjust(left=.085,right=.96,top=.85,bottom=.15,hspace=.38,wspace=.28)
    fig.text(.085,.95,'Scalar numerical histories',fontsize=20,fontweight='bold');fig.text(.085,.905,'Native scalar extrema, reconstructed global defects and each solver\'s own terminal scalar residual.',fontsize=9,color='#425e6d')
    for name,b in bundles.items():
        color=COLORS[name];t=b['time'];label=NAMES[name];axes[0,0].plot(t,b['stats']['maximum'],color=color,label=label+' maximum');axes[0,0].plot(t,b['stats']['minimum'],color=color,ls='--',label=label+' minimum');axes[0,1].semilogy(t[1:],np.where(b['defect']>0,b['defect'],np.nan),color=color,label=label)
        ax=axes[1,0] if name=='fluent' else axes[1,1];ax.semilogy(t[1:],np.where(b['residual']>0,b['residual'],np.nan),color=color,label=label);ax.set(xlabel='Actual solved time (s)',ylabel=b['residual_label'],title=label+' | own terminal scalar norm')
    axes[0,0].set(xlabel='Actual physical time (s)',ylabel='Unmodified dimensionless c',title='All saved cell extrema; no clipping');axes[0,1].set(xlabel='Actual solved time (s)',ylabel='Relative reconstructed Euler defect',title='Storage minus signed boundary transfer');plot_style(axes)
    fig.text(.085,.07,'Zero values are omitted only from logarithmic drawing and retained in CSVs. Initial state has no solved residual or step budget.',fontsize=8,color='#425e6d');fig.text(.085,.035,'Residual definitions differ. Numerical acceptance comes from selected full audits; visual agreement does not establish physical validation.',fontsize=8,color='#425e6d');save_figure(fig,out,PREFIX+'-numerical-history')
def export_csv(g,bundles,out):
    fields=out/'actual-native-scalar-cell-fields.csv.gz';header=['solver','step','time_s','source_cell_index_zero_based','native_cell_index_zero_based','x_m','y_m','area_m2','c'];rows=0
    with gzip.open(fields,'xt',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(header)
        for name,b in bundles.items():
            for step in range(STATES):
                for cell in range(N):w.writerow([name,step,b['time'][step],cell,int(b['inverse'][cell]),*g['centres'][cell],g['areas'][cell],b['c'][step,cell]])
                rows+=N
    final=out/'record-final-scalar-fields.csv'
    with final.open('x',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(header)
        for name,b in bundles.items():
            for cell in range(N):w.writerow([name,800,b['time'][-1],cell,int(b['inverse'][cell]),*g['centres'][cell],g['areas'][cell],b['c'][-1,cell]])
    history=out/'record-scalar-step-history.csv';columns=['solver','step','time_s','inventory_kg_per_m','scalar_min','scalar_max','negative_cell_count','above_one_cell_count','top_inlet_c','outlet_outflow_weighted_c','global_Euler_defect_relative','own_terminal_scalar_residual']+[patch+'_'+term+'_inward_kg_m_s' for patch in PATCHES for term in ('advective','diffusive','total')]
    with history.open('x',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(columns)
        for name,b in bundles.items():
            for step in range(STATES):
                initial=step==0;j=step-1;row=[name,step,b['time'][step],b['inventory'][step],b['stats']['minimum'][step],b['stats']['maximum'][step],int(b['stats']['negative'][step]),int(b['stats']['above_one'][step])]
                row+=['','','',''] if initial else [b['schedule'][j],'' if np.isnan(b['outlet'][j]) else b['outlet'][j],b['defect'][j],b['residual'][j]];row+=['']*12 if initial else b['boundary'][j].ravel().tolist();w.writerow(row)
    require(rows==7048800,'Actual scalar CSV row total differs')
    return dict(all_native_cell_rows=rows,final_cell_rows=2*N,summary_rows=2*STATES,full_field_file=fields.name,record_final_file=final.name,record_history_file=history.name,csv_exact_double_roundtrip_readback_pending=True)
def write_json(path,value):
    with path.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--selection',type=Path,required=True);p.add_argument('--selection-sha256',required=True);p.add_argument('--source-sha256',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    require(sha(__file__)==a.source_sha256,'Reviewed plotter SHA differs');before=resources(3.);out=within(a.out,ROOT/'exports/figures');require(not out.exists(),'Fresh versioned figure output required')
    inputs=Inputs();inputs.pin(dict(path=str(Path(__file__).resolve()),sha256=a.source_sha256));s=inputs.json(dict(path=str(a.selection),sha256=a.selection_sha256));require(s['schema_version']==1 and s['root_authorizes_offline_render'] is True,'Reviewed root visual selection required');require(s['actual_paired_physical_steps']==800 and s['dt_s']==DT,'Actual paired20s/800step selection required')
    for ref in s['context_inputs']:inputs.pin(ref)
    g,native,gproof=load_geometry(inputs);bundles={name:load_history(name,s['solvers'][name],g,native,inputs) for name in NAMES};require(np.array_equal(bundles['fluent']['schedule'],bundles['openfoam']['schedule']),'Paired inlet samples differ');raw={name:verify_native_frames(b,inputs) for name,b in bundles.items()}
    out.mkdir(parents=True);animation=render_animation(g,bundles,out);render_snapshots(g,bundles,out);render_histories(bundles,out);csv_info=export_csv(g,bundles,out);input_rows=inputs.finish();write_json(out/'input-manifest.json',dict(inputs_unchanged=True,files=input_rows))
    checks=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),kind='Actual accepted paired passive tracer offline visuals',selected_full_audits_not_repeated=True,numerical_acceptance_source='Explicit root-selected accepted full audits',raw_selected_native_readback=raw,geometry=gproof,animation=animation,csv=csv_info,scalar_extrema_counts={name:{k:v.tolist() for k,v in b['stats'].items()} for name,b in bundles.items()},initial_solved_budget_claim=False,input_arrays_modified=False,scalar_field_interpolation=False,physical_validation=False,mesh_independence=False,time_step_independence=False,complete_washout_established=False,resources=dict(before=before,after=resources(2.5)),input_hashes_stable=True,independent_visual_and_CSV_QA_pending=True);write_json(out/'checks.json',checks)
    files=[dict(path=str(p),relative_path=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.iterdir()) if p.is_file()];write_json(out/'output-manifest.json',dict(manifest_self_excluded=True,files=files));print(json.dumps(dict(output=str(out),actual_animation=str(out/(PREFIX+'-actual-concentration.gif')),frames=101,cell_csv_rows=csv_info['all_native_cell_rows'],independent_QA_pending=True)),flush=True)
if __name__=='__main__':main()
