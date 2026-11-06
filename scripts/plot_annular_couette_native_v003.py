"""Annular Couette native fields: actual states, polygonal hole, no invented motion.

Geometry-only preparation is available without solver results. Pilot/full export
requires separately accepted native audits and explicitly selected native files.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,importlib.util,json,os,re,sys,textwrap
os.environ['TEMP']=os.environ['TMP']='E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-annular'
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image,ImageDraw
import plot_cavity_native_v003 as common
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-annular'
require,sha=common.require,common.sha
NAMES={'fluent':'ANSYS Fluent','openfoam':'OpenFOAM'}
COLOURS={'fluent':'#008498','openfoam':'#d65d28'}
GEOMETRY_SHA='d9fed1a699d1558e32de2fce7cb81776bf7294da8e4f405274d930c2625b1a97'
TARGETS_SHA='4c3f6bdb49ffb208230a7c89f574d9c030fd0d97de800e256f791ecf64abfff9'
FROZEN_HELPERS={'cavity':'c2425654f684bc8de8d9219dffa872d8b3c80a605bb35fce2b9070362f388445','native':'4596d6a6b975abc6110d80fdb8b5b680878cdff18022ce969866b288aacb2f02'}


def geometry(path):
    require(sha(path)==GEOMETRY_SHA,'Selected source geometry changed');g=dict(np.load(path,allow_pickle=False))
    require(g['nodes'].shape==(1600,2) and g['quads'].shape==(1536,4) and np.array_equal(g['grid_shape'],[24,64]),'Actual annular mesh size differs')
    poly=g['nodes'][g['quads']];x,y=poly[:,:,0],poly[:,:,1];cross=x*np.roll(y,-1,axis=1)-np.roll(x,-1,axis=1)*y;area=cross.sum(axis=1)/2
    centres=np.column_stack([((poly[:,:,i]+np.roll(poly[:,:,i],-1,axis=1))*cross).sum(axis=1) for i in (0,1)])/(6*area[:,None])
    require(np.all(area>0) and np.max(abs(area-g['areas']))<1e-14 and np.max(abs(centres-g['centres']))<1e-12,'Source actual polygon areas/centres differ')
    require(len(np.unique(g['nodes'],axis=0))==1600 and np.array_equal(g['radial_index'],np.repeat(np.arange(24),64)),'Annulus contains duplicate seam or unexpected ordering')
    edges={tuple(sorted(e)) for q in g['quads'] for e in zip(q,np.roll(q,-1))}
    require(1600-len(edges)+1536==0 and len(edges)==3136,'Annular Euler characteristic/hole differs')
    g['polygons']=poly;g['actual_polygon_area_m2']=float(area.sum());g['relative_area_difference_from_circular_annulus']=float(area.sum()/(np.pi*(.1**2-.05**2))-1)
    return g


class AcceptedSources:
    """Consume selected native files from an accepted unchanged audit manifest.

The full native audit owns all-state physics/face-flux checks. This export only
rehashes the raw frames/topology it independently reads, plus selected metadata.
"""
    def __init__(self,name,audit,manifest,fields,expected_sha):
        self.audit_path=audit.resolve();self.manifest_path=manifest.resolve();self.record=json.loads(audit.read_text(encoding='utf-8-sig'));m=json.loads(manifest.read_text(encoding='utf-8-sig'))
        self.flags=self.record['checks'];require(self.flags and all(type(v) is bool for v in self.flags.values()),'Literal native audit boolean flags required')
        failed=[k for k,v in self.flags.items() if not v];aggregate=self.record['all_primary_checks_passed'] if name=='fluent' else self.record['all_checks_passed'];require(type(aggregate) is bool and aggregate==(not failed) and set(self.record.get('failed_checks',failed))==set(failed),'Raw audit aggregate/flags differ')
        require(not self.record.get('errors',{}),'Native/parser errors block data selection')
        rows=m['inputs'] if 'inputs' in m else m['sources'];self.bound={}
        for row in rows:
            p=str(Path(row['path']).resolve())
            if 'sha256_before' in row:
                require(row['sha256_before']==row['sha256_after'],'Native source changed during accepted audit');h=row['sha256_before']
            else:h=row['sha256']
            require(p not in self.bound,'Repeated native audit input');self.bound[p]=h
        stable=self.flags['source_bytes_unchanged_during_read'] if name=='fluent' else m['inputs_unchanged_during_analysis'];require(stable is True,'Native audit does not prove stable input hashes')
        require(sha(fields)==expected_sha==self.record['audited_fields_npz_sha256'],'Selected native fields not bound to accepted audit')
        output_binding=m['outputs'][fields.name] if name=='fluent' else m['audited_outputs'][fields.name]['sha256'];require(output_binding==expected_sha,'Source manifest/NPZ binding differs')
        self.used={str(p.resolve()):sha(p) for p in (audit,manifest,fields)}
        state_results=self.record.get('state_results',[]);state_misses=[]
        for row in state_results:
            require(all(type(v) is bool for v in row['checks'].values()),'Native state flags must retain literal booleans');state_misses.extend(f"step {row['step']}: {key}" for key,value in row['checks'].items() if not value)
        self.status=dict(checks=self.flags,failed_checks=failed,failed_state_checks=state_misses,all_declared_checks_passed=not failed,audit_path=str(audit),manifest_path=str(manifest),fields_sha256=expected_sha,native_state_results=state_results)
    def consume(self,path):
        p=Path(path).resolve();key=str(p);require(p.is_file() and key in self.bound,'Selected native input not audit-bound: '+key)
        if key not in self.used:require(sha(p)==self.bound[key],'Selected consumed native input changed: '+key);self.used[key]=self.bound[key]
        return p


def load_history(name,fields,sources,g,phase):
    z=np.load(fields,allow_pickle=False);n=6 if phase=='pilot' else 401
    b={k:z[k].copy() for k in ('time_s','step_index','U_source_history_m_s','p_source_history_Pa','source_centres_m','source_areas_m2','source_native_cell_index','dimensionless_divergence_max','dimensionless_divergence_RMS','angular_momentum_per_depth','solved_time_s')}
    require(b['time_s'].shape==b['step_index'].shape==(n,) and np.array_equal(b['step_index'],np.arange(n)) and np.max(abs(b['time_s']-np.arange(n)*.05))<1e-10,'Native saved step/time schedule differs')
    require(b['U_source_history_m_s'].shape==(n,1536,2) and b['p_source_history_Pa'].shape==(n,1536),'Actual native field history shape differs')
    require(all(np.isfinite(v).all() for v in b.values()) and b['source_centres_m'].shape==(1536,2) and b['source_areas_m2'].shape==(1536,),'Native common arrays must be finite and complete')
    require(np.max(abs(b['source_centres_m']-g['centres']))<1e-12 and np.max(abs(b['source_areas_m2']-g['areas']))<1e-14 and np.array_equal(np.sort(b['source_native_cell_index']),np.arange(1536)),'Actual native geometry/source mapping differs')
    require(np.array_equal(b['solved_time_s'],b['time_s'][1:]),'Unsolved initial state must not have a residual row')
    for key in ('dimensionless_divergence_max','dimensionless_divergence_RMS','angular_momentum_per_depth'):require(b[key].shape==(n,),'Native history summary shape differs')
    U=b['U_source_history_m_s'];xy=b['source_centres_m'];r=np.linalg.norm(xy,axis=1);b.update(name=name,sources=sources,speed=np.linalg.norm(U,axis=2),radius_m=r,vtheta=(-xy[None,:,1]*U[:,:,0]+xy[None,:,0]*U[:,:,1])/r[None,:],vr=np.sum(xy[None,:,:]*U,axis=2)/r[None,:])
    L=np.sum(b['source_areas_m2'][None,:]*(xy[None,:,0]*U[:,:,1]-xy[None,:,1]*U[:,:,0]),axis=1)
    require(np.max(abs(L-b['angular_momentum_per_depth']))<1e-14,'Actual native-centre angular momentum differs from audited summary')
    b['steady_assessment']=dict(assessed=False,numerically_steady_for_this_run=False)
    reference=-r/3+1/(300*r);difference=b['vtheta'][-1]-reference;b['circular_steady_limit_difference']=dict(max_abs_vtheta_difference_m_s=float(abs(difference).max()),area_weighted_RMS_vtheta_difference_m_s=float(np.sqrt(np.average(difference*difference,weights=b['source_areas_m2']))),accuracy_gate=False,reference_is_transient_solution=False)
    if phase=='full':
        delta=U[400]-U[360];den=np.sum(b['source_areas_m2'][:,None]*U[400]**2);Lden=abs(L[400]);du=float(np.sqrt(np.sum(b['source_areas_m2'][:,None]*delta**2)/den)) if den>0 else None;dl=float(abs(L[400]-L[360])/Lden) if Lden>0 else None
        b['steady_assessment']=dict(assessed=True,comparison_steps=[360,400],velocity_final_denominator_positive=bool(den>0),angular_momentum_final_denominator_positive=bool(Lden>0),weighted_velocity_relative_L2_change=du,relative_angular_momentum_change=dl,declared_limit=.001,numerically_steady_for_this_run=du is not None and dl is not None and du<=.001 and dl<=.001)
    if name=='fluent':b['terminal_scaled_residuals']=z['terminal_scaled_residuals'].copy();require(b['terminal_scaled_residuals'].shape==(n-1,3),'Fluent native scaled residual shape differs')
    if name=='openfoam':
        b['residual_component_names']=z['residual_component_names'].astype(str).tolist();require(b['residual_component_names'] in (['p','Ux','Uy'],['p','Ux','Uy','Uz']),'Actual printed native OpenFOAM residual component names required')
        for key in ('last_outer_initial_residuals','last_outer_final_linear_residuals'):b[key]=z[key].copy();require(b[key].shape==(n-1,len(b['residual_component_names'])) and np.isfinite(b[key]).all(),'Actual native OpenFOAM final-outer residual history differs')
        b['outer_corrector_count']=z['outer_corrector_count'].copy();require(b['outer_corrector_count'].shape==(n-1,),'Actual native outer-count shape differs')
    z.close();return b


def verify_fluent_frames(b,reader_path,g,indices):
    sources=b['sources'];path=sources.consume(reader_path);require(path.name.startswith('annular_couette_fluent_common_v'),'Explicit native read-only Fluent helper required')
    spec=importlib.util.spec_from_file_location('_annular_native_reader',path);reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    rows=sources.record['native_state_paths'];require(len(rows)==len(b['time_s']),'Literal native Fluent state rows required');proof=[]
    for step in indices:
        common.resource_guard();row=rows[int(step)];require(row['index']==step and abs(row['time_s']-b['time_s'][step])<1e-10 and row['source_kind']=='native_fluent_hdf','Native Fluent source row differs')
        case=sources.consume(row['case']);data=sources.consume(row['data']);s=reader.read_state(case,data,g);require(np.array_equal(np.argsort(s['mapping']),b['source_native_cell_index']),'Actual Fluent vertex/cell bijection differs from common fields')
        U=np.column_stack((s['common']['u'],s['common']['v']));require(np.array_equal(U,b['U_source_history_m_s'][step]) and np.array_equal(s['common']['p'],b['p_source_history_Pa'][step]),'Raw selected Fluent native fields differ from audited arrays')
        t=s['native_time'];require(t['time-step']==step and abs(t['flow-time']-b['time_s'][step])<1e-10 and abs(t['physical-time-step']-.05)<1e-12,'Raw selected Fluent time/step/dt differs')
        proof.append(dict(step=int(step),case=str(case),data=str(data),native_vertex_cell_bijection=True,raw_U_p_exact=True))
    b['native_readback']=dict(source_kind='native_fluent_hdf',selected_states=proof,selected_cell_velocity_values_verified=len(indices)*1536*2,selected_cell_pressure_values_verified=len(indices)*1536)


def native_labels(path):
    count,raw=common.native.foam_list(path);require(re.fullmatch(r'[\s\d+\-]*',raw) is not None,'Native integer list required');a=np.fromstring(raw,sep=' ',dtype=np.int64);require(len(a)==count,'Native label count differs');return a


def verify_openfoam_frames(b,g,indices):
    sources=b['sources'];rows=sources.record['native_state_paths'];require(len(rows)==len(b['time_s']),'Literal native rank state rows required');cached={};proof=[]
    node_lookup={tuple(np.round(x,12)):i for i,x in enumerate(g['nodes'])};quad_lookup={tuple(sorted(q)):i for i,q in enumerate(g['quads'])}
    for step in indices:
        common.resource_guard();row=rows[int(step)];require(row['source_kind']=='native_rank_reassembly' and row['index']==step and abs(row['time_s']-b['time_s'][step])<1e-10,'Actual native rank state row differs')
        case=Path(row['native_case']).resolve()
        if case not in cached:
            mesh=common.native.foam_mesh(case)
            for p in mesh['inputs']:sources.consume(p)
            nodeids=np.array([node_lookup[tuple(np.round(x[:2],12))] for x in mesh['nodes']]);mapping=np.array([quad_lookup[tuple(np.unique(nodeids[c]))] for c in mesh['cells']]);require(np.array_equal(np.sort(mapping),np.arange(1536)) and np.max(abs(mesh['nodes'][:,:2]-g['nodes'][nodeids]))<1e-12 and np.array_equal(np.argsort(mapping),b['source_native_cell_index']),'Actual native global/source vertex-set bijection differs')
            ranks=[]
            for rank,paths in enumerate(row['rank_mesh']):
                for p in paths.values():sources.consume(p)
                rm=common.native.foam_mesh(case/f'processor{rank}');ids=native_labels(Path(paths['cellProcAddressing']));points=native_labels(Path(paths['pointProcAddressing']))
                require(len(ids)==len(rm['cells']) and len(points)==len(rm['nodes']) and np.max(abs(rm['nodes']-mesh['nodes'][points]))<1e-12,'Native rank point/cell addressing differs')
                projected=np.array([quad_lookup[tuple(np.unique(nodeids[points[c]]))] for c in rm['cells']]);require(np.array_equal(projected,mapping[ids]),'Native rank eight-vertex cell differs from source quad')
                ranks.append(ids)
            require(len(ranks)==2 and np.array_equal(np.sort(np.concatenate(ranks)),np.arange(1536)),'Both native ranks must cover every cell exactly once');cached[case]=ranks
        ranks=cached[case];U=np.empty((1536,3));p=np.empty(1536)
        for rank,ids in enumerate(ranks):
            upath=sources.consume(row['rank_U'][rank]);ppath=sources.consume(row['rank_p'][rank]);sources.consume(row['rank_phi'][rank]);U[ids]=common.native.foam_field(upath,len(ids),3,[0,1,-1,0,0,0,0]);p[ids]=common.native.foam_field(ppath,len(ids),1,[0,2,-2,0,0,0,0])
            if step:
                path=sources.consume(row['rank_time_metadata'][rank]);raw=common.native.foam_text(path);values={key:float(re.search(r'\b'+key+r'\s+([-+0-9.eE]+)\s*;',raw)[1]) for key in ('value','index','deltaT')};require(abs(values['value']-step*.05)<1e-10 and values['index']==step and abs(values['deltaT']-.05)<1e-12,'Native rank time/dt/global index differs')
        order=b['source_native_cell_index'];require(np.array_equal(U[order,:2],b['U_source_history_m_s'][step]) and np.array_equal(p[order],b['p_source_history_Pa'][step]),'Actual selected two-rank native fields differ from common arrays')
        proof.append(dict(step=int(step),rank_U=row['rank_U'],rank_p=row['rank_p'],complete_cell_vertex_bijection=True,raw_U_p_exact=True))
    b['native_readback']=dict(source_kind='native_rank_reassembly',selected_states=proof,selected_cell_velocity_values_verified=len(indices)*1536*2,selected_cell_pressure_values_verified=len(indices)*1536)


def notes(bundles):
    lines=[]
    for name,b in bundles.items():
        status=b['sources'].status;failed=status['failed_checks']+status['failed_state_checks'];text=NAMES[name]+(': primary numerical audit passed' if not failed else ': primary numerical MISSED '+', '.join(failed));lines.extend(textwrap.wrap(text,width=150))
    return lines


def map_canvas(g,subtitle,bundles=None,state=0,geometry_only=False):
    status=[] if geometry_only else notes(bundles);height=7.2+.12*max(0,len(status)-2);fig=plt.figure(figsize=(12.8,height),dpi=100)
    fig.text(.065,.958,'Annular Couette spin-up',fontsize=21,fontweight='bold');fig.text(.065,.912,subtitle,fontsize=10,color='#425e6d');clock=fig.text(.065,.866,'Geometry only | no native solver field loaded' if geometry_only else '',fontsize=10,color='#304b59')
    artists={};axes={};bottom=.215+.018*max(0,len(status)-2);side=.34;mapheight=side*12.8/height
    for left,name in ((.09,'fluent'),(.56,'openfoam')):
        ax=fig.add_axes([left,bottom,side,mapheight]);kw=dict(facecolors='#f3f6f8',edgecolors='#b4bfc5',linewidths=.22) if geometry_only else dict(array=bundles[name]['speed'][state],cmap='viridis',norm=Normalize(0,.05),edgecolors='none',antialiased=False)
        collection=PolyCollection(g['polygons'],**kw);ax.add_collection(collection)
        for radius_index in (0,24):
            xy=g['nodes'][radius_index*64:(radius_index+1)*64];xy=np.vstack((xy,xy[0]));ax.plot(xy[:,0],xy[:,1],color='#243b49',lw=.8)
        ax.annotate('',xy=(0,.031),xytext=(.031,0),arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=.35',lw=1.2,color='#b45a33'));ax.text(0,-.008,'CCW wall\nfixed mesh',ha='center',va='center',fontsize=7.5,color='#805335')
        ax.set(xlim=(-.105,.105),ylim=(-.105,.105),xlabel='x (m)',ylabel='y (m)');ax.set_aspect('equal');ax.set_xticks([-.1,-.05,0,.05,.1]);ax.set_yticks([-.1,-.05,0,.05,.1]);ax.tick_params(labelsize=8);ax.set_title(NAMES[name]+(' | planned mesh layout' if geometry_only else ' | actual native speed'),loc='left',fontsize=11,pad=9,fontweight='bold');artists[name]=collection;axes[name]=ax
    cax=fig.add_axes([.29,.125,.46,.018]);fig.colorbar(ScalarMappable(norm=Normalize(0,.05),cmap='viridis'),cax=cax,orientation='horizontal',label='Speed (m/s) | fixed common 0-0.05 m/s scale')
    if geometry_only:fig.text(.065,.026,'Actual shared polygon mesh only. White cells are geometry, not a computed velocity field.',fontsize=8.5,color='#425e6d')
    else:common.add_footer(fig,status,bottom=.018)
    return fig,artists,clock,axes


def render_animation(g,bundles,out,phase):
    indices=np.arange(6) if phase=='pilot' else np.arange(0,401,4);title='Pilot | 5 actual steps to 0.25 s' if phase=='pilot' else '400 actual steps to 20 s | initial rest, inner wall 1 rad/s, outer wall fixed'
    fig,artists,clock,axes=map_canvas(g,title,bundles);images=[];crops={}
    for j,step in enumerate(indices):
        if j%10==0:common.resource_guard()
        ts=[b['time_s'][step] for b in bundles.values()];require(max(ts)-min(ts)<1e-10,'Actual paired frame times differ')
        for name,b in bundles.items():artists[name].set_array(b['speed'][step])
        clock.set_text(f'Actual physical time {ts[0]:.2f} s | native step {step} | playback 60 ms/frame')
        if j==0:common.verify_text_bounds(fig,out,'annular-couette-actual-speed')
        fig.canvas.draw();im=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy());images.append(im)
        for name,ax in axes.items():box=ax.get_window_extent();crops[name]=[round(box.x0)+2,round(im.height-box.y1)+2,round(box.x1)-2,round(im.height-box.y0)-2]
    images[-1].save(out/'annular-couette-final-speed.png');path=out/'annular-couette-actual-speed.gif';images[0].save(path,save_all=True,append_images=images[1:],duration=60,loop=0,optimize=False,disposal=2);plt.close(fig)
    full=[];maps={name:[] for name in bundles};durations=[]
    with Image.open(path) as gif:
        require(gif.n_frames==len(indices),'Decoded actual frame count differs')
        for j in range(gif.n_frames):
            gif.seek(j);im=gif.convert('RGB');durations.append(gif.info['duration']);full.append(hashlib.sha256(im.tobytes()).hexdigest())
            for name,crop in crops.items():maps[name].append(hashlib.sha256(im.crop(crop).tobytes()).hexdigest())
        require(durations==[60]*len(indices),'Decoded GIF timing differs');sheet=Image.new('RGB',(images[0].width,(images[0].height+26)*3),'#e9eff2');draw=ImageDraw.Draw(sheet)
        for slot,frame in enumerate((0,len(indices)//2,len(indices)-1)):
            gif.seek(frame);draw.text((12,slot*(images[0].height+26)+6),f'Decoded native step {indices[frame]} | physical time {bundles["fluent"]["time_s"][indices[frame]]:.2f} s',fill='#203c4b');sheet.paste(gif.convert('RGB'),(0,slot*(images[0].height+26)+26))
        sheet.save(out/'annular-couette-decoded-contact-sheet.png')
    return dict(frames=len(indices),frame_delay_ms=60,actual_native_step_indices=indices.tolist(),actual_physical_times_s=bundles['fluent']['time_s'][indices].tolist(),actual_solver_velocity_evolution=True,temporal_interpolation=False,particles=False,wall_arrow_annotation_only=True,mesh_motion=False,common_speed_limits_m_s=[0,.05],true_equal_xy_aspect=True,distinct_decoded_RGB_frames=len(set(full)),map_region_unique_frames={k:len(set(v)) for k,v in maps.items()},map_crop_pixel_bounds=crops,map_frame_uniqueness_requirement=False,plateau_policy='Genuine settled or quantized map repeats remain unaltered',sha256=sha(path))


def render_profiles(g,bundles,out):
    fig,axes=plt.subplots(1,2,figsize=(12.8,5.6));fig.subplots_adjust(left=.09,right=.96,top=.78,bottom=.25,wspace=.29);fig.suptitle('Final radial velocity profiles',x=.09,y=.96,ha='left',fontsize=20,fontweight='bold')
    fig.text(.09,.875,'Area-weighted angular means and angular min-max bands at actual native cell-centre radii.',fontsize=9,color='#425e6d')
    for name,b in bundles.items():
        radius=[];vt=[];vr=[];vtlo=[];vthi=[];vrlo=[];vrhi=[]
        for j in range(24):
            mask=g['radial_index']==j;weights=b['source_areas_m2'][mask];radius.append(np.average(b['radius_m'][mask],weights=weights));vt.append(np.average(b['vtheta'][-1,mask],weights=weights));vr.append(np.average(b['vr'][-1,mask],weights=weights));vtlo.append(b['vtheta'][-1,mask].min());vthi.append(b['vtheta'][-1,mask].max());vrlo.append(b['vr'][-1,mask].min());vrhi.append(b['vr'][-1,mask].max())
        style='--' if name=='openfoam' else '-';axes[0].plot(radius,vt,label=NAMES[name],color=COLOURS[name],ls=style);axes[0].fill_between(radius,vtlo,vthi,color=COLOURS[name],alpha=.12);axes[1].plot(radius,vr,label=NAMES[name],color=COLOURS[name],ls=style);axes[1].fill_between(radius,vrlo,vrhi,color=COLOURS[name],alpha=.12)
    r=np.linspace(.05,.1,250);axes[0].plot(r,-r/3+1/(300*r),'k--',lw=1,label='Circular steady limit')
    axes[0].set(xlabel='Actual radius r (m)',ylabel='Tangential velocity\nvθ (m/s)');axes[1].set(xlabel='Actual radius r (m)',ylabel='Radial velocity\nvr (m/s)')
    for ax in axes:ax.grid(alpha=.17);ax.legend(frameon=False,fontsize=9);ax.spines[['top','right']].set_visible(False)
    assessed=[NAMES[n]+(': numerically steady for this run' if b['steady_assessment']['numerically_steady_for_this_run'] else ': final saved transient') for n,b in bundles.items()];fig.text(.09,.13,'; '.join(assessed)+'. Circular steady reference is descriptive; no accuracy gate.',fontsize=8.1,color='#425e6d');common.add_footer(fig,notes(bundles));common.save(fig,out,'annular-couette-radial-profiles')


def render_numerics(bundles,out):
    fig,axes=plt.subplots(2,2,figsize=(12.8,9));fig.subplots_adjust(left=.1,right=.96,top=.84,bottom=.17,wspace=.32,hspace=.49);fig.suptitle('Angular momentum and native numerical histories',x=.1,y=.96,ha='left',fontsize=19,fontweight='bold')
    for name,b in bundles.items():
        style='--' if name=='openfoam' else '-';axes[0,0].plot(b['time_s'],b['angular_momentum_per_depth'],ls=style,color=COLOURS[name],label=NAMES[name]);axes[0,1].plot(b['time_s'][1:],b['dimensionless_divergence_max'][1:],ls=style,color=COLOURS[name],label=NAMES[name]+' max');axes[0,1].plot(b['time_s'][1:],b['dimensionless_divergence_RMS'][1:],ls=':',color=COLOURS[name],label=NAMES[name]+' RMS')
    axes[0,0].set(xlabel='Actual physical time (s)',ylabel='Angular momentum per depth\nLz (kg m/s)');axes[0,1].set(xlabel='Solved physical time (s)',ylabel='Dimensionless native\nlocal divergence');axes[0,1].set_yscale('symlog',linthresh=1e-12)
    axes[0,1].axhline(1e-5,color='#4b616f',lw=.8,ls='--',label='Max target 1e-5');axes[0,1].axhline(1e-6,color='#74828b',lw=.8,ls=':',label='RMS target 1e-6')
    b=bundles['fluent']
    for j,label in enumerate(('Continuity','x momentum','y momentum')):axes[1,0].plot(b['solved_time_s'],b['terminal_scaled_residuals'][:,j],label=label)
    axes[1,0].axhline(1e-6,color='#4b616f',lw=.8,ls=':',label='Own target 1e-6');axes[1,0].set(title='ANSYS Fluent | terminal scaled residuals',xlabel='Solved physical time (s)',ylabel='Native scaled\nresidual')
    b=bundles['openfoam'];outer=np.max(b['last_outer_initial_residuals'],axis=1);linear=np.max(b['last_outer_final_linear_residuals'],axis=1);axes[1,1].plot(b['solved_time_s'],outer,label='Final outer: max initial p/U');axes[1,1].plot(b['solved_time_s'],linear,ls='--',label='Final outer: max final linear p/U');axes[1,1].axhline(1e-6,color='#4b616f',lw=.8,ls=':',label='Outer target 1e-6');axes[1,1].axhline(1e-8,color='#74828b',lw=.8,ls='-.',label='Linear target 1e-8');axes[1,1].set(title='OpenFOAM | distinct outer and linear residuals',xlabel='Solved physical time (s)',ylabel='Native own-norm\nresidual')
    for ax in axes[1]:ax.set_yscale('symlog',linthresh=1e-14)
    for ax in axes.ravel():ax.grid(alpha=.17);ax.legend(frameon=False,fontsize=7.5);ax.spines[['top','right']].set_visible(False)
    fig.text(.1,.09,'Local divergence uses actual native face fluxes. Closed zero net flux alone is not proof; solver residual norms differ.',fontsize=8.2,color='#425e6d');common.add_footer(fig,notes(bundles));common.save(fig,out,'annular-couette-native-histories')


def export_csvs(g,bundles,out,indices):
    header=['solver','native_step','physical_time_s','source_cell_0based','x_m','y_m','u_m_s','v_m_s','native_gauge_pressure_Pa','pressure_minus_area_mean_Pa']
    def rows(name,b,states):
        for step in states:
            U=b['U_source_history_m_s'][step];p=b['p_source_history_Pa'][step];mean=np.average(p,weights=b['source_areas_m2'])
            for cell in range(1536):yield [name,int(step),b['time_s'][step],cell,*b['source_centres_m'][cell],*U[cell],p[cell],p[cell]-mean]
    with gzip.open(out/'actual-native-cell-fields.csv.gz','wt',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(header)
        for name,b in bundles.items():w.writerows(rows(name,b,range(len(b['time_s']))))
    selected=np.array([0,1,4,5]) if len(bundles['fluent']['time_s'])==6 else np.array([0,4,40,100,360,400])
    with (out/'record-selected-native-fields.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(header)
        for name,b in bundles.items():w.writerows(rows(name,b,selected))
    with (out/'actual-step-history.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','native_step','physical_time_s','angular_momentum_per_depth_kg_m_s','dimensionless_divergence_max','dimensionless_divergence_RMS','speed_min_m_s','speed_max_m_s'])
        for name,b in bundles.items():
            for step in range(len(b['time_s'])):w.writerow([name,step,b['time_s'][step],b['angular_momentum_per_depth'][step],b['dimensionless_divergence_max'][step],b['dimensionless_divergence_RMS'][step],b['speed'][step].min(),b['speed'][step].max()])
    with (out/'actual-native-residuals.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','native_step','physical_time_s','native_norm','component','actual_residual'])
        for name,b in bundles.items():
            groups=[('terminal_scaled',('continuity','x-velocity','y-velocity'),b['terminal_scaled_residuals'])] if name=='fluent' else [(key,b['residual_component_names'],b[key]) for key in ('last_outer_initial_residuals','last_outer_final_linear_residuals')]
            for norm,labels,values in groups:
                for i,t in enumerate(b['solved_time_s']):
                    for j,label in enumerate(labels):w.writerow([name,i+1,t,norm,label,values[i,j]])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--phase',choices=('geometry','pilot','full'),required=True);ap.add_argument('--geometry',type=Path,required=True);ap.add_argument('--targets',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--fluent-reader',type=Path)
    for name in NAMES:
        for suffix in ('audit','manifest','fields'):ap.add_argument('--'+name+'-'+suffix,type=Path)
        ap.add_argument('--'+name+'-fields-sha256')
    a=ap.parse_args();out=a.out.resolve();require(out.drive.upper()=='E:' and not out.exists(),'Fresh E: output directory required');require(sha(a.targets)==TARGETS_SHA,'Frozen model targets changed');require(sha(Path(common.__file__))==FROZEN_HELPERS['cavity'] and sha(Path(common.native.__file__))==FROZEN_HELPERS['native'],'Frozen native reader/layout helper changed');base_inputs={str(p.resolve()):sha(p) for p in (Path(__file__),Path(common.__file__),Path(common.native.__file__),a.geometry,a.targets)};g=geometry(a.geometry);resources=common.resource_guard(render_reserve=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.7})
    if a.phase=='geometry':
        out.mkdir(parents=True);fig,_,_,_=map_canvas(g,'Geometry-only preparation | actual 64 × 24 polar quad mesh',geometry_only=True);common.save(fig,out,'annular-couette-geometry-only-layout');(out/'geometry-only-proof.json').write_text(json.dumps(dict(actual_mesh_cells=1536,actual_mesh_nodes=1600,holes=1,geometry_sha256=sha(a.geometry),plotter_sha256=sha(Path(__file__)),solver_fields_loaded=False,actual_polygon_area_m2=g['actual_polygon_area_m2'],relative_area_difference_from_circular_annulus=g['relative_area_difference_from_circular_annulus'],resources=resources),indent=2)+'\n');return
    bundles={}
    for name in NAMES:
        audit,manifest,fields=[getattr(a,name+'_'+suffix) for suffix in ('audit','manifest','fields')];require(all(p is not None for p in (audit,manifest,fields)),'Explicit native selected audit/manifest/fields required');source=AcceptedSources(name,audit,manifest,fields,getattr(a,name+'_fields_sha256'));bundles[name]=load_history(name,fields,source,g,a.phase)
    from ansys.fluent.core.offline import _lispy
    base_inputs[str(Path(_lispy.__file__).resolve())]=sha(Path(_lispy.__file__))
    indices=np.arange(6) if a.phase=='pilot' else np.arange(0,401,4);require(a.fluent_reader is not None,'Explicit audit-bound Fluent native reader required');verify_fluent_frames(bundles['fluent'],a.fluent_reader,g,indices);verify_openfoam_frames(bundles['openfoam'],g,indices)
    out.mkdir(parents=True);animation=render_animation(g,bundles,out,a.phase);render_profiles(g,bundles,out);render_numerics(bundles,out);export_csvs(g,bundles,out,indices)
    inputs=dict(base_inputs)
    for b in bundles.values():inputs.update(b['sources'].used)
    require(all(sha(Path(p))==h for p,h in inputs.items()),'Consumed native/source input changed during export')
    checks=dict(phase=a.phase,actual_native_fields=True,numerical_audits={n:b['sources'].status for n,b in bundles.items()},native_mapping={n:b['native_readback'] for n,b in bundles.items()},steady_assessment={n:b['steady_assessment'] for n,b in bundles.items()},descriptive_circular_steady_limit_difference={n:b['circular_steady_limit_difference'] for n,b in bundles.items()},animation=animation,input_hashes_stable=True,accepted_full_audit_not_repeated=True,independently_consumed_source_count=len(inputs),actual_polygon_area_m2=g['actual_polygon_area_m2'],polygon_vs_circle_relative_area_difference=g['relative_area_difference_from_circular_annulus'],analytic_reference=dict(formula='vtheta=A*r+B/r',A_per_s=-1/3,B_m2_per_s=1/300,steady_limit_only=True,accuracy_gate=False),pressure_gauge_policy='Native raw Pa retained; any cross-code pressure comparison subtracts each actual field area/volume-weighted mean',torque_balance_claim=False,physical_validation=False,mesh_independence=False,time_step_independence=False,three_dimensional_Taylor_vortex_stability_claim=False,resources=dict(before=resources,after=common.resource_guard()))
    checks['diagnostic_comparisons']=bundles['openfoam']['sources'].record.get('diagnostic_comparisons',{})
    checks['execution_completion']=bundles['openfoam']['sources'].record.get('execution_completion',{})
    (out/'checks.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n',encoding='utf-8');(out/'input-hash-manifest.json').write_text(json.dumps(dict(inputs=[dict(path=p,sha256=h) for p,h in sorted(inputs.items())]),indent=2)+'\n',encoding='utf-8');(out/'executed-plotter.py').write_bytes(Path(__file__).read_bytes())
    caption=f"Actual native transient Couette velocity fields on a fixed polygonal annulus. This {a.phase} export selects {len(indices)} actual saved states from step {indices[0]} to {indices[-1]}, physical time {bundles['fluent']['time_s'][indices[0]]:.2f} to {bundles['fluent']['time_s'][indices[-1]]:.2f} s, at 60 ms per playback frame. No interpolated times or particles. The CCW arrow annotates the prescribed wall rotation on a stationary mesh. Genuine settled or quantized map repeats are preserved. Analytical A*r+B/r is a circular steady limit, never a startup field. Final numerically steady wording requires both actual 360-to-400 weighted velocity change with final denominator and relative angular-momentum change <=1e-3. Native local continuity uses audited actual face fluxes; closed zero net boundary flux alone is not proof. Raw native pressure and separately mean-subtracted pressure are labeled in the CSV. No analytical or cross-code accuracy gate, physical validation, mesh/time independence or 3D stability claim.\n\n"+'\n'.join(notes(bundles))+'\n'
    d=checks['diagnostic_comparisons'];completion=checks['execution_completion']
    if d:
        mutations=d['derived_restart_output_mutations'];maximum=max(row['native_internal_max_abs_difference_s_inverse'] for row in mutations)
        caption+=f"\nOpenFOAM artifact diagnostics remain separate from primary numerical checks: original_strict_all_prepared_inputs_unchanged={d['original_strict_all_prepared_inputs_unchanged']}; raw_wrapper_postcheck_passed={d['raw_wrapper_postcheck_passed']}. {d['unchanged_prepared_input_count']}/{d['prepared_input_count']} prepared inputs remained byte-identical; {len(mutations)} copied derived restart vorticity files changed, with maximum native internal difference {maximum:.17g} s^-1. original_accepted_pilot_unchanged={d['original_accepted_pilot_unchanged']}; primary_restart_inputs_unchanged={d['primary_restart_inputs_unchanged']}. Actual native solver completion={completion['actual_solver_completed']}; raw wrapper stage={completion['raw_wrapper_status']['stage']}, exit code={completion['raw_wrapper_status']['exit_code']}. Exact raw diagnostic flags, paths, hashes, and completion evidence are retained in checks.json.\n"
    (out/'CAPTIONS.md').write_text(caption,encoding='utf-8');files=sorted(p for p in out.iterdir() if p.is_file());(out/'output-hash-manifest.json').write_text(json.dumps([dict(path=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in files],indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(output=str(out),animation=animation,steady_assessment=checks['steady_assessment']),indent=2))


if __name__=='__main__':main()
