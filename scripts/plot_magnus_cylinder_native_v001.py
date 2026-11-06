"""Rotating-cylinder startup: actual native states and facet pressures.

Geometry mode reads only the shared mesh. Pilot/full modes require explicit
accepted audit selections; no solver, reconstruction or native write is invoked.
"""
from pathlib import Path
from typing import TypedDict
import argparse,csv,gzip,hashlib,importlib.util,json,os,re,sys,textwrap
os.environ['TEMP']=os.environ['TMP']='E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-magnus'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Polygon
from PIL import Image,ImageDraw
import plot_cavity_native_v003 as common
from plot_annular_couette_native_v004 import AcceptedSources,native_labels
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-magnus'
require,sha=common.require,common.sha
N=8192
NAMES={'fluent':'ANSYS Fluent','openfoam':'OpenFOAM'}
COLOURS={'fluent':'#008498','openfoam':'#d65d28'}
GEOMETRY_SHA='7490cd34096043a789a8a7abdc01ea806f4f533a9b9d437ac63a8158c190092a'
TARGETS_SHA='87dfdeb12300ba0d6dfaefdedf58d0c08ef984673522f4432b1fe90807f0dbbe'
FLUENT_READER_SHA='1a61e08a2f3624e6ccd36b3c5f244c67d3cbd89eee4865c2a67d030378905820'
OPENFOAM_READER_SHA='587991b8159edce429e2133170938bbc7c40ece58a0197f1ef63ba0147c3490f'
PLAN_SHA='0e811eb689652575cdb036560b65a2a3f3c0cc2abd604b0c6073bac8a390d49d'
SETUP_SHA='936c01fe314266f1acd9c3c26640801b7c5f784ed2fb6337fe990ef33f1711f8'
HELPERS={'plot_cavity_native_v003.py':'c2425654f684bc8de8d9219dffa872d8b3c80a605bb35fce2b9070362f388445','plot_nozzle_native_v004.py':'4596d6a6b975abc6110d80fdb8b5b680878cdff18022ce969866b288aacb2f02','plot_annular_couette_native_v004.py':'d5c4414beeed4bc6cb78795a06bf8d7e60445a476cd43eaa2ea9f818ce92884b'}
NEAR_VIEW=(-.35,.65,-.30,.30)
PREFIX='magnus-cylinder'


class NativeHistory(TypedDict):
    """Root-agreed common arrays; actual native residual norms remain separate."""
    time_s:np.ndarray
    step_index:np.ndarray
    U_source_history_m_s:np.ndarray
    p_source_history_Pa:np.ndarray
    source_centres_m:np.ndarray
    source_areas_m2:np.ndarray
    source_native_cell_index:np.ndarray
    dimensionless_divergence_max:np.ndarray
    dimensionless_divergence_RMS:np.ndarray
    relative_boundary_imbalance:np.ndarray
    native_flux_Courant_max:np.ndarray
    solved_time_s:np.ndarray
    body_pressure_force_coefficient_history:np.ndarray
    body_viscous_force_coefficient_history:np.ndarray
    body_total_force_coefficient_history:np.ndarray


def geometry(path):
    require(sha(path)==GEOMETRY_SHA,'Selected source geometry changed')
    with np.load(path,allow_pickle=False) as z:g={k:z[k] for k in z.files}
    require(g['nodes'].shape==(8320,2) and g['quads'].shape==(N,4) and np.array_equal(g['grid_shape'],[64,128]),'Actual rotating-cylinder O-grid size differs')
    poly=g['nodes'][g['quads']];x,y=poly[:,:,0],poly[:,:,1];cross=x*np.roll(y,-1,axis=1)-np.roll(x,-1,axis=1)*y;area=cross.sum(axis=1)/2
    centres=np.column_stack([((poly[:,:,i]+np.roll(poly[:,:,i],-1,axis=1))*cross).sum(axis=1) for i in (0,1)])/(6*area[:,None])
    require(np.all(area>0) and np.max(abs(area-g['areas']))<1e-14 and np.max(abs(centres-g['centres']))<1e-12,'Actual polygon area/centre differs')
    require(g['face_nodes'].shape==(16512,2) and 8320-16512+N==0 and len(np.unique(g['nodes'],axis=0))==8320,'Actual O-grid seam/connectivity differs')
    require(np.array_equal(g['radial_index'],np.repeat(np.arange(64),128)) and np.array_equal(g['sector_index'],np.tile(np.arange(128),64)),'Actual radial/sector source ordering differs')
    wall=np.flatnonzero(g['face_zone']==4);order=np.argsort(g['face_source_ids'][wall]);wall=wall[order]
    require(np.array_equal(g['face_source_ids'][wall],np.arange(16257,16385)),'Actual cylinder source face IDs differ')
    g.update(polygons=poly,body_outline=np.vstack((g['nodes'][:128],g['nodes'][0])),outer_outline=np.vstack((g['nodes'][-128:],g['nodes'][-128])),surface_source_indices=wall,actual_polygon_area_m2=float(area.sum()))
    return g

def load_history(name,fields,sources,g,phase)->NativeHistory:
    n=6 if phase=='pilot' else 501
    keys=tuple(NativeHistory.__annotations__)
    with np.load(fields,allow_pickle=False) as z:
        b={k:z[k] for k in keys}
        require(b['time_s'].shape==b['step_index'].shape==(n,) and np.array_equal(b['step_index'],np.arange(n)) and np.max(abs(b['time_s']-np.arange(n)*.02))<1e-10,'Actual native state/time schedule differs')
        require(b['U_source_history_m_s'].shape==(n,N,2) and b['p_source_history_Pa'].shape==(n,N),'Actual native U/p history shape differs')
        require(b['source_centres_m'].shape==(N,2) and b['source_areas_m2'].shape==(N,) and b['source_native_cell_index'].shape==(N,),'Actual native geometry array shape differs')
        require(np.array_equal(np.sort(b['source_native_cell_index']),np.arange(N)) and np.max(abs(b['source_centres_m']-g['centres']))<1e-12 and np.max(abs(b['source_areas_m2']-g['areas']))<1e-14,'Actual native/source geometry or bijection differs')
        for key in keys:
            if key!='relative_boundary_imbalance':require(np.isfinite(b[key]).all(),'Nonfinite physical common array: '+key)
        for key in ('dimensionless_divergence_max','dimensionless_divergence_RMS','relative_boundary_imbalance','native_flux_Courant_max'):require(b[key].shape==(n,),'Native summary history shape differs: '+key)
        ratio=b['relative_boundary_imbalance'];require(np.isfinite(ratio[1:]).all() and (np.isfinite(ratio[0]) or np.isnan(ratio[0])),'Solved native boundary balance unavailable; only unsolved initial ratio may be NaN')
        b['initial_boundary_ratio']=dict(available=bool(np.isfinite(ratio[0])),value=float(ratio[0]) if np.isfinite(ratio[0]) else None,solved_gate_applied=False)
        require(np.array_equal(b['solved_time_s'],b['time_s'][1:]),'Unsolved initialization must not acquire residual/force rows')
        for term in ('pressure','viscous','total'):require(b['body_'+term+'_force_coefficient_history'].shape==(n-1,1,2),'Actual solved per-body force shape differs')
        require(z['body_names'].astype(str).tolist()==['cylinder'] and z['force_component_names'].astype(str).tolist()==['x','y'],'Native body/component ordering differs')
        if name=='fluent':
            b['terminal_scaled_residuals']=z['terminal_scaled_residuals'];b['terminal_residual_equations']=z['terminal_residual_equations'].astype(str).tolist()
            require(b['terminal_residual_equations']==['continuity','x-velocity','y-velocity'] and b['terminal_scaled_residuals'].shape==(n-1,3) and np.isfinite(b['terminal_scaled_residuals']).all(),'Actual Fluent own residual history differs')
        else:
            b['residual_component_names']=z['residual_component_names'].astype(str).tolist();require(b['residual_component_names'] in (['p','Ux','Uy'],['p','Ux','Uy','Uz']),'Actual printed native OpenFOAM components required')
            for key in ('last_outer_initial_residuals','last_outer_final_linear_residuals'):
                b[key]=z[key];require(b[key].shape==(n-1,len(b['residual_component_names'])) and np.isfinite(b[key]).all(),'Actual OpenFOAM residual shape differs')
            b['outer_corrector_count']=z['outer_corrector_count'];require(b['outer_corrector_count'].shape==(n-1,) and np.isfinite(b['outer_corrector_count']).all(),'Actual outer count history differs')
    with np.load(fields,allow_pickle=False) as z:
        for key in ('body_surface_source_face_ids','body_surface_centres_m','body_surface_pressure_final_Pa','body_surface_native_face_ids'):b[key]=z[key]
    require(np.array_equal(b['body_surface_source_face_ids'],np.arange(16257,16385)) and b['body_surface_centres_m'].shape==(128,2) and b['body_surface_pressure_final_Pa'].shape==b['body_surface_native_face_ids'].shape==(128,),'Actual native surface array shape/source IDs differ')
    require(all(np.isfinite(b[key]).all() for key in ('body_surface_centres_m','body_surface_pressure_final_Pa','body_surface_native_face_ids')) and len(np.unique(b['body_surface_native_face_ids']))==128 and np.max(abs(b['body_surface_centres_m']-g['face_centres'][g['surface_source_indices']]))<1e-12,'Actual native surface geometry/face bijection differs')
    b['surface_Cp']=b['body_surface_pressure_final_Pa']/.005;b['surface_theta_deg']=np.mod(np.degrees(np.arctan2(b['body_surface_centres_m'][:,1],b['body_surface_centres_m'][:,0])),360.)
    b.update(name=name,sources=sources,speed=np.linalg.norm(b['U_source_history_m_s'],axis=2))
    b['speed_saturation_cell_counts']=np.count_nonzero(b['speed']>.25,axis=1)
    return b

def verify_fluent_frames(b,reader_path,g,indices):
    sources=b['sources'];path=sources.consume(reader_path);require(sha(path)==FLUENT_READER_SHA,'Selected native Fluent reader differs from frozen reviewed source')
    spec=importlib.util.spec_from_file_location('_magnus_native_fluent_reader',path);reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    rows=sources.record['native_state_paths'];require(len(rows)==len(b['time_s']),'Literal native Fluent state rows required');proof=[];mesh=None
    for step in indices:
        common.resource_guard();row=rows[int(step)];require(row['index']==step and abs(row['time_s']-b['time_s'][step])<1e-10 and row['source_kind']=='native_fluent_hdf','Actual native Fluent source row differs')
        case=sources.consume(row['case']);data=sources.consume(row['data'])
        if mesh is None:mesh=reader.map_native_mesh(case,g)
        s=reader.read_state(case,data,g,mesh=mesh)
        require(np.array_equal(np.argsort(s['mapping']),b['source_native_cell_index']),'Actual Fluent native/source topology mapping differs')
        U=np.column_stack((s['common']['u'],s['common']['v']));require(np.array_equal(U,b['U_source_history_m_s'][step]) and np.array_equal(s['common']['p'],b['p_source_history_Pa'][step]),'Raw selected Fluent U/p differs from audited arrays')
        t=s['native_time'];require(t['time-step']==step and abs(t['flow-time']-b['time_s'][step])<1e-10 and abs(t['physical-time-step']-.02)<1e-12,'Raw Fluent time/step/dt differs')
        if step==len(b['time_s'])-1:
            actual=reader.final_surface_pressure(s)
            for key,value in actual.items():require(np.array_equal(value,b[key]),'Actual Fluent final surface pressure/face map differs: '+key)
            b['surface_native_readback']=dict(all128_native_facet_pressures_exact=True,actual_saved_face_pressure=True,native_face_id_convention='Fluent1-based face IDs',source_face_ids_1based=True,cell_interpolation=False)
        proof.append(dict(step=int(step),case=str(case),data=str(data),native_vertex_cell_bijection=True,actual_saved_mesh_fingerprint_verified=True,raw_U_p_exact=True))
        del s,U
    b['native_readback']=dict(source_kind='native_fluent_hdf',selected_states=proof,source_quad_cells=N,vertex_set_bijection=True,selected_velocity_values_verified=len(indices)*N*2,selected_pressure_values_verified=len(indices)*N,streaming_mesh_policy='One mapped mesh; each selected saved-case mesh fingerprint verified before reuse')

def verify_openfoam_frames(b,g,indices):
    sources=b['sources'];rows=sources.record['native_state_paths'];require(len(rows)==len(b['time_s']),'Literal native rank state rows required');cached={};proof=[]
    node_lookup={tuple(np.round(x,12)):i for i,x in enumerate(g['nodes'])};quad_lookup={tuple(sorted(q)):i for i,q in enumerate(g['quads'])}
    for step in indices:
        common.resource_guard();row=rows[int(step)];require(row['source_kind']=='native_rank_reassembly' and row['index']==step and abs(row['time_s']-b['time_s'][step])<1e-10,'Actual native rank source row differs')
        case=Path(row['native_case']).resolve()
        if case not in cached:
            mesh=common.native.foam_mesh(case)
            for p in mesh['inputs']:sources.consume(p)
            nodeids=np.array([node_lookup[tuple(np.round(x[:2],12))] for x in mesh['nodes']]);mapping=np.array([quad_lookup[tuple(np.unique(nodeids[c]))] for c in mesh['cells']])
            require(np.array_equal(np.sort(mapping),np.arange(N)) and np.max(abs(mesh['nodes'][:,:2]-g['nodes'][nodeids]))<1e-12 and np.array_equal(np.argsort(mapping),b['source_native_cell_index']),'Actual native global/source vertex-set bijection differs')
            ranks=[]
            for rank,paths in enumerate(row['rank_mesh']):
                for p in paths.values():sources.consume(p)
                rm=common.native.foam_mesh(case/f'processor{rank}');ids=native_labels(Path(paths['cellProcAddressing']));points=native_labels(Path(paths['pointProcAddressing']))
                require(len(ids)==len(rm['cells']) and len(points)==len(rm['nodes']) and np.max(abs(rm['nodes']-mesh['nodes'][points]))<1e-12,'Actual native rank addressing differs')
                projected=np.array([quad_lookup[tuple(np.unique(nodeids[points[c]]))] for c in rm['cells']]);require(np.array_equal(projected,mapping[ids]),'Actual native rank eight-vertex cell differs from source quad');ranks.append(ids)
            require(len(ranks)==2 and np.array_equal(np.sort(np.concatenate(ranks)),np.arange(N)),'Two native ranks must cover all cells once');cached[case]=ranks
        ranks=cached[case];U=np.empty((N,3));p=np.empty(N)
        require(len(row['rank_U'])==len(row['rank_p'])==len(row['rank_phi'])==2,'Both literal native rank paths required')
        for rank,ids in enumerate(ranks):
            upath=sources.consume(row['rank_U'][rank]);ppath=sources.consume(row['rank_p'][rank]);sources.consume(row['rank_phi'][rank]);U[ids]=common.native.foam_field(upath,len(ids),3,[0,1,-1,0,0,0,0]);p[ids]=common.native.foam_field(ppath,len(ids),1,[0,2,-2,0,0,0,0])
            if step:
                path=sources.consume(row['rank_time_metadata'][rank]);raw=common.native.foam_text(path);values={key:float(re.search(r'\b'+key+r'\s+([-+0-9.eE]+)\s*;',raw)[1]) for key in ('value','index','deltaT')};require(abs(values['value']-step*.02)<1e-10 and values['index']==step and abs(values['deltaT']-.02)<1e-12,'Actual native rank clock/dt/index differs')
        order=b['source_native_cell_index'];require(np.array_equal(U[order,:2],b['U_source_history_m_s'][step]) and np.array_equal(p[order],b['p_source_history_Pa'][step]),'Actual selected two-rank U/p differs from audited arrays')
        proof.append(dict(step=int(step),rank_U=row['rank_U'],rank_p=row['rank_p'],complete_cell_vertex_bijection=True,raw_U_p_exact=True))
    b['native_readback']=dict(source_kind='native_rank_reassembly',selected_states=proof,source_quad_cells=N,vertex_set_bijection=True,selected_velocity_values_verified=len(indices)*N*2,selected_pressure_values_verified=len(indices)*N)


def verify_openfoam_surface(b,g,reader_path):
    sources=b['sources'];path=sources.consume(reader_path);require(sha(path)==OPENFOAM_READER_SHA,'Reviewed native OpenFOAM reader changed')
    spec=importlib.util.spec_from_file_location('_magnus_native_surface_reader',path);reader=importlib.util.module_from_spec(spec);sys.modules[spec.name]=reader;spec.loader.exec_module(reader)
    row=sources.record['native_state_paths'][-1];lookup={tuple(np.round(x,12)):i for i,x in enumerate(g['nodes'])};edge={tuple(sorted(e)):i for i,e in enumerate(g['face_nodes'])};result={};semantics=[]
    for rank,paths in enumerate(row['rank_mesh']):
        for value in paths.values():sources.consume(value)
        case=Path(row['native_case'])/f'processor{rank}';rm=common.native.foam_mesh(case);pointmap=np.array([lookup[tuple(np.round(x[:2],12))] for x in rm['nodes']]);addr=native_labels(Path(paths['faceProcAddressing']));owner=native_labels(Path(paths['owner']));patch=reader.patch_table(Path(paths['boundary']).parent)['cylinder'];count,raw=reader.counted(Path(paths['faces']));faces=np.array([list(map(int,s.split())) for n,s in re.findall(r'(\d+)\s*\(([^()]*)\)',raw)]);require(faces.shape==(count,4),'Actual rank face connectivity required')
        ppath=sources.consume(row['rank_p'][rank]);text=common.native.foam_text(ppath);block=reader.top(reader.blocks(text)['cylinder']);require(reader.kind(block)=='zeroGradient','Actual native cylinder pressure BC differs')
        if re.search(r'\bvalue\s+',block):pressure=np.asarray(reader.values(block,patch['count'],1,'value'));kind='actual saved cylinder patch value'
        else:
            cells=native_labels(Path(paths['cellProcAddressing']));internal=common.native.foam_field(ppath,len(cells),1,[0,2,-2,0,0,0,0]);pressure=internal[owner[patch['start']:patch['start']+patch['count']]];kind='actual native zeroGradient pressure equals native owner value'
        semantics.append(kind)
        for j,lf in enumerate(range(patch['start'],patch['start']+patch['count'])):
            projected=tuple(np.unique(pointmap[faces[lf]]));require(len(projected)==2 and projected in edge,'Actual rank cylinder facet/source edge differs');source=edge[projected];fid=int(g['face_source_ids'][source]);require(fid not in result and 16257<=fid<=16384 and addr[lf]>0,'Actual body surface face coverage/orientation differs');result[fid]=(float(pressure[j]),abs(int(addr[lf]))-1,rm['nodes'][faces[lf],:2].mean(axis=0))
    require(sorted(result)==list(range(16257,16385)),'Both ranks must cover all128 actual cylinder facets once');actual_p=np.array([result[i][0] for i in sorted(result)]);actual_ids=np.array([result[i][1] for i in sorted(result)]);actual_xy=np.array([result[i][2] for i in sorted(result)])
    require(np.array_equal(actual_p,b['body_surface_pressure_final_Pa']) and np.array_equal(actual_ids,b['body_surface_native_face_ids']) and np.max(abs(actual_xy-b['body_surface_centres_m']))<1e-12,'Actual rank surface pressure/native face IDs/centres differ from accepted arrays')
    b['surface_native_readback']=dict(all128_native_facet_pressures_exact=True,source_face_ids_1based=True,native_face_id_convention='OpenFOAM0-based global face indices',actual_pressure_semantics=semantics,cell_interpolation=False)

def notes(bundles):
    result=[]
    for name,b in bundles.items():
        s=b['sources'].status;failed=s['failed_checks']+s['failed_state_checks'];result.extend(textwrap.wrap(NAMES[name]+(': primary numerical checks passed' if not failed else ': primary numerical MISSED '+', '.join(failed)),width=155))
    return result


def map_canvas(g,subtitle,bundles=None,state=0,geometry_only=False,quantity='speed',limits=None):
    fig=plt.figure(figsize=(12.8,6.4),dpi=100);fig.text(.075,.955,'Rotating-cylinder startup',fontsize=21,fontweight='bold');fig.text(.075,.89,subtitle,fontsize=10,color='#425e6d');clock=fig.text(.075,.832,'Geometry only | no solver field loaded' if geometry_only else '',fontsize=9,color='#304b59')
    artists={};axes={};norm=Normalize(*(limits or (0,.25)));cmap='viridis' if quantity=='speed' else 'RdBu_r'
    for left,name in ((.075,'fluent'),(.555,'openfoam')):
        ax=fig.add_axes([left,.27,.4,.48]);array=None if geometry_only else bundles[name]['speed'][state] if quantity=='speed' else bundles[name]['p_source_history_Pa'][state];kw=dict(facecolors='#f3f6f8',edgecolors='#adbdc5',linewidths=.18) if geometry_only else dict(array=array,cmap=cmap,norm=norm,edgecolors='none',antialiased=False);col=PolyCollection(g['polygons'],**kw);ax.add_collection(col);ax.add_patch(Polygon(g['body_outline'],facecolor='#f0f2f3',edgecolor='#263f4c',lw=.9,zorder=3));ax.set(xlim=NEAR_VIEW[:2],ylim=NEAR_VIEW[2:],xlabel='x (m)',ylabel='y (m)');ax.set_aspect('equal');ax.tick_params(labelsize=8);ax.set_title(NAMES[name]+(' | source mesh only' if geometry_only else ' | actual native '+quantity),loc='left',fontsize=10.5,pad=9,fontweight='bold');artists[name]=col;axes[name]=ax
    if geometry_only:fig.text(.075,.07,'Geometry only. Exact 128-facet solid hole on a stationary mesh; no physical field is loaded.',fontsize=9,color='#425e6d')
    else:
        cax=fig.add_axes([.28,.136,.44,.02]);fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal',label='Speed (m/s) | fixed common scale');common.add_footer(fig,notes(bundles),bottom=.018)
    return fig,artists,clock,axes

def gif_palette():
    rgb=np.rint(matplotlib.colormaps['viridis'](np.linspace(0,1,128))[:,:3]*255).astype('uint8').tolist()
    rgb += [[v,v,v] for v in range(0,256,4)]
    rgb += [[255,255,255],[38,63,76],[66,94,109],[48,75,89],[240,242,243],[0,132,152],[214,93,40]]
    rgb += [[255,255,255]]*(256-len(rgb));palette=Image.new('P',(1,1));palette.putpalette(np.asarray(rgb,dtype=np.uint8).ravel().tolist());return palette


def render_animation(g,bundles,out,phase):
    indices=np.arange(6) if phase=='pilot' else np.arange(0,501,5);subtitle='Pilot | five actual 0.02 s steps; no full-run completion claim' if phase=='pilot' else 'ReD = 20 | clockwise wall: -2 rad/s | fixed mesh | 500 actual 0.02 s steps'
    fig,artists,clock,axes=map_canvas(g,subtitle,bundles);images=[];crops={};palette=gif_palette()
    for j,step in enumerate(indices):
        if j%10==0:common.resource_guard()
        ts=[b['time_s'][step] for b in bundles.values()];require(max(ts)-min(ts)<1e-10,'Actual paired native frame times differ')
        for name,b in bundles.items():artists[name].set_array(b['speed'][step])
        counts=[int(b['speed_saturation_cell_counts'][step]) for b in bundles.values()];saturation=f' | >0.25 m/s cells F/OF: {counts[0]}/{counts[1]}' if any(counts) else ''
        clock.set_text(f'Actual physical time {ts[0]:.2f} s | native step {step} | playback 60 ms/frame'+saturation)
        if j in (0,len(indices)-1):common.verify_text_bounds(fig,out,PREFIX+'-actual-speed-'+str(step))
        fig.canvas.draw();im=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy());images.append(im.quantize(palette=palette,dither=Image.Dither.NONE))
        for name,ax in axes.items():box=ax.get_window_extent();crops[name]=[round(box.x0)+2,round(im.height-box.y1)+2,round(box.x1)-2,round(im.height-box.y0)-2]
    fig.savefig(out/(PREFIX+'-final-speed.png'),dpi=170);fig.savefig(out/(PREFIX+'-final-speed.svg'));path=out/(PREFIX+'-actual-speed.gif');images[0].save(path,save_all=True,append_images=images[1:],duration=60,loop=0,optimize=False,disposal=2);plt.close(fig)
    require(path.stat().st_size<=15*1024**2,'Actual GIF exceeds the 15 MiB public asset target; preserve this attempt')
    full=[];maps={name:[] for name in bundles};durations=[]
    with Image.open(path) as gif:
        require(gif.n_frames==len(indices),'Decoded actual GIF state count differs')
        for j in range(gif.n_frames):
            gif.seek(j);frame=gif.convert('RGB');durations.append(gif.info['duration']);full.append(hashlib.sha256(frame.tobytes()).hexdigest())
            for name,crop in crops.items():maps[name].append(hashlib.sha256(frame.crop(crop).tobytes()).hexdigest())
        require(durations==[60]*len(indices),'Decoded actual GIF timing differs');sheet=Image.new('RGB',(images[0].width,(images[0].height+26)*3),'#e9eff2');draw=ImageDraw.Draw(sheet)
        for slot,frame in enumerate((0,len(indices)//2,len(indices)-1)):
            gif.seek(frame);draw.text((12,slot*(images[0].height+26)+6),f'Decoded actual step {indices[frame]} | physical time {bundles["fluent"]["time_s"][indices[frame]]:.2f} s',fill='#203c4b');sheet.paste(gif.convert('RGB'),(0,slot*(images[0].height+26)+26))
        sheet.save(out/(PREFIX+'-decoded-contact-sheet.png'))
    proof=dict(frames=len(indices),frame_delay_ms=60,actual_native_step_indices=indices.tolist(),actual_physical_times_s=bundles['fluent']['time_s'][indices].tolist(),actual_solver_velocity_evolution=True,temporal_interpolation=False,particles=False,streamlines=False,signed_vorticity_used=False,common_speed_limits_m_s=[0,.25],speed_above_scale_cell_counts={n:b['speed_saturation_cell_counts'][indices].tolist() for n,b in bundles.items()},true_equal_xy_aspect=True,exact_polygon_body_hole=True,distinct_decoded_RGB_frames=len(set(full)),map_region_unique_frames={k:len(set(v)) for k,v in maps.items()},map_crop_pixel_bounds=crops,map_frame_uniqueness_requirement=False,plateau_policy='Actual settled or quantized map repeats are retained',palette=dict(fixed=True,speed_colour_levels=128,dither=False),pixel_dimensions=list(images[0].size),bytes=path.stat().st_size,maximum_public_bytes=15*1024**2,sha256=sha(path))
    print(json.dumps(dict(actual_animation_ready=str(path),animation=proof)),flush=True);return proof

def decorate_axes(axes):
    for ax in np.asarray(axes).ravel():ax.grid(alpha=.17);ax.spines[['top','right']].set_visible(False);ax.legend(frameon=True,facecolor='white',framealpha=1,edgecolor='none',fontsize=7.2)


def render_pressure_surface_context(g,bundles,out):
    lower=min(float(b['p_source_history_Pa'][-1].min()) for b in bundles.values());upper=max(float(b['p_source_history_Pa'][-1].max()) for b in bundles.values());require(upper>lower,'Actual final native pressure range is degenerate');fig=plt.figure(figsize=(12.8,11));fig.text(.075,.965,'Final native pressure and cylinder surface Cp',fontsize=20,fontweight='bold');fig.text(.075,.923,f'Actual saved final time {bundles["fluent"]["time_s"][-1]:.2f} s | common range from both complete final native fields | outlet gauge 0 Pa',fontsize=9.5,color='#425e6d')
    for left,name in ((.075,'fluent'),(.555,'openfoam')):
        ax=fig.add_axes([left,.59,.4,.4*12.8*.6/11]);ax.add_collection(PolyCollection(g['polygons'],array=bundles[name]['p_source_history_Pa'][-1],norm=Normalize(lower,upper),cmap='RdBu_r',edgecolors='none',antialiased=False));ax.add_patch(Polygon(g['body_outline'],facecolor='#f0f2f3',edgecolor='#263f4c',lw=.9));ax.set(xlim=NEAR_VIEW[:2],ylim=NEAR_VIEW[2:],xlabel='x (m)',ylabel='y (m)');ax.set_aspect('equal');ax.set_title(NAMES[name]+' | actual native gauge pressure',loc='left',fontsize=10,pad=9,fontweight='bold');ax.tick_params(labelsize=8)
    cax=fig.add_axes([.28,.529,.44,.012]);fig.colorbar(ScalarMappable(norm=Normalize(lower,upper),cmap='RdBu_r'),cax=cax,orientation='horizontal',label='Gauge pressure (Pa) | shared complete-final-field range')
    cp=fig.add_axes([.09,.17,.48,.24])
    for name,b in bundles.items():
        order=np.argsort(b['surface_theta_deg']);cp.plot(b['surface_theta_deg'][order],b['surface_Cp'][order],ls='--' if name=='openfoam' else '-',color=COLOURS[name],marker='.',ms=2.7,label=NAMES[name])
    cp.set(xlabel='Actual facet-centre angle from +x, CCW (degrees)',ylabel='Surface pressure coefficient Cp',title='128 actual native surface facets | no appended endpoints');cp.set_xticks([0,90,180,270,360]);decorate_axes([cp])
    context=fig.add_axes([.72,.17,.2,.2*12.8/11]);draw_geometry(context,g,full=True);context.set_title('Full-domain geometry context',fontsize=9,pad=9);context.tick_params(labelsize=7)
    fig.text(.075,.095,'Cp = (native surface pressure - outlet 0 Pa) / 0.005 Pa. Actual computed lift direction is an outcome; no ideal-flow target.',fontsize=8.5,color='#425e6d');common.add_footer(fig,notes(bundles));common.save(fig,out,PREFIX+'-pressure-surface-context')
    return dict(common_complete_final_native_range_Pa=[lower,upper],mean_subtraction=False,surface_facets_per_solver=128,Cp_reference_pressure_Pa=0,Cp_dynamic_pressure_Pa=.005,appended_wall_points=False,ideal_flow_accuracy_gate=False)


def draw_geometry(ax,g,full=False):
    ax.add_collection(PolyCollection(g['polygons'],facecolors='#f4f7f8',edgecolors='#aebec7',linewidths=.16));ax.plot(*g['body_outline'].T,color='#263f4c',lw=.8);ax.plot(*g['outer_outline'].T,color='#263f4c',lw=.8);ax.set_aspect('equal')
    if full:
        ax.set(xlim=(-2.08,2.08),ylim=(-2.08,2.08),xlabel='x (m)',ylabel='y (m)');x0,x1,y0,y1=NEAR_VIEW;ax.plot([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],ls='--',lw=1,color='#bd6e3f')
    else:ax.set(xlim=NEAR_VIEW[:2],ylim=NEAR_VIEW[2:],xlabel='x (m)',ylabel='y (m)')


def render_geometry(g,out):
    fig=plt.figure(figsize=(12.8,7.8));fig.text(.075,.96,'Rotating cylinder | actual shared O-grid',fontsize=21,fontweight='bold');fig.text(.075,.903,'Geometry only — no CFD results | 8,192 quads | 128 exact cylinder facets | fixed stationary mesh',fontsize=10,color='#425e6d');a=fig.add_axes([.09,.2,.35,.35*12.8/7.8]);draw_geometry(a,g,full=True);a.set_title('Full circular computational boundary',fontsize=10,pad=9);b=fig.add_axes([.57,.31,.36,.36*12.8*.6/7.8]);draw_geometry(b,g,full=False);b.set_title('Planned near-cylinder field view',fontsize=10,pad=9);fig.text(.075,.07,'Dashed rectangle locates the near-cylinder view. Polygon areas and native facet centres are used; no circle-area substitute.',fontsize=8.8,color='#425e6d');common.save(fig,out,PREFIX+'-geometry-only')

def render_forces(bundles,out):
    fig,axes=plt.subplots(2,2,figsize=(12.8,9.4));fig.subplots_adjust(left=.09,right=.965,top=.85,bottom=.18,wspace=.28,hspace=.42);fig.suptitle('Native cylinder drag and lift components',x=.09,y=.96,ha='left',fontsize=20,fontweight='bold');fig.text(.09,.9,'Actual pressure, viscous and total coefficients; native depth-normalized reference 0.5 rho Uin² D depth.',fontsize=9,color='#425e6d');colours={'pressure':'#3975a9','viscous':'#b87e38','total':'#008498'}
    for row,(name,b) in enumerate(bundles.items()):
        for j,label in enumerate(('Drag coefficient Cd','Lift coefficient Cl')):
            ax=axes[row,j]
            for term in ('pressure','viscous','total'):ax.plot(b['solved_time_s'],b['body_'+term+'_force_coefficient_history'][:,0,j],color=colours[term],lw=1.2,label=term.capitalize())
            ax.set(title=NAMES[name]+' | '+('x' if j==0 else 'y')+' force',xlabel='Solved physical time (s)',ylabel=label)
    decorate_axes(axes);fig.text(.09,.105,'All solved samples retained; native time-zero outputs are separate. No prescribed lift sign, plateau or periodicity target.',fontsize=8.5,color='#425e6d');common.add_footer(fig,notes(bundles));common.save(fig,out,PREFIX+'-force-histories')

def render_numerics(bundles,out):
    fig,axes=plt.subplots(3,2,figsize=(12.8,12));fig.subplots_adjust(left=.105,right=.965,top=.86,bottom=.14,wspace=.32,hspace=.57);fig.suptitle('Native numerical histories | solved steps only',x=.105,y=.965,ha='left',fontsize=19,fontweight='bold');fig.text(.105,.918,'Initialization is prescribed rest with actual boundary flux; it has no solved conservation or residual claim.',fontsize=9,color='#425e6d')
    for name,b in bundles.items():
        t=b['solved_time_s'];style='--' if name=='openfoam' else '-';kw=dict(ls=style,color=COLOURS[name]);axes[0,0].plot(t,b['relative_boundary_imbalance'][1:],label=NAMES[name],**kw);axes[0,1].plot(t,b['dimensionless_divergence_max'][1:],label=NAMES[name]+' max',**kw);axes[0,1].plot(t,b['dimensionless_divergence_RMS'][1:],color=COLOURS[name],ls=':',label=NAMES[name]+' RMS');axes[1,0].plot(t,b['native_flux_Courant_max'][1:],label=NAMES[name],**kw)
    axes[0,0].axhline(1e-5,color='#596c76',ls=':',lw=.8,label='Own target 1e-5');axes[0,0].set(ylabel='Relative native\nboundary imbalance');axes[0,1].axhline(1e-5,color='#596c76',ls='--',lw=.8,label='Max target 1e-5');axes[0,1].axhline(1e-6,color='#87939a',ls=':',lw=.8,label='RMS target 1e-6');axes[0,1].set(ylabel='Dimensionless native\nlocal divergence');axes[1,0].set(ylabel='Maximum cell flux\nCourant number',title='Matched native flux definition | descriptive')
    b=bundles['fluent']
    for j,label in enumerate(b['terminal_residual_equations']):axes[1,1].plot(b['solved_time_s'],b['terminal_scaled_residuals'][:,j],label=label)
    axes[1,1].axhline(1e-6,color='#596c76',ls=':',lw=.8,label='Own target 1e-6');axes[1,1].set(ylabel='Fluent terminal\nscaled residual',title='ANSYS Fluent | native scaled norms')
    b=bundles['openfoam']
    for ax,key,target,title in ((axes[2,0],'last_outer_initial_residuals',1e-6,'OpenFOAM | final outer initial residuals'),(axes[2,1],'last_outer_final_linear_residuals',1e-8,'OpenFOAM | final outer linear residuals')):
        for j,label in enumerate(b['residual_component_names']):ax.plot(b['solved_time_s'],b[key][:,j],label=label)
        ax.axhline(target,color='#596c76',ls=':',lw=.8,label=f'Own target {target:.0e}');ax.set(ylabel='Native own-norm\nresidual',title=title)
    for ax in axes.ravel():
        ax.set_xlabel('Solved physical time (s)')
        if ax is not axes[1,0]:ax.set_yscale('symlog',linthresh=1e-14);ax.set_ylim(bottom=0)
    decorate_axes(axes);fig.text(.105,.075,'Local divergence uses actual all-face fluxes. Solver residual norms differ. Courant number is not a time-accuracy claim.',fontsize=8.5,color='#425e6d');common.add_footer(fig,notes(bundles));common.save(fig,out,PREFIX+'-native-numerical-histories')


def history_columns(bundles):
    header=['native_step'];values=[]
    for name,b in bundles.items():
        prefix=name+'_';parts=[b['solved_time_s']];labels=['physical_time_s']
        for key in ('relative_boundary_imbalance','dimensionless_divergence_max','dimensionless_divergence_RMS','native_flux_Courant_max'):parts.append(b[key][1:]);labels.append(key)
        for term in ('pressure','viscous','total'):
            for body in range(1):
                for j,component in enumerate(('Cd','Cl')):parts.append(b['body_'+term+'_force_coefficient_history'][:,body,j]);labels.append(f'cylinder_{term}_{component}')
        groups=[('terminal_scaled',b['terminal_residual_equations'],b['terminal_scaled_residuals'])] if name=='fluent' else [(key,b['residual_component_names'],b[key]) for key in ('last_outer_initial_residuals','last_outer_final_linear_residuals')]
        for norm,names,array in groups:
            for j,label in enumerate(names):parts.append(array[:,j]);labels.append(norm+'_'+label)
        if name=='openfoam':parts.append(b['outer_corrector_count']);labels.append('outer_corrector_count')
        header.extend(prefix+label for label in labels);values.extend(parts)
    return header,np.column_stack([bundles['fluent']['step_index'][1:]]+values)

def export_csvs(g,bundles,out):
    header=['solver','native_step','physical_time_s','source_cell_0based','x_m','y_m','u_m_s','v_m_s','native_gauge_pressure_Pa']
    def rows(name,b,states):
        for step in states:
            U=b['U_source_history_m_s'][step];p=b['p_source_history_Pa'][step]
            for cell in range(N):yield [name,int(step),b['time_s'][step],cell,*b['source_centres_m'][cell],*U[cell],p[cell]]
    with gzip.open(out/'actual-native-cell-fields.csv.gz','wt',encoding='utf-8',newline='',compresslevel=1) as stream:
        w=csv.writer(stream);w.writerow(header)
        for name,b in bundles.items():
            for state in range(len(b['time_s'])):
                if state%50==0:common.resource_guard();print(json.dumps(dict(external_CSV_solver=name,native_state=state)),flush=True)
                w.writerows(rows(name,b,[state]))
    with (out/'record-final-native-fields.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(header)
        for name,b in bundles.items():w.writerows(rows(name,b,[len(b['time_s'])-1]))
    h,v=history_columns(bundles)
    with (out/'record-solved-step-history.csv').open('x',encoding='utf-8',newline='') as stream:w=csv.writer(stream);w.writerow(h);w.writerows(v)
    with (out/'record-final-surface-pressure.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','native_step','physical_time_s','source_face_ID_1based','native_face_ID','native_face_ID_convention','x_m','y_m','theta_CCW_degrees','native_surface_pressure_Pa','Cp_outlet0_q005Pa'])
        for name,b in bundles.items():
            convention='Fluent1based' if name=='fluent' else 'OpenFOAM0based_global'
            for i in range(128):w.writerow([name,int(b['step_index'][-1]),b['time_s'][-1],int(b['body_surface_source_face_ids'][i]),int(b['body_surface_native_face_ids'][i]),convention,*b['body_surface_centres_m'][i],b['surface_theta_deg'][i],b['body_surface_pressure_final_Pa'][i],b['surface_Cp'][i]])
    return dict(full_external_field_rows=sum(len(b['time_s'])*N for b in bundles.values()),record_final_field_rows=2*N,record_solved_history_rows=len(bundles['fluent']['solved_time_s']),record_final_surface_rows=256,full_external_field_file='actual-native-cell-fields.csv.gz',compact_record_files=['record-final-native-fields.csv','record-solved-step-history.csv','record-final-surface-pressure.csv'])

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--phase',choices=('geometry','pilot','full'),required=True)
    for key in ('geometry','targets','plan','setup','out'):ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--source-sha256',required=True);ap.add_argument('--fluent-reader',type=Path);ap.add_argument('--openfoam-reader',type=Path)
    for name in NAMES:
        for suffix in ('audit','manifest','fields'):ap.add_argument('--'+name+'-'+suffix,type=Path)
        ap.add_argument('--'+name+'-fields-sha256')
    a=ap.parse_args();out=a.out.resolve();require(sha(__file__)==a.source_sha256,'Executed plotter differs from reviewed source hash');require(out.drive.upper()=='E:' and not out.exists(),'Fresh E output required');require(sha(a.targets)==TARGETS_SHA and sha(a.plan)==PLAN_SHA and sha(a.setup)==SETUP_SHA,'Frozen model targets/plan/setup changed');paths=[Path(__file__),a.geometry,a.targets,a.plan,a.setup]
    for filename,digest in HELPERS.items():p=Path(__file__).with_name(filename);require(sha(p)==digest,'Frozen reused helper changed');paths.append(p)
    inputs={str(p.resolve()):sha(p) for p in paths};g=geometry(a.geometry);resources=common.resource_guard(render_reserve=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.5,'axes.linewidth':.7})
    if a.phase=='geometry':
        out.mkdir(parents=True);render_geometry(g,out);(out/'geometry-only-proof.json').write_text(json.dumps(dict(solver_fields_loaded=False,native_solver_launched=False,cells=N,nodes=8320,holes=1,actual_polygon_area_m2=g['actual_polygon_area_m2'],near_view_xy_m=list(NEAR_VIEW),source_sha256=sha(__file__),geometry_sha256=sha(a.geometry),resources=resources),indent=2)+'\n',encoding='utf-8');return
    bundles={}
    for name in NAMES:
        audit,manifest,fields=[getattr(a,name+'_'+key) for key in ('audit','manifest','fields')];digest=getattr(a,name+'_fields_sha256');require(all(p is not None for p in (audit,manifest,fields)) and digest is not None,'Explicit accepted native selections required');s=AcceptedSources(name,audit,manifest,fields,digest);bundles[name]=load_history(name,fields,s,g,a.phase)
    require(a.fluent_reader is not None and a.openfoam_reader is not None,'Explicit audit-bound native readers required');indices=np.arange(6) if a.phase=='pilot' else np.arange(0,501,5)
    from ansys.fluent.core.offline import _lispy
    inputs[str(Path(_lispy.__file__).resolve())]=sha(_lispy.__file__);verify_fluent_frames(bundles['fluent'],a.fluent_reader,g,indices);verify_openfoam_frames(bundles['openfoam'],g,indices);verify_openfoam_surface(bundles['openfoam'],g,a.openfoam_reader)
    out.mkdir(parents=True);animation=render_animation(g,bundles,out,a.phase);pressure=render_pressure_surface_context(g,bundles,out);render_forces(bundles,out);render_numerics(bundles,out);print(json.dumps(dict(science_PNGs_ready=str(out))),flush=True);csv_proof=export_csvs(g,bundles,out)
    for b in bundles.values():inputs.update(b['sources'].used)
    require(all(sha(Path(p))==digest for p,digest in inputs.items()),'Consumed native/source input changed during export')
    checks=dict(phase=a.phase,actual_native_fields=True,numerical_claim_scope='Primary native numerical checks; artifact diagnostics separate',numerical_audits={n:b['sources'].status for n,b in bundles.items()},native_mapping={n:b['native_readback'] for n,b in bundles.items()},native_surface_readback={n:b['surface_native_readback'] for n,b in bundles.items()},animation=animation,pressure_surface=pressure,input_hashes_stable=True,accepted_full_audit_not_repeated=True,independently_consumed_source_count=len(inputs),common_cell_count=N,exact_body_holes=1,actual_polygon_area_m2=g['actual_polygon_area_m2'],near_view_xy_m=list(NEAR_VIEW),initial_boundary_ratio={n:b['initial_boundary_ratio'] for n,b in bundles.items()},initial_solved_conservation_claim=False,initial_solved_forces_available=False,solved_only_force_history=True,raw_native_gauge_pressure_Pa=True,pressure_mean_subtraction=False,speed_scale_saturation_counts_all_native_states={n:b['speed_saturation_cell_counts'].tolist() for n,b in bundles.items()},csv=csv_proof,openfoam_diagnostic_comparisons=bundles['openfoam']['sources'].record.get('diagnostic_comparisons',{}),openfoam_execution_completion=bundles['openfoam']['sources'].record.get('execution_completion',{}),physical_validation=False,mesh_independence=False,time_step_independence=False,periodicity_established=False,steady_force_plateau_claim=False,prescribed_lift_direction_gate=False,signed_vorticity_used=False,resources=dict(before=resources,after=common.resource_guard()))
    (out/'checks.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n',encoding='utf-8');(out/'input-hash-manifest.json').write_text(json.dumps(dict(inputs=[dict(path=p,sha256=h) for p,h in sorted(inputs.items())]),indent=2)+'\n',encoding='utf-8');(out/'executed-plotter.py').write_bytes(Path(__file__).read_bytes())
    caption=f"Actual native viscous rotating-cylinder startup, Re_D20, clockwise wall Omega_z=-2 rad/s on a stationary mesh. This {a.phase} export selects {len(indices)} actual states, steps {indices[0]} to {indices[-1]}, physical times {bundles['fluent']['time_s'][indices[0]]:.2f} to {bundles['fluent']['time_s'][indices[-1]]:.2f} s, at60 ms playback per frame. Near-cylinder view x[-0.35,0.65], y[-0.30,0.30] m uses equal xy aspect, exact128-facet polygon hole and fixed common speed0..0.25 m/s. Saturated native cell counts are retained/labeled. No temporal interpolation, particles, synthetic signed vorticity or forced map uniqueness. Genuine repeated maps remain. Fixed palette affects preview colors only.\n\nFinal native gauge-pressure maps share one range from both complete final native fields, with no mean subtraction. Surface Cp uses actual saved/native-boundary pressure at128 native facet centres, outlet0 Pa and q=0.005 Pa; no interpolated or appended wall points. Fluent native face IDs are1-based; OpenFOAM global native face indices are0-based; canonical source face IDs are1-based16257..16384. Any saved OpenFOAM zeroGradient patch without a serialized value uses its exact native owner-pressure semantics, not interpolation. Full-domain uncolored geometry provides context.\n\nInitialization is prescribed rest with actual active inlet/rotating-wall boundary fields, not solved conservation. Unavailable initial ratios stay unavailable; mass/divergence/residual gates and force histories begin at solved step1. Raw native time-zero function outputs remain separate. Pressure/viscous/total force coefficients use each code's actual depth and shared rho/U/D reference. Actual computed lift direction is an outcome, never a gate. Native residual norms differ; cell flux Courant is descriptive. No ideal-flow equality, plateau, periodicity, physical validation or mesh/time independence is claimed.\n\n"+'\n'.join(notes(bundles))+'\n'
    for name,b in bundles.items():
        for key in ('diagnostic_comparisons','execution_completion'):
            if b['sources'].record.get(key):caption+='\n'+NAMES[name]+' exact native '+key+':\n```json\n'+json.dumps(b['sources'].record[key],indent=2,allow_nan=False)+'\n```\n'
    (out/'CAPTIONS.md').write_text(caption,encoding='utf-8');files=sorted(p for p in out.iterdir() if p.is_file());(out/'output-hash-manifest.json').write_text(json.dumps([dict(path=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in files],indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(output=str(out),animation=animation,csv=csv_proof),indent=2))

if __name__=='__main__':main()
