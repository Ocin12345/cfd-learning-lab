"""Offset tandem-square wake: actual native states, two holes, no invented motion.

Geometry mode reads only the shared mesh. Pilot/full modes require explicit
accepted audit selections; no solver, reconstruction or native write is invoked.
"""
from pathlib import Path
from typing import TypedDict
import argparse,csv,gzip,hashlib,importlib.util,json,os,re,sys,textwrap
os.environ['TEMP']=os.environ['TMP']='E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-tandem'
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
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-tandem'
require,sha=common.require,common.sha
N=9472
NAMES={'fluent':'ANSYS Fluent','openfoam':'OpenFOAM'}
COLOURS={'fluent':'#008498','openfoam':'#d65d28'}
GEOMETRY_SHA='099165a86f06dd6566c6d76860f5a0935fba4b956e16d6908f6f1b735270b178'
TARGETS_SHA='426448b7e6909e75dbf08f065d93f90699b98f12be007891631d9ed9b10f71d7'
FLUENT_READER_SHA='46a5b53aad554a6072a7407a78dc21a9c6a2c687b89b5295cbb85d6b9b3d832a'
HELPERS={'plot_cavity_native_v003.py':'c2425654f684bc8de8d9219dffa872d8b3c80a605bb35fce2b9070362f388445','plot_nozzle_native_v004.py':'4596d6a6b975abc6110d80fdb8b5b680878cdff18022ce969866b288aacb2f02','plot_annular_couette_native_v004.py':'d5c4414beeed4bc6cb78795a06bf8d7e60445a476cd43eaa2ea9f818ce92884b'}
PROFILE_COLUMNS=(52,84,112)
PREFIX='tandem-offset-wake'


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
    require(g['nodes'].shape==(9723,2) and g['quads'].shape==(N,4) and np.array_equal(g['grid_shape'],[60,160]),'Actual tandem mesh size differs')
    poly=g['nodes'][g['quads']];x,y=poly[:,:,0],poly[:,:,1];cross=x*np.roll(y,-1,axis=1)-np.roll(x,-1,axis=1)*y;area=cross.sum(axis=1)/2
    centres=np.column_stack([((poly[:,:,i]+np.roll(poly[:,:,i],-1,axis=1))*cross).sum(axis=1) for i in (0,1)])/(6*area[:,None])
    require(np.all(area>0) and abs(area.sum()-1.48)<1e-12 and np.max(abs(area-g['areas']))<1e-14 and np.max(abs(centres-g['centres']))<1e-12,'Source polygon areas/centres differ')
    require(g['face_nodes'].shape==(19196,2) and 9723-19196+N==-1 and int(g['fluid_mask'].sum())==N and int(g['solid_mask'].sum())==128,'Actual two-hole connectivity/mask differs')
    require(np.array_equal(g['slot_to_source_cell'][g['iy'],g['ix']],np.arange(N)),'Source sparse grid/source-cell bijection differs')
    outlines={}
    for zone in (8,9):
        edges=g['face_nodes'][g['face_zone']==zone];adj={}
        require(edges.shape==(32,2),'Actual square boundary edge count differs')
        for a,b in edges:adj.setdefault(int(a),[]).append(int(b));adj.setdefault(int(b),[]).append(int(a))
        require(len(adj)==32 and all(len(v)==2 for v in adj.values()),'Actual body boundary is not a closed loop')
        loop=[min(adj)];previous=None
        for unused in range(32):
            nxt=next(v for v in adj[loop[-1]] if v!=previous);previous=loop[-1];loop.append(nxt)
        require(loop[-1]==loop[0] and len(set(loop[:-1]))==32,'Actual body loop does not close exactly once');outlines[zone]=g['nodes'][loop]
    for ix in PROFILE_COLUMNS:require(np.count_nonzero(g['ix']==ix)==60,'Native profile column intersects a hole or leaves the domain')
    g.update(polygons=poly,body_outlines=outlines,actual_polygon_area_m2=float(area.sum()))
    return g


def load_history(name,fields,sources,g,phase)->NativeHistory:
    n=6 if phase=='pilot' else 801
    keys=tuple(NativeHistory.__annotations__)
    with np.load(fields,allow_pickle=False) as z:
        b={k:z[k] for k in keys}
        require(b['time_s'].shape==b['step_index'].shape==(n,) and np.array_equal(b['step_index'],np.arange(n)) and np.max(abs(b['time_s']-np.arange(n)*.05))<1e-10,'Actual native state/time schedule differs')
        require(b['U_source_history_m_s'].shape==(n,N,2) and b['p_source_history_Pa'].shape==(n,N),'Actual native U/p history shape differs')
        require(b['source_centres_m'].shape==(N,2) and b['source_areas_m2'].shape==(N,) and b['source_native_cell_index'].shape==(N,),'Actual native geometry array shape differs')
        require(np.array_equal(np.sort(b['source_native_cell_index']),np.arange(N)) and np.max(abs(b['source_centres_m']-g['centres']))<1e-12 and np.max(abs(b['source_areas_m2']-g['areas']))<1e-14,'Actual native/source geometry or bijection differs')
        for key in keys:
            if key!='relative_boundary_imbalance':require(np.isfinite(b[key]).all(),'Nonfinite physical common array: '+key)
        for key in ('dimensionless_divergence_max','dimensionless_divergence_RMS','relative_boundary_imbalance','native_flux_Courant_max'):require(b[key].shape==(n,),'Native summary history shape differs: '+key)
        ratio=b['relative_boundary_imbalance'];require(np.isfinite(ratio[1:]).all() and (np.isfinite(ratio[0]) or np.isnan(ratio[0])),'Solved native boundary balance unavailable; only unsolved initial ratio may be NaN')
        b['initial_boundary_ratio']=dict(available=bool(np.isfinite(ratio[0])),value=float(ratio[0]) if np.isfinite(ratio[0]) else None,solved_gate_applied=False)
        require(np.array_equal(b['solved_time_s'],b['time_s'][1:]),'Unsolved initialization must not acquire residual/force rows')
        for term in ('pressure','viscous','total'):require(b['body_'+term+'_force_coefficient_history'].shape==(n-1,2,2),'Actual solved per-body force shape differs')
        require(z['body_names'].astype(str).tolist()==['body1','body2'] and z['force_component_names'].astype(str).tolist()==['x','y'],'Native body/component ordering differs')
        if name=='fluent':
            b['terminal_scaled_residuals']=z['terminal_scaled_residuals'];b['terminal_residual_equations']=z['terminal_residual_equations'].astype(str).tolist()
            require(b['terminal_residual_equations']==['continuity','x-velocity','y-velocity'] and b['terminal_scaled_residuals'].shape==(n-1,3) and np.isfinite(b['terminal_scaled_residuals']).all(),'Actual Fluent own residual history differs')
        else:
            b['residual_component_names']=z['residual_component_names'].astype(str).tolist();require(b['residual_component_names'] in (['p','Ux','Uy'],['p','Ux','Uy','Uz']),'Actual printed native OpenFOAM components required')
            for key in ('last_outer_initial_residuals','last_outer_final_linear_residuals'):
                b[key]=z[key];require(b[key].shape==(n-1,len(b['residual_component_names'])) and np.isfinite(b[key]).all(),'Actual OpenFOAM residual shape differs')
            b['outer_corrector_count']=z['outer_corrector_count'];require(b['outer_corrector_count'].shape==(n-1,) and np.isfinite(b['outer_corrector_count']).all(),'Actual outer count history differs')
    b.update(name=name,sources=sources,speed=np.linalg.norm(b['U_source_history_m_s'],axis=2))
    b['speed_saturation_cell_counts']=np.count_nonzero(b['speed']>.2,axis=1)
    return b


def verify_fluent_frames(b,reader_path,g,indices):
    sources=b['sources'];path=sources.consume(reader_path);require(sha(path)==FLUENT_READER_SHA,'Selected native Fluent reader differs from frozen reviewed source')
    spec=importlib.util.spec_from_file_location('_tandem_native_fluent_reader',path);reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    rows=sources.record['native_state_paths'];require(len(rows)==len(b['time_s']),'Literal native Fluent state rows required');proof=[];mesh=None
    for step in indices:
        common.resource_guard();row=rows[int(step)];require(row['index']==step and abs(row['time_s']-b['time_s'][step])<1e-10 and row['source_kind']=='native_fluent_hdf','Actual native Fluent source row differs')
        case=sources.consume(row['case']);data=sources.consume(row['data'])
        if mesh is None:mesh=reader.map_native_mesh(case,g)
        s=reader.read_state(case,data,g,mesh=mesh)
        require(np.array_equal(np.argsort(s['mapping']),b['source_native_cell_index']),'Actual Fluent native/source topology mapping differs')
        U=np.column_stack((s['common']['u'],s['common']['v']));require(np.array_equal(U,b['U_source_history_m_s'][step]) and np.array_equal(s['common']['p'],b['p_source_history_Pa'][step]),'Raw selected Fluent U/p differs from audited arrays')
        t=s['native_time'];require(t['time-step']==step and abs(t['flow-time']-b['time_s'][step])<1e-10 and abs(t['physical-time-step']-.05)<1e-12,'Raw Fluent time/step/dt differs')
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
                path=sources.consume(row['rank_time_metadata'][rank]);raw=common.native.foam_text(path);values={key:float(re.search(r'\b'+key+r'\s+([-+0-9.eE]+)\s*;',raw)[1]) for key in ('value','index','deltaT')};require(abs(values['value']-step*.05)<1e-10 and values['index']==step and abs(values['deltaT']-.05)<1e-12,'Actual native rank clock/dt/index differs')
        order=b['source_native_cell_index'];require(np.array_equal(U[order,:2],b['U_source_history_m_s'][step]) and np.array_equal(p[order],b['p_source_history_Pa'][step]),'Actual selected two-rank U/p differs from audited arrays')
        proof.append(dict(step=int(step),rank_U=row['rank_U'],rank_p=row['rank_p'],complete_cell_vertex_bijection=True,raw_U_p_exact=True))
    b['native_readback']=dict(source_kind='native_rank_reassembly',selected_states=proof,source_quad_cells=N,vertex_set_bijection=True,selected_velocity_values_verified=len(indices)*N*2,selected_pressure_values_verified=len(indices)*N)


def notes(bundles):
    result=[]
    for name,b in bundles.items():
        s=b['sources'].status;failed=s['failed_checks']+s['failed_state_checks'];result.extend(textwrap.wrap(NAMES[name]+(': primary numerical checks passed' if not failed else ': primary numerical MISSED '+', '.join(failed)),width=155))
    return result


def map_canvas(g,subtitle,bundles=None,state=0,geometry_only=False):
    fig=plt.figure(figsize=(10.24,10.4),dpi=100)
    fig.text(.08,.972,'Offset tandem-square wake',fontsize=20,fontweight='bold');fig.text(.08,.939,subtitle,fontsize=9.5,color='#425e6d');clock=fig.text(.08,.911,'Geometry only | no solver field loaded' if geometry_only else '',fontsize=9,color='#304b59')
    artists={};axes={};width=.87;h=width*10.24*.75/2/10.4
    for bottom,name in ((.55,'fluent'),(.15,'openfoam')):
        ax=fig.add_axes([.08,bottom,width,h]);kw=dict(facecolors='#f3f6f8',edgecolors='#adbdc5',linewidths=.15) if geometry_only else dict(array=bundles[name]['speed'][state],cmap='viridis',norm=Normalize(0,.2),edgecolors='none',antialiased=False)
        col=PolyCollection(g['polygons'],**kw);ax.add_collection(col)
        for zone,centre,label in ((8,(.5,-.025),'1'),(9,(.9,.025),'2')):
            ax.add_patch(Polygon(g['body_outlines'][zone],facecolor='#f0f2f3',edgecolor='#263f4c',lw=.9,zorder=3));ax.text(*centre,label,ha='center',va='center',fontsize=7.5,color='#263f4c',zorder=4)
        ax.set(xlim=(0,2),ylim=(-.375,.375),xlabel='x (m)',ylabel='y (m)');ax.set_aspect('equal');ax.set_xticks([0,.25,.5,.75,1,1.25,1.5,1.75,2]);ax.set_yticks([-.375,-.25,0,.25,.375]);ax.tick_params(labelsize=8);ax.set_title(NAMES[name]+(' | planned native-field layout' if geometry_only else ' | actual native speed'),loc='left',fontsize=10.5,pad=8,fontweight='bold');artists[name]=col;axes[name]=ax
    if geometry_only:fig.text(.08,.04,'Exact 9,472-cell source mesh; both solid holes remain empty of flow data. Planned speed scale: 0-0.2 m/s.',fontsize=8,color='#425e6d')
    else:
        cax=fig.add_axes([.25,.077,.54,.012]);fig.colorbar(ScalarMappable(norm=Normalize(0,.2),cmap='viridis'),cax=cax,orientation='horizontal',label='Speed (m/s) | fixed common scale');common.add_footer(fig,notes(bundles),bottom=.014)
    return fig,artists,clock,axes


def gif_palette():
    rgb=np.rint(matplotlib.colormaps['viridis'](np.linspace(0,1,128))[:,:3]*255).astype('uint8').tolist()
    rgb += [[v,v,v] for v in range(0,256,4)]
    rgb += [[255,255,255],[38,63,76],[66,94,109],[48,75,89],[240,242,243],[0,132,152],[214,93,40]]
    rgb += [[255,255,255]]*(256-len(rgb));palette=Image.new('P',(1,1));palette.putpalette(np.asarray(rgb,dtype=np.uint8).ravel().tolist());return palette


def render_animation(g,bundles,out,phase):
    indices=np.arange(6) if phase=='pilot' else np.arange(0,801,8);subtitle='Pilot | five actual 0.05 s steps; no full-run completion claim' if phase=='pilot' else 'ReD = 100 | two stationary squares | 800 actual 0.05 s steps'
    fig,artists,clock,axes=map_canvas(g,subtitle,bundles);images=[];crops={};palette=gif_palette()
    for j,step in enumerate(indices):
        if j%10==0:common.resource_guard()
        ts=[b['time_s'][step] for b in bundles.values()];require(max(ts)-min(ts)<1e-10,'Actual paired native frame times differ')
        for name,b in bundles.items():artists[name].set_array(b['speed'][step])
        counts=[int(b['speed_saturation_cell_counts'][step]) for b in bundles.values()];saturation=f' | >0.2 m/s cells F/OF: {counts[0]}/{counts[1]}' if any(counts) else ''
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
    proof=dict(frames=len(indices),frame_delay_ms=60,actual_native_step_indices=indices.tolist(),actual_physical_times_s=bundles['fluent']['time_s'][indices].tolist(),actual_solver_velocity_evolution=True,temporal_interpolation=False,particles=False,streamlines=False,signed_vorticity_used=False,common_speed_limits_m_s=[0,.2],speed_above_scale_cell_counts={n:b['speed_saturation_cell_counts'][indices].tolist() for n,b in bundles.items()},true_equal_xy_aspect=True,exact_two_body_holes=True,distinct_decoded_RGB_frames=len(set(full)),map_region_unique_frames={k:len(set(v)) for k,v in maps.items()},map_crop_pixel_bounds=crops,map_frame_uniqueness_requirement=False,plateau_policy='Actual settled or quantized map repeats are retained',palette=dict(fixed=True,speed_colour_levels=128,dither=False),pixel_dimensions=list(images[0].size),bytes=path.stat().st_size,maximum_public_bytes=15*1024**2,sha256=sha(path))
    print(json.dumps(dict(actual_animation_ready=str(path),animation=proof)),flush=True);return proof


def decorate_axes(axes):
    for ax in np.asarray(axes).ravel():ax.grid(alpha=.17);ax.spines[['top','right']].set_visible(False);ax.legend(frameon=True,facecolor='white',framealpha=1,edgecolor='none',fontsize=7.2)


def render_profiles(g,bundles,out):
    fig,axes=plt.subplots(1,3,figsize=(12.8,5.8));fig.subplots_adjust(left=.08,right=.965,top=.79,bottom=.24,wspace=.34);fig.suptitle('Final native wake profiles',x=.08,y=.96,ha='left',fontsize=20,fontweight='bold');fig.text(.08,.88,'Actual cell-centre columns downstream of the squares; no interpolated or appended wall samples.',fontsize=9.5,color='#425e6d')
    for ax,ix in zip(axes,PROFILE_COLUMNS):
        cells=np.flatnonzero(g['ix']==ix);cells=cells[np.argsort(g['iy'][cells])]
        for name,b in bundles.items():
            ax.plot(b['U_source_history_m_s'][-1,cells,0]/.1,b['source_centres_m'][cells,1],color=COLOURS[name],ls='--' if name=='openfoam' else '-',label=NAMES[name])
        ax.axvline(1,color='#71828b',lw=.7,ls=':');ax.set(title=f'x = {g["centres"][cells[0],0]:.5f} m | native ix {ix}',xlabel='Streamwise velocity Ux / Uin',ylabel='Actual y (m)',ylim=(-.375,.375))
    decorate_axes(axes);fig.text(.08,.12,'Uin = 0.1 m/s. These are final saved profiles; no steady, periodic or benchmark agreement gate is imposed.',fontsize=8.5,color='#425e6d');common.add_footer(fig,notes(bundles));common.save(fig,out,PREFIX+'-final-wake-profiles')


def render_forces(bundles,out):
    fig,axes=plt.subplots(2,2,figsize=(12.8,9.4));fig.subplots_adjust(left=.09,right=.965,top=.85,bottom=.18,wspace=.28,hspace=.42);fig.suptitle('Native force coefficients on each square',x=.09,y=.96,ha='left',fontsize=20,fontweight='bold');fig.text(.09,.9,'Pressure + viscous forces; each code uses its own actual depth and the common per-body coefficient reference.',fontsize=9,color='#425e6d')
    for body in range(2):
        for component,label in enumerate(('Drag coefficient Cd','Lift coefficient Cl')):
            ax=axes[body,component]
            for name,b in bundles.items():ax.plot(b['solved_time_s'],b['body_total_force_coefficient_history'][:,body,component],ls='--' if name=='openfoam' else '-',color=COLOURS[name],lw=1.2,label=NAMES[name])
            ax.set(title=f'Body {body+1} | total {"x" if component==0 else "y"} force',xlabel='Solved physical time (s)',ylabel=label)
    decorate_axes(axes);fig.text(.09,.105,'Every solved startup sample is retained. Initial forces are unavailable. No smoothing, plateau, zero-lift or periodicity gate.',fontsize=8.5,color='#425e6d');common.add_footer(fig,notes(bundles));common.save(fig,out,PREFIX+'-per-body-force-histories')


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
            for body in range(2):
                for j,component in enumerate(('Cd','Cl')):parts.append(b['body_'+term+'_force_coefficient_history'][:,body,j]);labels.append(f'body{body+1}_{term}_{component}')
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
    with (out/'record-final-wake-profiles.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','native_step','physical_time_s','source_cell_0based','native_ix','native_iy','x_m','y_m','u_m_s','v_m_s','Ux_over_Uin'])
        for name,b in bundles.items():
            for ix in PROFILE_COLUMNS:
                for cell in np.flatnonzero(g['ix']==ix):w.writerow([name,int(b['step_index'][-1]),b['time_s'][-1],cell,ix,g['iy'][cell],*b['source_centres_m'][cell],*b['U_source_history_m_s'][-1,cell],b['U_source_history_m_s'][-1,cell,0]/.1])
    return dict(full_external_field_rows=sum(len(b['time_s'])*N for b in bundles.values()),record_final_field_rows=2*N,record_solved_history_rows=len(bundles['fluent']['solved_time_s']),record_final_profile_rows=2*3*60,full_external_field_file='actual-native-cell-fields.csv.gz',compact_record_files=['record-final-native-fields.csv','record-solved-step-history.csv','record-final-wake-profiles.csv'])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--phase',choices=('geometry','pilot','full'),required=True);ap.add_argument('--geometry',type=Path,required=True);ap.add_argument('--targets',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--source-sha256',required=True);ap.add_argument('--fluent-reader',type=Path)
    for name in NAMES:
        for suffix in ('audit','manifest','fields'):ap.add_argument('--'+name+'-'+suffix,type=Path)
        ap.add_argument('--'+name+'-fields-sha256')
    a=ap.parse_args();out=a.out.resolve();require(sha(__file__)==a.source_sha256,'Executed plotter differs from reviewed source hash');require(out.drive.upper()=='E:' and not out.exists(),'Fresh E output directory required');require(sha(a.targets)==TARGETS_SHA,'Frozen model targets changed')
    paths=[Path(__file__),a.geometry,a.targets]
    for filename,digest in HELPERS.items():p=Path(__file__).with_name(filename);require(sha(p)==digest,'Frozen reused helper changed');paths.append(p)
    base_inputs={str(p.resolve()):sha(p) for p in paths};g=geometry(a.geometry);resources=common.resource_guard(render_reserve=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.5,'axes.linewidth':.7})
    if a.phase=='geometry':
        out.mkdir(parents=True);fig,_,_,_=map_canvas(g,'Geometry-only preparation | exact source mesh and planned map layout',geometry_only=True);common.save(fig,out,PREFIX+'-geometry-only-layout');(out/'geometry-only-proof.json').write_text(json.dumps(dict(solver_fields_loaded=False,native_solver_launched=False,actual_cells=N,actual_nodes=9723,holes=2,actual_fluid_area_m2=g['actual_polygon_area_m2'],profile_columns=list(PROFILE_COLUMNS),all_profile_columns_have60_actual_cells=True,geometry_sha256=sha(a.geometry),source_sha256=sha(__file__),resources=resources),indent=2)+'\n',encoding='utf-8');return
    bundles={}
    for name in NAMES:
        audit,manifest,fields=[getattr(a,name+'_'+suffix) for suffix in ('audit','manifest','fields')];digest=getattr(a,name+'_fields_sha256');require(all(p is not None for p in (audit,manifest,fields)) and digest is not None,'Explicit accepted audit/manifest/field paths and field hash required');source=AcceptedSources(name,audit,manifest,fields,digest);bundles[name]=load_history(name,fields,source,g,a.phase)
    require(a.fluent_reader is not None,'Explicit audit-bound native Fluent reader required');indices=np.arange(6) if a.phase=='pilot' else np.arange(0,801,8)
    from ansys.fluent.core.offline import _lispy
    base_inputs[str(Path(_lispy.__file__).resolve())]=sha(_lispy.__file__);verify_fluent_frames(bundles['fluent'],a.fluent_reader,g,indices);verify_openfoam_frames(bundles['openfoam'],g,indices)
    out.mkdir(parents=True);animation=render_animation(g,bundles,out,a.phase);render_profiles(g,bundles,out);render_forces(bundles,out);render_numerics(bundles,out);print(json.dumps(dict(science_PNGs_ready=str(out))),flush=True);csv_proof=export_csvs(g,bundles,out)
    inputs=dict(base_inputs)
    for b in bundles.values():inputs.update(b['sources'].used)
    require(all(sha(Path(p))==h for p,h in inputs.items()),'Consumed native/source input changed during export')
    checks=dict(phase=a.phase,actual_native_fields=True,numerical_claim_scope='Primary native numerical audit checks; artifact diagnostics remain separate',numerical_audits={n:b['sources'].status for n,b in bundles.items()},native_mapping={n:b['native_readback'] for n,b in bundles.items()},animation=animation,input_hashes_stable=True,accepted_full_audit_not_repeated=True,independently_consumed_source_count=len(inputs),common_cell_count=N,actual_polygon_area_m2=g['actual_polygon_area_m2'],exact_body_holes=2,initial_boundary_ratio={n:b['initial_boundary_ratio'] for n,b in bundles.items()},initial_solved_conservation_claim=False,initial_forces_available=False,solved_only_force_history=True,raw_native_gauge_pressure_Pa=True,pressure_mean_subtraction=False,speed_scale_saturation_counts_all_native_states={n:b['speed_saturation_cell_counts'].tolist() for n,b in bundles.items()},profile_native_columns=list(PROFILE_COLUMNS),profile_native_x_m=[float(g['x_edges'][ix]+.00625) for ix in PROFILE_COLUMNS],csv=csv_proof,openfoam_diagnostic_comparisons=bundles['openfoam']['sources'].record.get('diagnostic_comparisons',{}),openfoam_execution_completion=bundles['openfoam']['sources'].record.get('execution_completion',{}),physical_validation=False,mesh_independence=False,time_step_independence=False,periodicity_established=False,steady_force_plateau_claim=False,mean_lift_zero_gate=False,signed_vorticity_used=False,resources=dict(before=resources,after=common.resource_guard()))
    (out/'checks.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n',encoding='utf-8');(out/'input-hash-manifest.json').write_text(json.dumps(dict(inputs=[dict(path=p,sha256=h) for p,h in sorted(inputs.items())]),indent=2)+'\n',encoding='utf-8');(out/'executed-plotter.py').write_bytes(Path(__file__).read_bytes())
    caption=f"Actual native offset tandem-square velocity fields, two stationary 0.1 m squares, Re_D100. This {a.phase} export selects {len(indices)} actual saved states, steps {indices[0]} to {indices[-1]}, physical times {bundles['fluent']['time_s'][indices[0]]:.2f} to {bundles['fluent']['time_s'][indices[-1]]:.2f} s, at 60 ms per playback frame. The exact 9472-cell, 1.48 m2 fluid mesh has two solid holes. Maps use equal xy aspect and a fixed common 0 to 0.2 m/s speed scale; any saturation is labeled and actual counts retained. No temporal interpolation, particles, synthesized vorticity or forced map uniqueness. Genuine settled/quantized repeats are preserved. Fixed palette quantization affects preview colors only; static maps and CSV retain native values.\n\nInitialization is prescribed cell rest with actual native boundary flux, not a solved conservation state. Undefined initial relative boundary ratio remains unavailable, never replaced with zero. Mass/divergence histories and residual targets begin at solved step1; forces are solved-only with no invented initial row. Per-body pressure, viscous and total coefficients use each code's native depth and the same reference D/U/rho; different-depth raw forces are not compared. Native residual norms remain separate. Courant number is descriptive. Profiles use actual columns ix52/84/112 and native cell centres without appended wall values. Native gauge pressure is exported in Pa without mean subtraction.\n\nNo steady-force plateau, zero mean lift, periodicity, benchmark agreement, physical validation, mesh independence or time-step independence is claimed. Full native field CSV remains external; compact final maps and all solved-step history rows are provided for the record.\n\n"+'\n'.join(notes(bundles))+'\n'
    for name,b in bundles.items():
        record=b['sources'].record
        for key in ('diagnostic_comparisons','execution_completion'):
            if record.get(key):caption+='\n'+NAMES[name]+' exact native '+key+' (separate from primary numerical checks):\n```json\n'+json.dumps(record[key],indent=2,allow_nan=False)+'\n```\n'
    (out/'CAPTIONS.md').write_text(caption,encoding='utf-8');files=sorted(p for p in out.iterdir() if p.is_file());(out/'output-hash-manifest.json').write_text(json.dumps([dict(path=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in files],indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(output=str(out),animation=animation,csv=csv_proof),indent=2))


if __name__=='__main__':main()
