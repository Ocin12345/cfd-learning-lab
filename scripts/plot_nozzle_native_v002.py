"""Actual native nozzle fields and fixed-snapshot tracers; never launch CFD.
Native vertex sets map bijectively to common source quads. No imposed shock,
synthetic velocity, appended boundary samples or solver-time animation.
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, math, os, re, sys
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-nozzle'
os.environ['TEMP']=os.environ['TMP']='E:/CFD-Learning-Lab/temp'
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import Normalize
from PIL import Image, ImageDraw
NAMES={'fluent':'ANSYS Fluent','openfoam':'OpenFOAM'}
TOL=1e-8
GAMMA,R_GAS,P_TOTAL,T_TOTAL,P_BACK=1.4,287.05,200000.,300.,140000.

def require(ok,message):
    if not ok: raise ValueError(message)

def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''): result.update(block)
    return result.hexdigest()

def indexed(group,count,width=1,partial=False):
    values=np.full((count,width),np.nan); seen=np.zeros(count,bool)
    for dataset in group.values():
        require(isinstance(dataset,h5py.Dataset),'Unsupported nested native HDF block')
        first=int(np.asarray(dataset.attrs['minId']).ravel()[0])-1
        last=int(np.asarray(dataset.attrs['maxId']).ravel()[0])
        block=np.asarray(dataset[:]).reshape(-1,width)
        require(0<=first<last<=count and len(block)==last-first and not seen[first:last].any(),'Invalid native HDF IDs')
        values[first:last],seen[first:last]=block,True
    require(partial or seen.all(),'Incomplete native HDF ID coverage')
    return values[:,0] if width==1 else values

def foam_text(path):
    return re.sub(r'/\*.*?\*/|//[^\n]*','',Path(path).read_text(encoding='utf-8'),flags=re.S)

def foam_list(path):
    raw=foam_text(path); header=re.search(r'\bFoamFile\s*\{([^}]*)\}',raw,re.S)
    require(header is not None and re.search(r'\bformat\s+ascii\s*;',header[1]),'Saved ASCII list required: '+str(path))
    match=re.fullmatch(r'\s*(\d+)\s*\((.*)\)\s*;?\s*',raw[header.end():],re.S)
    require(match is not None,'Unsupported native list: '+str(path))
    return int(match[1]),match[2]

def foam_field(path,count,width,dimensions):
    raw=foam_text(path)
    require(re.search(r'\bformat\s+ascii\s*;',raw),'Saved ASCII field required')
    units=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',raw)
    require(units is not None and np.array_equal(np.fromstring(units[1],sep=' '),dimensions),'Unexpected field dimensions: '+str(path))
    kind='vector' if width==3 else 'scalar'
    match=re.search(r'\binternalField\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',raw,re.S)
    if match:
        require(int(match[1])==count,'Native field cell count differs')
        values=np.fromstring(match[2].replace('(',' ').replace(')',' '),sep=' ')
    else:
        match=re.search(r'\binternalField\s+uniform\s+([^;]+);',raw)
        require(match is not None,'No supported saved internalField')
        value=np.fromstring(match[1].replace('(',' ').replace(')',' '),sep=' ')
        require(value.size==width,'Invalid uniform field')
        values=np.tile(value,count)
    require(values.size==count*width and np.isfinite(values).all(),'Nonfinite or wrong-sized native field')
    return values.reshape(count,width) if width>1 else values

class QuadGeometry:
    def __init__(self,path):
        self.path=Path(path)
        with np.load(path,allow_pickle=False) as s:
            self.nodes=s['nodes'].copy(); self.quads=s['quads'].astype(np.int64)
            self.centres=s['centres'].copy(); self.areas=s['areas'].copy()
            self.ny,self.nx=map(int,s['grid_shape']); self.ix=s['ix'].astype(int); self.iy=s['iy'].astype(int)
            self.x_edges=s['x_edges'].copy(); self.h_edges=s['half_height_edges'].copy()
        n,count=len(self.nodes),len(self.quads)
        require(self.nodes.shape==(n,2) and self.quads.shape==(count,4) and count==self.nx*self.ny and self.centres.shape==(count,2),'Invalid source geometry shape')
        require(np.isfinite(self.nodes).all() and np.isfinite(self.centres).all() and np.isfinite(self.areas).all() and self.areas.min()>0 and np.all((self.quads>=0)&(self.quads<n)),'Invalid source geometry values')
        self.polygons=self.nodes[self.quads]; x,y=self.polygons[:,:,0],self.polygons[:,:,1]
        cross=x*np.roll(y,-1,axis=1)-np.roll(x,-1,axis=1)*y; signed=cross.sum(axis=1)/2
        centres=np.column_stack((((x+np.roll(x,-1,axis=1))*cross).sum(axis=1),((y+np.roll(y,-1,axis=1))*cross).sum(axis=1)))/(6*signed[:,None])
        require(np.all(np.abs(signed)>0) and np.max(np.abs(np.abs(signed)-self.areas))<=1e-12 and np.max(np.abs(centres-self.centres))<=TOL,'Source polygon centroids or areas disagree')
        slots=self.iy*self.nx+self.ix
        require(np.array_equal(np.sort(slots),np.arange(count)),'Source cell indices not bijective')
        self.order=np.argsort(slots); c=self.centres[self.order].reshape(self.ny,self.nx,2)
        require(np.ptp(c[:,:,0],axis=0).max()<TOL,'Expected common centroid x per cell column')
        self.xc=c[0,:,0]
        require(np.all(np.diff(self.x_edges)>0) and np.all(np.diff(self.xc)>0) and self.x_edges.shape==(self.nx+1,) and self.h_edges.shape==(self.nx+1,),'Invalid wall stations')
        self.eta=c[:,:,1]/self.height(self.xc)[None,:]
        require(np.all(np.diff(self.eta,axis=0)>0),'Transverse cell rows not ordered')
        self.node_lookup={tuple(np.rint(p/TOL).astype(np.int64)):i for i,p in enumerate(self.nodes)}
        self.cell_lookup={tuple(row):i for i,row in enumerate(np.sort(self.quads,axis=1))}
        require(len(self.node_lookup)==n and len(self.cell_lookup)==count,'Repeated or quantization-colliding source entities')
    def height(self,x): return np.interp(x,self.x_edges,self.h_edges)
    def native_node_map(self,points):
        require(points.ndim==2 and points.shape[1] in (2,3) and np.isfinite(points).all(),'Invalid native nodes')
        mapped=np.empty(len(points),np.int64); maximum=0.
        for i,p in enumerate(points[:,:2]):
            node=self.node_lookup.get(tuple(np.rint(p/TOL).astype(np.int64)))
            if node is None or np.linalg.norm(self.nodes[node]-p)>TOL:
                distance=np.linalg.norm(self.nodes-p,axis=1); node=int(np.argmin(distance))
                require(distance[node]<=TOL,'Native node absent from common source geometry')
            mapped[i]=node; maximum=max(maximum,float(np.linalg.norm(self.nodes[node]-p)))
        return mapped,maximum
    def match(self,native):
        nodes,error=self.native_node_map(native['nodes']); ids=np.empty(len(native['cells']),np.int64)
        for i,vertices in enumerate(native['cells']):
            key=tuple(np.unique(nodes[vertices]))
            require(len(key)==4 and key in self.cell_lookup,'Native cell does not project to a source quad')
            ids[i]=self.cell_lookup[key]
        require(np.array_equal(np.sort(ids),np.arange(len(self.quads))),'Native cell map not bijective')
        order=np.argsort(ids); result={key:native[key][order] for key in ('U','p','T')}
        if 'rho' in native:
            result['rho']=native['rho'][order]; eos=result['p']/(result['rho']*result['T'])
            require(np.allclose(eos,R_GAS,rtol=3e-4,atol=.03),'Native p/(rho T) differs from R287.05')
            result['maximum_EOS_R_relative_error']=float(np.max(np.abs(eos/R_GAS-1)))
        require(all(np.isfinite(result[k]).all() for k in ('U','p','T')) and result['p'].min()>0 and result['T'].min()>0,'Nonphysical native absolute pressure or temperature')
        result['speed']=np.linalg.norm(result['U'],axis=1); result['Mach']=result['speed']/np.sqrt(GAMMA*R_GAS*result['T'])
        result['gridU']=result['U'][self.order].reshape(self.ny,self.nx,2)
        for key in ('Mach','p','T'): result['grid'+key]=result[key][self.order].reshape(self.ny,self.nx)
        result['mapping']=dict(native_cell_to_common_quad_bijective=True,vertex_set_mapping=True,maximum_native_node_error_m=error,projected_single_layer=native['layer'])
        if 'C' in native:
            e=float(np.max(np.linalg.norm(native['C'][order,:2]-self.centres,axis=1)))
            require(e<=TOL,'Saved C differs from source polygon centroids')
            result['mapping']['maximum_saved_C_centroid_error_m']=e
        result['solve_evidence']=native['solve_evidence']; return result
    def interpolate(self,point,grid):
        x,y=point
        if not self.xc[0]<=x<=self.xc[-1] or abs(y)>=self.height(x): return None
        eta=y/self.height(x); col=max(0,min(int(np.searchsorted(self.xc,x,side='right'))-1,self.nx-2))
        if not max(self.eta[0,col],self.eta[0,col+1])<=eta<=min(self.eta[-1,col],self.eta[-1,col+1]): return None
        values=[]
        for j in (col,col+1):
            row=max(0,min(int(np.searchsorted(self.eta[:,j],eta,side='right'))-1,self.ny-2))
            f=(eta-self.eta[row,j])/(self.eta[row+1,j]-self.eta[row,j]); values.append((1-f)*grid[row,j]+f*grid[row+1,j])
        f=(x-self.xc[col])/(self.xc[col+1]-self.xc[col]); return (1-f)*values[0]+f*values[1]
    def trace(self,seed,grid):
        p=np.array(seed,float); time=0.; points=[p.copy()]; times=[0.]; speeds=[]; stop='maximum integration steps'
        for _ in range(6000):
            k1=self.interpolate(p,grid)
            if k1 is None: stop='outside cell-centre interpolation envelope'; break
            speed=float(np.linalg.norm(k1))
            if not speeds: speeds.append(speed)
            if speed<1e-8: stop='stagnation below1e-8m/s'; break
            dt=min(float(np.min(np.diff(self.x_edges)))/(3*speed),1e-5)
            k2=self.interpolate(p+dt*k1/2,grid); k3=None if k2 is None else self.interpolate(p+dt*k2/2,grid)
            k4=None if k3 is None else self.interpolate(p+dt*k3,grid)
            if k4 is None: stop='RK4 stage reached interpolation envelope'; break
            next_p=p+dt*(k1+2*k2+2*k3+k4)/6; value=self.interpolate(next_p,grid); segment=np.linspace(p,next_p,5)
            if value is None or np.any(np.abs(segment[:,1])>=self.height(segment[:,0])): stop='segment reached envelope or actual wall'; break
            time+=dt; p=next_p; points.append(p.copy()); times.append(time); speeds.append(float(np.linalg.norm(value)))
        require(len(points)>=8 and len(points)==len(times)==len(speeds),'Unusable native tracer trajectory')
        return dict(seed=list(seed),points=np.array(points),times=np.array(times),speeds=np.array(speeds),stop=stop)

def fluent_native(case,data,operating_pressure):
    with h5py.File(case,'r') as handle:
        mesh=handle['meshes/1']; top=mesh['cells/zoneTopology']
        require(np.asarray(top['minId']).tolist()==[1] and len(np.asarray(top['maxId']))==1,'Expected one native fluid zone')
        count=int(top['maxId'][0]); groups=mesh['nodes/coords']
        total=max(int(np.asarray(d.attrs['maxId']).ravel()[0]) for d in groups.values()); width=next(iter(groups.values())).shape[1]
        nodes=indexed(groups,total,width)
        require(width in (2,3) and (width==2 or np.ptp(nodes[:,2])<=TOL),'Nonplanar Fluent mesh')
        blocks=mesh['faces/nodes']; nf=max(int(np.asarray(b['nnodes'].attrs.get('maxId',b.attrs.get('maxId'))).ravel()[0]) for b in blocks.values())
        faces=np.full((nf,2),-1,np.int64)
        for block in blocks.values():
            sizes=block['nnodes']; start=int(np.asarray(sizes.attrs.get('minId',block.attrs.get('minId',[1]))).ravel()[0])-1
            require(np.all(sizes[:]==2),'Expected native two-node Fluent faces')
            values=np.asarray(block['nodes'][:],np.int64).reshape(-1,2)-1
            require(np.all(faces[start:start+len(values)]==-1),'Repeated native face IDs'); faces[start:start+len(values)]=values
        require(np.all((faces>=0)&(faces<total)),'Invalid native face nodes')
        owner=indexed(mesh['faces/c0'],nf).astype(np.int64)-1; other=indexed(mesh['faces/c1'],nf,partial=True)
        valid=np.isfinite(other); neighbour=other[valid].astype(np.int64)-1
        require(np.all((owner>=0)&(owner<count)) and np.all((neighbour>=0)&(neighbour<count)),'Invalid native cell IDs')
        keys=np.unique(np.concatenate((np.repeat(owner,2)*total+faces.ravel(),np.repeat(neighbour,2)*total+faces[valid].ravel())))
        ids,sizes=np.unique(keys//total,return_counts=True)
        require(np.array_equal(ids,np.arange(count)) and np.all(sizes==4),'Native Fluent cells require four vertices')
        cells=(keys%total).reshape(count,4)
    with h5py.File(data,'r') as handle:
        phase=handle['results/1/phase-1/cells']; U=np.column_stack([indexed(phase[k],count) for k in ('SV_U','SV_V')])
        pressure=indexed(phase['SV_P'],count)+operating_pressure; temp=indexed(phase['SV_T'],count)
        rho=indexed(phase['SV_DENSITY'],count) if 'SV_DENSITY' in phase else None
        require('results/residuals/phase-1' in handle,'Initialized data without actual Fluent residual history refused')
        residuals={}; terminals=[]
        for name,group in handle['results/residuals/phase-1'].items():
            if 'iterations' not in group or 'data' not in group or not len(group['iterations']): continue
            iteration=float(group['iterations'][-1]); row=np.asarray(group['data'][-1]).ravel()
            require(np.isfinite(iteration) and iteration>0 and iteration==int(iteration) and len(row)>=2 and np.isfinite(row[:2]).all() and row[1]>0,'Invalid actual Fluent residual evidence')
            residuals[name]=float(row[0]/row[1]); terminals.append(int(iteration))
        require(len(residuals)>=3 and len(set(terminals))==1,'Incomplete actual Fluent solve evidence')
    result=dict(nodes=nodes,cells=cells,U=U,p=pressure,T=temp,layer='native planar2D quads',solve_evidence=dict(actual_terminal_iteration=terminals[0],terminal_scaled_residuals=residuals,operating_pressure_Pa=operating_pressure,initialized_only_data_refused=True))
    if rho is not None: result['rho']=rho
    return result

def foam_mesh(case):
    mesh=next((p for p in (case/'constant/polyMesh',case/'0/polyMesh') if (p/'points').is_file()),None)
    require(mesh is not None,'No native OpenFOAM polyMesh')
    files=[mesh/name for name in ('points','faces','owner','neighbour')]; n,body=foam_list(files[0])
    nodes=np.fromstring(body.replace('(',' ').replace(')',' '),sep=' ').reshape(-1,3)
    require(len(nodes)==n and np.isfinite(nodes).all() and len(np.unique(nodes[:,2]))==2,'Expected one native extruded cell layer')
    nf,body=foam_list(files[1]); matches=list(re.finditer(r'(\d+)\s*\(([^()]*)\)',body))
    require(len(matches)==nf and all(int(m[1])==4 for m in matches),'Expected native hexahedron quad faces')
    faces=np.array([np.fromstring(m[2],sep=' ',dtype=np.int64) for m in matches]); labels=[]
    for path in files[2:]:
        n,body=foam_list(path); values=np.fromstring(body,sep=' ',dtype=np.int64)
        require(len(values)==n,'Native label count differs'); labels.append(values)
    owner,neighbour=labels; count=int(max(owner.max(),neighbour.max(initial=-1)))+1
    require(len(owner)==nf and len(neighbour)<=nf and np.all((faces>=0)&(faces<len(nodes))) and np.all((owner>=0)&(owner<count)) and np.all((neighbour>=0)&(neighbour<count)),'Invalid native adjacency')
    require(np.all(np.bincount(owner,minlength=count)+np.bincount(neighbour,minlength=count)==6),'Native hexahedra require six faces')
    keys=np.unique(np.concatenate((np.repeat(owner,4)*len(nodes)+faces.ravel(),np.repeat(neighbour,4)*len(nodes)+faces[:len(neighbour)].ravel())))
    ids,sizes=np.unique(keys//len(nodes),return_counts=True)
    require(np.array_equal(ids,np.arange(count)) and np.all(sizes==8),'Native hexahedra require eight vertices')
    return dict(nodes=nodes,cells=(keys%len(nodes)).reshape(count,8),inputs=files)

def foam_saved_folder(case,selected):
    if selected is not None:
        require(float(selected)>0,'Prepared OpenFOAM timezero refused'); return case/selected
    candidates=[]
    for path in case.iterdir():
        try: value=float(path.name)
        except ValueError: continue
        if path.is_dir() and np.isfinite(value) and value>0 and all((path/name).is_file() for name in ('U','p','T')): candidates.append((value,path))
    require(bool(candidates),'No positive computed OpenFOAM snapshot')
    return max(candidates,key=lambda item:item[0])[1]

def foam_solve_evidence(log,folder):
    raw=Path(log).read_text(encoding='utf-8',errors='replace'); value=float(folder.name)
    events=list(re.finditer(r'^\s*Time\s*=\s*([-+0-9.eE]+)\s*s?\s*$',raw,re.M))
    for index,event in enumerate(events):
        if abs(float(event[1])-value)>1e-8: continue
        block=raw[event.start():events[index+1].start() if index+1<len(events) else len(raw)]
        if re.search(r'Solving for (?:U(?:x|y|z)?|rhoU(?:x|y|z)?)\s*,',block) and re.search(r'Solving for (?:p|rho|rhoE|e|h)\s*,',block):
            return dict(saved_pseudo_iteration=value,solver_time_is_physical=False,native_solver_log=str(log),initialized_only_data_refused=True,observed_solve_block=block)
    raise ValueError('No actual OpenFOAM flow solve evidence at saved pseudo-iteration '+folder.name)

def foam_native(case,folder,log):
    result=foam_mesh(case); count=len(result['cells']); U=foam_field(folder/'U',count,3,[0,1,-1,0,0,0,0])
    require(np.max(np.abs(U[:,2]))<=1e-7,'Unexpected native spanwise velocity')
    result.update(U=U[:,:2],p=foam_field(folder/'p',count,1,[1,-1,-2,0,0,0,0]),T=foam_field(folder/'T',count,1,[0,0,0,1,0,0,0]),layer='one native hexahedron layer projected to source quads',solve_evidence=foam_solve_evidence(log,folder))
    if (folder/'C').is_file(): result['C']=foam_field(folder/'C',count,3,[0,1,0,0,0,0,0])
    if (folder/'rho').is_file(): result['rho']=foam_field(folder/'rho',count,1,[1,-3,0,0,0,0,0])
    return result

def numerical_status(name,path):
    require(path is not None,'Explicit completed numerical analysis required for '+name)
    record=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    checks=record.get('checks',record.get('check_results'))
    require(isinstance(checks,dict) and checks and all(isinstance(v,bool) for v in checks.values()),'No actual numerical check flags in '+str(path))
    failed=[key for key,value in checks.items() if not value]
    require(set(failed)==set(record.get('failed_checks',failed)),'Numerical analysis flags and failure list disagree')
    passed=not failed
    require(record.get('all_declared_checks_passed')==passed,'Numerical analysis aggregate differs from actual flags')
    note=NAMES[name]+(': declared numerical checks passed' if passed else ': '+str(len(failed))+' numerical checks missed')
    manifest_path=Path(path).parent/'source-manifest.json'
    manifest=json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    rows=manifest.get('source_files',manifest.get('sources'))
    require(isinstance(rows,list) and rows,'Actual numerical analysis source manifest required')
    native_sources={str(Path(row['path']).resolve()):row['sha256'] for row in rows}
    require(all(sha(Path(p))==digest for p,digest in native_sources.items()),'Numerical analysis native source hash changed')
    return dict(path=str(Path(path).resolve()),source_manifest_path=str(manifest_path.resolve()),native_source_hashes=native_sources,all_declared_checks_passed=passed,failed_checks=failed,short_note=note)


def common_limits(fields):
    maximum=max(float(f['Mach'].max()) for f in fields.values()); mmax=max(1.5,math.ceil(maximum*5)/5)
    pmin=min(float(f['p'].min()) for f in fields.values())/1000; pmax=max(float(f['p'].max()) for f in fields.values())/1000
    return [0.,mmax],[math.floor(pmin/10)*10,math.ceil(pmax/10)*10]

def decorate(ax,geometry,label):
    ax.plot(geometry.x_edges,geometry.h_edges,color='#223644',lw=1.4,zorder=5)
    ax.plot(geometry.x_edges,-geometry.h_edges,color='#223644',lw=1.4,zorder=5)
    ax.set(xlim=(0,1),ylim=(-.108,.108),xlabel='Axial coordinate, x (m)',ylabel='y (m)')
    ax.set_aspect('equal',adjustable='box'); ax.set_title(label,loc='left',fontsize=11,pad=9,fontweight='semibold')
    ax.tick_params(labelsize=9); ax.spines[['top','right']].set_visible(False); ax.set_facecolor('#eef2f4')

def collection(ax,geometry,values,limits,cmap):
    artist=PolyCollection(geometry.polygons,array=values,cmap=cmap,norm=Normalize(*limits),edgecolors='none',linewidths=0,antialiased=False,rasterized=True)
    ax.add_collection(artist); return artist

def save(fig,out,name):
    fig.canvas.draw(); fig.savefig(out/(name+'.png'),dpi=180,facecolor='white'); fig.savefig(out/(name+'.svg'),dpi=180,facecolor='white'); plt.close(fig)

def render_fields(geometry,fields,out,mach_limits,pressure_limits,status_note):
    fig,axes=plt.subplots(len(fields),2,figsize=(14,7.2),squeeze=False)
    fig.subplots_adjust(left=.065,right=.96,top=.82,bottom=.29,hspace=.65,wspace=.22)
    fig.suptitle('Converging-diverging nozzle | actual saved CFD fields',x=.065,y=.96,ha='left',fontsize=17,fontweight='semibold')
    fig.text(.065,.903,'Inviscid ideal gas | inlet total pressure 200 kPa, total temperature 300 K | back pressure 140 kPa',fontsize=10,color='#405866')
    for row,(name,field) in enumerate(fields.items()):
        for col,(key,limits,cmap,label) in enumerate((('Mach',mach_limits,'viridis','Mach number'),('p',pressure_limits,'magma','Static pressure'))):
            artist=collection(axes[row,col],geometry,field[key]/(1000 if key=='p' else 1),limits,cmap)
            decorate(axes[row,col],geometry,NAMES[name]+' | '+label)
            if row==len(fields)-1:
                cbax=fig.add_axes([.09 if col==0 else .57,.18,.34,.021])
                fig.colorbar(artist,cax=cbax,orientation='horizontal',label='Mach number, M = |U| / sqrt(gamma R T)' if key=='Mach' else 'Absolute static pressure (kPa)')
    fig.text(.065,.07,status_note,fontsize=9,color='#6d3e2c')
    fig.text(.065,.035,'Common geometry and colour scales; native cell values. Appearance does not establish convergence or physical validation.',fontsize=9,color='#405866')
    save(fig,out,'nozzle-Mach-pressure')

def centreline_data(geometry,field):
    values={key:np.array([np.interp(0.,geometry.eta[:,j],field['grid'+key][:,j]) for j in range(geometry.nx)]) for key in ('Mach','p','T')}
    values['x']=geometry.xc; mach,pressure=values['Mach'],values['p']
    candidates=np.flatnonzero((mach[:-1]>1)&(mach[1:]<1)&(np.diff(pressure)>0))
    values['shock_candidates']=[dict(x_left_m=float(geometry.xc[i]),x_right_m=float(geometry.xc[i+1]),upstream_Mach=float(mach[i]),downstream_Mach=float(mach[i+1]),pressure_jump_Pa=float(pressure[i+1]-pressure[i]),description='Adjacent centreline M1 crossing and positive pressure jump; candidate only, not shock validation') for i in candidates]
    return values

def render_centreline(geometry,fields,out,status_note):
    curves={name:centreline_data(geometry,field) for name,field in fields.items()}
    fig,axes=plt.subplots(2,1,figsize=(11.2,6.8),sharex=True); fig.subplots_adjust(left=.095,right=.96,top=.82,bottom=.22,hspace=.22)
    fig.suptitle('Nozzle centreline | native-field comparison',x=.095,y=.96,ha='left',fontsize=17,fontweight='semibold')
    fig.text(.095,.902,'Interpolated between adjacent cell rows at y = 0; no one-dimensional solution imposed',fontsize=10,color='#405866')
    colours={'fluent':'#087c8e','openfoam':'#c65c24'}
    for name,v in curves.items():
        axes[0].plot(v['x'],v['Mach'],label=NAMES[name],color=colours[name],lw=1.8)
        axes[1].plot(v['x'],v['p']/P_TOTAL,label=NAMES[name],color=colours[name],lw=1.8)
    axes[0].axhline(1,color='#7b8990',lw=.85,ls='--'); axes[1].axhline(P_BACK/P_TOTAL,color='#7b8990',lw=.85,ls=':',label='Specified outlet p / p0')
    for ax in axes:
        ax.axvline(.4,color='#a8b2b7',lw=.9,ls='--'); ax.grid(alpha=.16,lw=.6); ax.spines[['top','right']].set_visible(False); ax.set_xlim(0,1)
    axes[0].set_ylabel('Mach number'); axes[1].set(ylabel='Absolute pressure / inlet total pressure',xlabel='Axial coordinate, x (m)')
    axes[0].legend(loc='best',frameon=False,ncol=2); axes[1].legend(loc='best',frameon=False)
    fig.text(.095,.095,status_note,fontsize=9,color='#6d3e2c')
    fig.text(.095,.055,'Dashed vertical line: geometric throat. Computed compression is retained; shock location is not forced or validated.',fontsize=9,color='#405866')
    save(fig,out,'nozzle-centreline'); return curves

def render_animation(geometry,fields,paths,out,mach_limits,frames,status_note):
    fig,axes=plt.subplots(len(fields),1,figsize=(12,8.1),squeeze=False)
    fig.subplots_adjust(left=.075,right=.945,top=.81,bottom=.29,hspace=.45)
    fig.suptitle('Inside a converging-diverging nozzle',x=.075,y=.958,ha='left',fontsize=19,fontweight='semibold')
    fig.text(.075,.905,'ANSYS Fluent + OpenFOAM | inviscid ideal gas | 200 kPa total inlet / 140 kPa back pressure',fontsize=10,color='#405866')
    lines,dots={},{}
    for row,(name,field) in enumerate(fields.items()):
        ax=axes[row,0]; artist=collection(ax,geometry,field['Mach'],mach_limits,'viridis')
        decorate(ax,geometry,NAMES[name]+' | Mach number in saved snapshot')
        for path in paths[name]: ax.plot(path['points'][:,0],path['points'][:,1],color='white',alpha=.10,lw=.55,zorder=3)
        lines[name]=LineCollection([],colors='#fcfffa',linewidths=1.35,alpha=.87,zorder=6); ax.add_collection(lines[name])
        dots[name]=ax.scatter([],[],s=12,c='#fffdf2',edgecolors='#173743',linewidths=.25,zorder=7)
    cbax=fig.add_axes([.225,.17,.60,.024])
    fig.colorbar(artist,cax=cbax,orientation='horizontal',label='Mach number, M = |U| / sqrt(gamma R T) | common scale')
    fig.text(.075,.095,status_note,fontsize=9,color='#6d3e2c')
    fig.text(.075,.055,'Particles move through fixed computed fields; this is not time-resolved CFD evolution.',fontsize=10,color='#304b5a')
    fig.text(.075,.025,'Tracer spacing and playback time are presentation choices. Native mesh, temperature and velocity determine the paths.',fontsize=8.5,color='#506774')
    images=[]
    for frame in range(frames):
        for name in fields:
            segments,heads=[],[]
            for seed_index,path in enumerate(paths[name]):
                total=path['times'][-1]
                for offset in (.0,.33,.66):
                    phase=(frame/frames+seed_index*.071+offset)%1; t=phase*total
                    heads.append([np.interp(t,path['times'],path['points'][:,k]) for k in (0,1)])
                    sample=np.linspace(max(0,t-total*.055),t,12)
                    segments.append(np.column_stack([np.interp(sample,path['times'],path['points'][:,k]) for k in (0,1)]))
            lines[name].set_segments(segments); dots[name].set_offsets(np.array(heads))
        fig.canvas.draw(); image=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy(),'RGB'); images.append(image)
        if frame==frames//2: image.save(out/'nozzle-steady-tracers.png')
    path=out/'nozzle-steady-field-tracers.gif'
    images[0].save(path,save_all=True,append_images=images[1:],duration=50,loop=0,optimize=False,disposal=2); plt.close(fig)
    contact=Image.new('RGB',(images[0].width,(images[0].height+34)*3),'#e8eef1'); draw=ImageDraw.Draw(contact)
    with Image.open(path) as readback:
        for row,index in enumerate((0,frames//2,frames-1)):
            top=row*(images[0].height+34); draw.text((18,top+9),'Decoded GIF frame '+str(index+1)+' / '+str(frames),fill='#193b4b')
            readback.seek(index); contact.paste(readback.convert('RGB'),(0,top+34))
    contact.save(out/'nozzle-steady-tracer-contact-sheet.png')
    with Image.open(path) as readback:
        durations,hashes=[],[]
        for index in range(readback.n_frames):
            readback.seek(index); durations.append(readback.info.get('duration')); hashes.append(hashlib.sha256(readback.convert('RGB').tobytes()).hexdigest())
        require(readback.n_frames==frames and durations==[50]*frames and len(set(hashes))==frames,'GIF frame/timing or distinct decoded frame readback failed')
    return dict(frames=frames,frame_delay_ms=50,playback_fps=20,playback_duration_s=frames/20,actual_saved_field=True,fixed_snapshot=True,solver_transient_time=False,native_velocity_RK4=True,shared_Mach_limits=mach_limits,distinct_decoded_RGB_frames=frames,interpolation='Bilinear in axial and normalized-height coordinates between actual mapped cell-centre velocities; no appended wall or boundary velocity',loops='Fixed phase offsets and recycling along actual integrated trajectories; not flow unsteadiness')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True); parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--geometry',type=Path); parser.add_argument('--fluent-case',type=Path); parser.add_argument('--fluent-data',type=Path)
    parser.add_argument('--fluent-operating-pressure-Pa',type=float,default=0.)
    parser.add_argument('--openfoam-case',type=Path); parser.add_argument('--openfoam-time'); parser.add_argument('--openfoam-log',type=Path)
    parser.add_argument('--fluent-checks',type=Path); parser.add_argument('--openfoam-checks',type=Path)
    parser.add_argument('--frames',type=int,default=48); parser.add_argument('--solver',choices=('both','fluent','openfoam'),default='both'); parser.add_argument('--no-animation',action='store_true')
    args=parser.parse_args(); run,out,script=args.run.resolve(),args.out.resolve(),Path(__file__).resolve()
    require(out.drive.upper()=='E:' and not out.exists(),'Fresh output directory on E: required')
    require(args.frames==48 and args.fluent_operating_pressure_Pa==0.,'This model requires48 frames and operating pressure0Pa')
    geometry_path=(args.geometry or run/'mesh/geometry.npz').resolve(); setup_path=geometry_path.parent/'setup.json'
    metadata=json.loads(setup_path.read_text(encoding='utf-8-sig'))
    for key,value in dict(gamma=GAMMA,gas_constant_J_kg_K=R_GAS,inlet_total_pressure_Pa=P_TOTAL,inlet_total_temperature_K=T_TOTAL,outlet_static_pressure_Pa=P_BACK,operating_pressure_Pa=0.,inviscid=True).items():
        require(metadata.get(key)==value,'Unsupported declared setup value: '+key)
    geometry=QuadGeometry(geometry_path)
    require((geometry.nx,geometry.ny)==(240,40),'Expected240x40 shared quads')
    require(abs(geometry.x_edges[0])<TOL and abs(geometry.x_edges[-1]-1)<TOL and np.allclose([geometry.h_edges[0],geometry.h_edges[96],geometry.h_edges[-1]],[.1,.04,.08],atol=TOL),'Unexpected nozzle walls or throat')
    case=(args.fluent_case or run/'fluent/nozzle-final.cas.h5').resolve(); data=(args.fluent_data or run/'fluent/nozzle-final.dat.h5').resolve()
    foam=(args.openfoam_case or run/'openfoam').resolve(); selected=tuple(NAMES) if args.solver=='both' else (args.solver,)
    for name in selected:
        require((args.fluent_case is not None and args.fluent_data is not None) if name=='fluent' else (args.openfoam_case is not None and args.openfoam_time is not None and args.openfoam_log is not None),'Explicit native case/data/time/log selection required for '+name)
    statuses={name:numerical_status(name,args.fluent_checks if name=='fluent' else args.openfoam_checks) for name in selected}
    status_note=' | '.join(value['short_note'] for value in statuses.values())
    for name,value in statuses.items():
        required_native=[case,data] if name=='fluent' else [foam/args.openfoam_time/p for p in ('U','p','T')]
        require(all(str(p.resolve()) in value['native_source_hashes'] for p in required_native),'Selected fields are not bound to numerical analysis for '+name)
    inputs=[geometry_path,setup_path,script,*(Path(value['path']) for value in statuses.values()),*(Path(value['source_manifest_path']) for value in statuses.values())]; folder,log=None,None
    if 'fluent' in selected: inputs.extend([case,data])
    if 'openfoam' in selected:
        folder=foam_saved_folder(foam,args.openfoam_time)
        log=args.openfoam_log or next((p for p in (foam/'solver.log',foam/'log.solver',foam/'log.foamRun') if p.is_file()),None)
        require(log is not None,'Native OpenFOAM solver log required')
        mesh=next((p for p in (foam/'constant/polyMesh',foam/'0/polyMesh') if (p/'points').is_file()),None)
        require(mesh is not None,'No native OpenFOAM mesh')
        inputs.extend([log,*(mesh/p for p in ('points','faces','owner','neighbour','boundary')),*(folder/p for p in ('U','p','T','C','rho') if (folder/p).is_file())])
    inputs=sorted(set(inputs),key=str); before={str(path):sha(path) for path in inputs}; executed_source=script.read_bytes()
    natives={name:(fluent_native(case,data,args.fluent_operating_pressure_Pa) if name=='fluent' else foam_native(foam,folder,log)) for name in selected}
    fields={name:geometry.match(native) for name,native in natives.items()}; limits=common_limits(fields)
    seeds=[[geometry.xc[0]+1e-9,eta*geometry.height(geometry.xc[0]+1e-9)] for eta in np.linspace(-.84,.84,12)]
    traces={} if args.no_animation else {name:[geometry.trace(seed,field['gridU']) for seed in seeds] for name,field in fields.items()}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.7}); out.mkdir(parents=True)
    render_fields(geometry,fields,out,*limits,status_note); curves=render_centreline(geometry,fields,out,status_note)
    animation=None if args.no_animation else render_animation(geometry,fields,traces,out,limits[0],args.frames,status_note)
    with (out/'actual-native-cell-fields.csv').open('x',encoding='utf-8',newline='') as stream:
        writer=csv.writer(stream); writer.writerow(['solver','common_quad_0based','x_centroid_m','y_centroid_m','u_m_s','v_m_s','absolute_static_pressure_Pa','static_temperature_K','Mach'])
        for name,field in fields.items():
            for index in range(len(geometry.quads)): writer.writerow([name,index,*geometry.centres[index],*field['U'][index],field['p'][index],field['T'][index],field['Mach'][index]])
    with (out/'actual-centreline.csv').open('x',encoding='utf-8',newline='') as stream:
        writer=csv.writer(stream); writer.writerow(['solver','x_m','Mach_interpolated_y0','absolute_pressure_Pa_interpolated_y0','temperature_K_interpolated_y0'])
        for name,values in curves.items(): writer.writerows(zip([name]*geometry.nx,values['x'],values['Mach'],values['p'],values['T']))
    if traces:
        with gzip.open(out/'actual-steady-RK4-traces.csv.gz','wt',encoding='utf-8',newline='') as stream:
            writer=csv.writer(stream); writer.writerow(['solver','seed_0based','integration_travel_time_s','x_m','y_m','speed_m_s'])
            for name,paths in traces.items():
                for index,path in enumerate(paths): writer.writerows([name,index,time,*point,speed] for time,point,speed in zip(path['times'],path['points'],path['speeds']))
    after={str(path):sha(path) for path in inputs}; require(before==after,'Native source or script changed during visualization')
    checks=dict(model='Planar converging-diverging inviscid ideal-gas nozzle',solvers=list(selected),paired=len(selected)==2,actual_native_fields=True,native_cell_to_source_quad_bijective=True,common_cells=len(geometry.quads),common_nodes=len(geometry.nodes),input_hashes_stable_after_loading=True,gamma=GAMMA,R_J_kg_K=R_GAS,inlet_total_pressure_Pa=P_TOTAL,inlet_total_temperature_K=T_TOTAL,outlet_back_pressure_Pa=P_BACK,Mach_shared_limits=limits[0],absolute_static_pressure_shared_limits_kPa=limits[1],animation=animation,numerical_analysis=statuses,native_mapping={name:field['mapping'] for name,field in fields.items()},native_solve_evidence={name:{key:value for key,value in native['solve_evidence'].items() if key!='observed_solve_block'} for name,native in natives.items()},actual_field_ranges={name:{key:[float(field[key].min()),float(field[key].max())] for key in ('Mach','p','T','speed')} for name,field in fields.items()},centreline_shock_candidates={name:values['shock_candidates'] for name,values in curves.items()},tracer_diagnostics={name:[dict(seed=path['seed'],samples=len(path['times']),integration_time_s=float(path['times'][-1]),x_endpoint_m=float(path['points'][-1,0]),stop=path['stop']) for path in paths] for name,paths in traces.items()},imposed_shock=False,synthetic_flow=False,physical_validation=False,mesh_independence=False,quantitative_acceptance='Separate native numerical checks; visual appearance is not convergence evidence')
    (out/'checks.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if 'openfoam' in natives: (out/'observed-openfoam-solve-block.txt').write_text(natives['openfoam']['solve_evidence']['observed_solve_block'],encoding='utf-8')
    (out/'executed-plotter.py').write_bytes(executed_source)
    (out/'input-hash-manifest.json').write_text(json.dumps(dict(command_arguments=sys.argv[1:],plotter_sha256=before[str(script)],inputs=[dict(path=str(path),bytes=path.stat().st_size,sha256=before[str(path)]) for path in inputs]),indent=2)+'\n',encoding='utf-8')
    caption=('''# Actual native nozzle visualization

Planar converging-diverging nozzle, length 1 m, full inlet/throat/exit heights 0.20/0.08/0.16 m, throat x = 0.4 m. Both codes solve inviscid ideal-gas Euler flow with gamma = 1.4, R = 287.05 J/(kg K), total inlet pressure 200 kPa, total inlet temperature 300 K and absolute outlet back pressure 140 kPa. Shared 240 x 40 wall-fitted quads. Numerical acceptance is recorded separately.

Native Fluent cells are reconstructed from saved face connectivity. Native OpenFOAM hexahedra are reconstructed from saved faces/owner/neighbour and projected to source quads. Sorted vertex sets give a bijective common-cell map; imported ordering is not assumed. Saved OpenFOAM C, when available, independently checks true polygon centroids. Native pressure is absolute (Fluent operating pressure 0 Pa). Native temperature and velocity give Mach = |U| / sqrt(gamma R T). No one-dimensional solution, shock position, wall velocity or synthetic flow is inserted. Both codes use common unclipped colour scales.

Centreline curves interpolate between the two adjacent actual cell rows at y = 0. An adjacent M > 1 to M < 1 transition and positive pressure jump is only a computed compression/shock candidate. It does not validate shock location or strength. The CSV preserves all samples.

Animation particles follow RK4 integration through velocity interpolated bilinearly in axial and normalized-height coordinates using mapped cell centres. No boundary velocity values are appended. Paths stop at the native cell-centre interpolation envelope or actual wall. Recycled phases and playback timing move particles through each fixed saved snapshot. This is not physical solver-time evolution; OpenFOAM localEuler is a pseudo-time steady-convergence procedure. Exported integration travel times differ from playback and solver pseudo-iterations.

ANSYS Fluent and OpenFOAM computed the fields; Codex prepared mesh, workflow and visualization. No user reflection, independent user execution, mesh/domain independence, shock accuracy, benchmark agreement or physical validation is claimed. Source hashes and executed source accompany the figures; all inputs are rehashed after output.
''')
    caption+='\nActual numerical status: '+status_note+'.\n'
    (out/'CAPTIONS.md').write_text(caption,encoding='utf-8')
    outputs=sorted((path for path in out.rglob('*') if path.is_file()),key=str)
    (out/'output-hash-manifest.json').write_text(json.dumps([dict(path=str(path.relative_to(out)),bytes=path.stat().st_size,sha256=sha(path)) for path in outputs],indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(out),cells=len(geometry.quads),solvers=list(fields),animation=animation,actual_field_ranges=checks['actual_field_ranges'],shock_candidates=checks['centreline_shock_candidates']),indent=2))
    return 0

if __name__=='__main__': sys.exit(main())
