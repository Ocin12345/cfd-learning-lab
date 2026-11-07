"""Read native wall fluxes, integrate geometry and cross-check solver reports."""
from pathlib import Path
import argparse,json,re,hashlib
import h5py,numpy as np
from plot_nozzle_native_v004 import indexed,foam_text,foam_list
from plume_production_common_v001 import ROOT,K,DEPTH
from render_native_plume_v002 import native_fluent_time,native_foam_time

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def scalar_patch(path,name,n,dimensions):
    raw=foam_text(path);d=re.search(r'dimensions\s*\[([^]]+)\]',raw)
    assert d and np.array_equal(np.fromstring(d[1],sep=' '),dimensions),(path,'dimensions')
    m=re.search(r'\b'+re.escape(name)+r'\s*\{([^{}]*)\}',raw,re.S);assert m,(path,name)
    body=m[1];v=re.search(r'\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;',body,re.S)
    if v:
        assert int(v[1])==n;values=np.fromstring(v[2],sep=' ')
    else:
        v=re.search(r'\bvalue\s+uniform\s+([-+0-9.eE]+)\s*;',body);assert v,(path,name)
        values=np.full(n,float(v[1]))
    assert values.shape==(n,) and np.isfinite(values).all()
    return values

def fluent_walls(F,names,t):
    case=F/'qualification-final.cas.h5';data=F/'qualification-final.dat.h5';report=F/'wall-heat-final.txt'
    assert abs(native_fluent_time(data)['time_s']-t)<1e-8
    raw=report.read_text();result={}
    with h5py.File(case) as c,h5py.File(data) as d:
        mesh=c['meshes/1'];z=mesh['faces/zoneTopology'];nf=int(max(z['maxId']))
        xyz=mesh['nodes/coords'];nn=max(int(np.asarray(a.attrs['maxId']).ravel()[0]) for a in xyz.values());nodes=indexed(xyz,nn,2)
        faces=np.full((nf,2),-1,dtype=int)
        for group in mesh['faces/nodes'].values():
            a=group['nnodes'];assert np.all(a[:]==2)
            start=int(np.asarray(a.attrs.get('minId',group.attrs.get('minId',[1]))).ravel()[0])-1
            v=group['nodes'][:].reshape(-1,2).astype(int)-1;assert np.all(faces[start:start+len(v)]==-1);faces[start:start+len(v)]=v
        assert np.all(faces>=0) and np.all(faces<nn)
        lengths=np.linalg.norm(nodes[faces[:,1]]-nodes[faces[:,0]],axis=1);assert np.all(lengths>0)
        q=indexed(d['results/1/phase-1/faces/SV_HEAT_FLUX'],nf,partial=True);mass=indexed(d['results/1/phase-1/faces/SV_FLUX'],nf)
        zn=z['name'][0].decode().split(';');assert set(names).issubset(zn)
        for name,start,end in zip(zn,z['minId'],z['maxId']):
            if name not in names:continue
            sl=slice(int(start)-1,int(end));v=q[sl];m=mass[sl];assert np.isfinite(v).all()
            reported=re.findall(r'^\s*'+re.escape(name)+r'\s+([-+0-9.eE]+)\s*$',raw,re.M);assert len(reported)==1
            Q=float(v.sum());assert np.isclose(Q,float(reported[0]),rtol=1e-7,atol=1e-10),(name,Q,reported)
            result[name]=dict(face_count=len(v),wall_length_m=float(lengths[sl].sum()),heat_per_depth_W_m=Q,report_heat_per_depth_W_m=float(reported[0]),absolute_mass_flux_per_depth_kg_ms=float(abs(m).sum()),net_mass_flux_per_depth_kg_ms=float(m.sum()),native_heat_representation='per-face heat rate, independently matched to Fluent Flux Report')
    return result,[case,data,report]

def foam_walls(O,names,t):
    mesh=O/'constant/polyMesh'
    if not mesh.exists():mesh=O/'0/polyMesh'
    n,body=foam_list(mesh/'points');points=np.fromstring(body.replace('(',' ').replace(')',' '),sep=' ').reshape(n,3)
    nf,body=foam_list(mesh/'faces');matches=re.findall(r'(\d+)\s*\(([^)]*)\)',body);assert len(matches)==nf
    area=[]
    for count,ids in matches:
        ids=np.fromstring(ids,sep=' ',dtype=int);assert len(ids)==int(count)
        assert np.all(ids>=0) and np.all(ids<n)
        v=points[ids];area.append(float(np.linalg.norm(np.cross(v,np.roll(v,-1,axis=0)).sum(axis=0))*.5))
    area=np.array(area);assert np.all(area>0)
    boundary=foam_text(mesh/'boundary');p=O/f'{t:g}';native_foam_time(p/'uniform/time',t)
    report=O/'wall-heat-postprocess.log';raw=report.read_text();assert re.search(r'^End\s*$',raw,re.M)
    result={}
    for name in names:
        block=re.search(r'\b'+re.escape(name)+r'\s*\{([^{}]*)\}',boundary,re.S);assert block
        count=int(re.search(r'nFaces\s+(\d+)',block[1])[1]);start=int(re.search(r'startFace\s+(\d+)',block[1])[1]);a=area[start:start+count];assert len(a)==count
        q=scalar_patch(p/'wallHeatFlux',name,count,[1,0,-3,0,0,0,0]);flux=scalar_patch(p/'phi',name,count,[1,0,-1,0,0,0,0])
        log=re.findall(r'for patch '+re.escape(name)+r' = ([-+0-9.eE]+), ([-+0-9.eE]+), ([-+0-9.eE]+), ([-+0-9.eE]+)',raw);assert len(log)==1
        Q=float(q@a);reported=float(log[0][2]);assert np.isclose(Q,reported,rtol=1e-9,atol=1e-13),(name,Q,reported)
        assert np.isclose(Q/a.sum(),float(log[0][3]),rtol=1e-9,atol=1e-10)
        result[name]=dict(face_count=count,wall_length_m=float(a.sum()/DEPTH),heat_per_depth_W_m=Q/DEPTH,report_heat_per_depth_W_m=reported/DEPTH,absolute_mass_flux_per_depth_kg_ms=float(abs(flux).sum()/DEPTH),net_mass_flux_per_depth_kg_ms=float(flux.sum()/DEPTH),native_heat_representation='W/m2 integrated over native polygon areas and divided by extrusion depth')
    return result,[mesh/'points',mesh/'faces',mesh/'boundary',p/'wallHeatFlux',p/'phi',p/'uniform/time',report]

def main():
    p=argparse.ArgumentParser();p.add_argument('--reference',action='store_true');p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    assert args.out.is_relative_to(ROOT) and not args.out.exists()
    base=ROOT/'native/runs'/('natural-convection-dt001-v001' if args.reference else 'hot-patch-plume-20s-v001')
    F=base/('fluent-v002' if args.reference else 'fluent');O=base/'openfoam';names=['top','bottom','left','right']+([] if args.reference else ['heater'])
    hot='left' if args.reference else 'heater';deltaT=.8 if args.reference else 1.
    assert json.loads((F/'status.json').read_text())['stage']=='closed'
    assert (O/'run-status.txt').read_text().strip()=='complete'
    f,fs=fluent_walls(F,names,20.);o,os=foam_walls(O,names,20.)
    for name in names:
        assert f[name]['face_count']==o[name]['face_count']
        assert abs(f[name]['wall_length_m']-o[name]['wall_length_m'])<1e-10
    result={}
    for solver,walls in [('ANSYS Fluent',f),('OpenFOAM',o)]:
        leak=sum(w['absolute_mass_flux_per_depth_kg_ms'] for w in walls.values());assert leak<1e-12,leak
        L=walls[hot]['wall_length_m'];Q=walls[hot]['heat_per_depth_W_m'];assert Q>0
        result[solver]=dict(walls=walls,total_absolute_wall_mass_flux_per_depth_kg_ms=leak,heated_wall_length_m=L,heated_wall_heat_per_depth_W_m=Q,heated_wall_mean_flux_W_m2=Q/L,Nu_based_on_heated_wall_length=Q/(K*deltaT),net_wall_heat_per_depth_W_m=sum(w['heat_per_depth_W_m'] for w in walls.values()))
    record=dict(study='completed square-cavity reference' if args.reference else 'localized floor-heater plume',time_s=20.,checks_pass=True,results=result,sign_convention='positive heat into fluid, verified against native solver reports',Nu_definition='q_mean * heated_wall_length / (conductivity * imposed_hot_cold_temperature_difference)',conductivity_W_mK=K,deltaT_K=deltaT,openfoam_depth_m=DEPTH,fluent_planar_depth_m=1.,energy_balance_claimed=False,note='Instantaneous wall heat rates and no-penetration checks; not a transient global energy balance or physical validation.',native_source_sha256={str(v):digest(v) for v in fs+os+[Path(__file__)]})
    args.out.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='native_source_sha256'},indent=2))
if __name__=='__main__':main()
