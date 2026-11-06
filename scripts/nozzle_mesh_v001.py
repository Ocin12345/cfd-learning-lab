"""Generate one small full-planar CD nozzle; no solver launch."""
from pathlib import Path
import argparse, hashlib, json, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,required=True); args=ap.parse_args()
    out=args.out.resolve(); repo=Path(__file__).resolve().parents[1]
    if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo): raise ValueError('Fresh external E: output required')
    if shutil.disk_usage(out.anchor).free<10*1024**3: raise RuntimeError('Require10GiB free E:')
    nx,ny=240,40
    x=np.linspace(0.,1.,nx+1)
    h=np.where(x<=.4,.04+(.1-.04)*.5*(1+np.cos(np.pi*x/.4)),.04+(.08-.04)*.5*(1-np.cos(np.pi*(x-.4)/.6)))
    eta=np.linspace(-1.,1.,ny+1)
    nodes=np.array([(xx,ee*hh) for ee in eta for xx,hh in zip(x,h)],dtype=float)
    def nid(i,j): return j*(nx+1)+i
    quads=np.array([[nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)] for j in range(ny) for i in range(nx)],dtype=np.int64)
    polys=nodes[quads]; nxt=np.roll(polys,-1,axis=1)
    cross=polys[:,:,0]*nxt[:,:,1]-nxt[:,:,0]*polys[:,:,1]
    areas=.5*cross.sum(axis=1)
    centres=((polys+nxt)*cross[:,:,None]).sum(axis=1)/(6*areas[:,None])
    assert np.all(areas>0) and abs(areas.sum()-.128)<1e-13
    assert np.allclose(nodes[:nx+1,1],-h,rtol=0,atol=1e-15)
    faces={}
    for cid,quad in enumerate(quads,1):
        for a,b in zip(quad,np.roll(quad,-1)):
            key=tuple(sorted((int(a),int(b))))
            if key in faces:
                assert faces[key][3]==0; faces[key][3]=cid
            else: faces[key]=[int(a)+1,int(b)+1,cid,0]
    groups={z:[] for z in (3,4,5,6,7)}; orth=[]
    for a,b,c0,c1 in faces.values():
        p0,p1=nodes[a-1],nodes[b-1]; midpoint=.5*(p0+p1)
        norm=np.array([p1[1]-p0[1],-(p1[0]-p0[0])])
        if c1:
            zone=3; delta=centres[c1-1]-centres[c0-1]
        else:
            delta=midpoint-centres[c0-1]
            if p0[0]==p1[0]==0: zone=6
            elif p0[0]==p1[0]==1: zone=7
            elif min(p0[1],p1[1])>0: zone=4
            else: assert max(p0[1],p1[1])<0; zone=5
        dot=float(np.dot(norm,delta)/(np.linalg.norm(norm)*np.linalg.norm(delta)))
        assert dot>0; orth.append(dot); groups[zone].append((a,b,c0,c1))
    assert len(nodes)-len(faces)+len(quads)==1
    assert [len(groups[z]) for z in (4,5,6,7)]==[nx,nx,ny,ny]
    out.mkdir(parents=True)
    mesh=out/'nozzle.msh'
    with mesh.open('w',encoding='ascii',newline='\n') as f:
        f.write('(0 "Full planar ideal-gas CD nozzle; 9600 shared quadrilaterals")\n(2 2)\n')
        f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
        for xx,yy in nodes: f.write(f'{xx:.17g} {yy:.17g}\n')
        f.write('))\n'); f.write(f'(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
        start=1
        for zone,bc in [(3,2),(4,3),(5,3),(6,4),(7,5)]:
            rows=groups[zone]; f.write(f'(13 ({zone:x} {start:x} {start+len(rows)-1:x} {bc:x} 2)(\n')
            for row in rows: f.write(' '.join(f'{v:x}' for v in row)+'\n')
            f.write('))\n'); start+=len(rows)
        for z,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','upper_wall'),(5,'wall','lower_wall'),(6,'pressure-inlet','inlet'),(7,'pressure-outlet','outlet')]: f.write(f'(45 ({z} {kind} {name})())\n')
    np.savez_compressed(out/'geometry.npz',nodes=nodes,quads=quads,centres=centres,areas=areas,grid_shape=np.array([ny,nx]),ix=np.tile(np.arange(nx),ny),iy=np.repeat(np.arange(ny),nx),x_edges=x,half_height_edges=h,eta_edges=eta)
    setup=dict(model='Full planar converging-diverging nozzle; illustrative inviscid ideal-gas Euler flow',nx=nx,ny=ny,cells=nx*ny,length_m=1.,throat_x_m=.4,inlet_half_height_m=.1,throat_half_height_m=.04,exit_half_height_m=.08,reference_depth_m=1.,openfoam_extrusion_depth_m=.01,gamma=1.4,gas_constant_J_kg_K=287.05,cp_J_kg_K=1004.675,molecular_weight_kg_kmol=8314.46261815324/287.05,inlet_total_pressure_Pa=200000.,inlet_total_temperature_K=300.,outlet_static_pressure_Pa=140000.,operating_pressure_Pa=0.,initial_pressure_Pa=140000.,initial_temperature_K=300.,initial_velocity_m_s=[50.,0.],inviscid=True,wall_condition='impermeable slip, adiabatic',boundary_names=['inlet','outlet','upper_wall','lower_wall'],mesh_sha256=sha(mesh),geometry_npz_sha256=sha(out/'geometry.npz'),geometry_checks=dict(positive_areas=True,total_area_m2=float(areas.sum()),minimum_area_m2=float(areas.min()),euler_characteristic=1,face_orthogonality_proxy_min=float(min(orth)),boundary_face_counts={name:len(groups[z]) for z,name in [(4,'upper_wall'),(5,'lower_wall'),(6,'inlet'),(7,'outlet')]}),generated_by='Codex',physical_validation=False)
    (out/'setup.json').write_text(json.dumps(setup,indent=2)+'\n',encoding='utf-8')
    (out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
    fig,ax=plt.subplots(figsize=(12,3.6),layout='constrained')
    for j in range(0,ny+1,4): ax.plot(x,eta[j]*h,color='#8fb0c2',lw=.45)
    for i in range(0,nx+1,4): ax.plot([x[i],x[i]],[-h[i],h[i]],color='#8fb0c2',lw=.45)
    ax.plot(x,h,color='#172a3a',lw=1.8);ax.plot(x,-h,color='#172a3a',lw=1.8)
    ax.axvline(.4,color='#b44428',ls='--',lw=1);ax.annotate('Throat',(0.4,.045),xytext=(.43,.13),arrowprops=dict(arrowstyle='->',lw=.8))
    ax.set(xlabel='x (m)',ylabel='y (m)',title='Converging-diverging nozzle · actual shared mesh',xlim=(0,1),ylim=(-.14,.16));ax.set_aspect('equal');ax.spines[['top','right']].set_visible(False)
    fig.text(.5,.01,'Geometry only · 9,600 planar cells · displayed grid every fourth line · no computed flow',ha='center',fontsize=9,color='#4d5963')
    fig.savefig(out/'geometry-preview.png',dpi=160);plt.close(fig)
    manifest={p.name:sha(p) for p in out.iterdir() if p.is_file()};(out/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(setup,indent=2))

if __name__=='__main__': main()