"""Unchanged dam-break geometry/physics with a smaller 0.05ms fixed time step."""
from pathlib import Path
import json,re
import numpy as np
L=H=.584
OB_LEFT=.292;OB_RIGHT=.31599948;OB_HEIGHT=.04799896
WATER_WIDTH=.1461;WATER_HEIGHT=.292
DEPTH=.0146;DT=.00005;G=9.81;SIGMA=.07
MATERIALS={'air':(1.,1.48e-5),'water-liquid':(1000.,.001)}
NX=92;NY=100;NCELLS=9072
WATER_AREA=WATER_WIDTH*WATER_HEIGHT
FLUID_AREA=L*H-(OB_RIGHT-OB_LEFT)*OB_HEIGHT
FLAT_HEIGHT=(WATER_AREA+(OB_RIGHT-OB_LEFT)*OB_HEIGHT)/L
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def geometry():
    x=np.r_[np.linspace(0,OB_LEFT,47),np.linspace(OB_LEFT,OB_RIGHT,9)[1:],np.linspace(OB_RIGHT,L,39)[1:]]
    y=np.r_[np.linspace(0,OB_HEIGHT,17),np.linspace(OB_HEIGHT,H,85)[1:]]
    original=np.array([(xx,yy) for yy in y for xx in x]);quads=[];ij=[]
    for j in range(NY):
        for i in range(NX):
            if 46<=i<54 and j<16:continue
            k=j*(NX+1)+i;quads.append((k,k+1,k+NX+2,k+NX+1));ij.append((i,j))
    quads=np.asarray(quads);used=np.unique(quads);remap=np.full(len(original),-1,int);remap[used]=np.arange(len(used));nodes=original[used];quads=remap[quads];ij=np.asarray(ij)
    centres=nodes[quads].mean(axis=1);bounds=np.column_stack([x[ij[:,0]],x[ij[:,0]+1],y[ij[:,1]],y[ij[:,1]+1]])
    areas=(bounds[:,1]-bounds[:,0])*(bounds[:,3]-bounds[:,2]);assert len(areas)==NCELLS and len(nodes)==9281 and np.all(areas>0)
    assert abs(areas.sum()-FLUID_AREA)<1e-13
    ids=np.full((NY,NX),-1,int);ids[ij[:,1],ij[:,0]]=np.arange(NCELLS)
    blocks=[(0,46,0,16),(54,92,0,16),(0,46,16,100),(46,54,16,100),(54,92,16,100)]
    foam_order=np.concatenate([ids[j0:j1,i0:i1].ravel() for i0,i1,j0,j1 in blocks]);assert np.array_equal(np.sort(foam_order),np.arange(NCELLS))
    return dict(nodes=nodes,quads=quads,cell_ij=ij,centres=centres,areas=areas,bounds=bounds,x_edges=x,y_edges=y,openfoam_to_common=foam_order)
def initial_fields(g,mode):
    assert mode in ['flat','collapse','production'];b=g['bounds'];dx=b[:,1]-b[:,0];dy=b[:,3]-b[:,2]
    height=FLAT_HEIGHT if mode=='flat' else WATER_HEIGHT
    xf=np.ones(NCELLS) if mode=='flat' else np.clip((WATER_WIDTH-b[:,0])/dx,0,1)
    yf=np.clip((height-b[:,2])/dy,0,1);alpha=xf*yf
    rho=MATERIALS['air'][0]+alpha*(MATERIALS['water-liquid'][0]-MATERIALS['air'][0])
    # Identical gauge-pressure input in both solvers; rho_operating=0 in Fluent.
    pressure=MATERIALS['air'][0]*G*(H-g['centres'][:,1])+(MATERIALS['water-liquid'][0]-MATERIALS['air'][0])*G*xf*np.maximum(height-g['centres'][:,1],0)
    reduced=pressure+rho*G*g['centres'][:,1]
    area=float(alpha@g['areas']);assert abs(area-WATER_AREA)<1e-13
    return dict(alpha=alpha,velocity=np.zeros((NCELLS,3)),pressure=pressure,p_rgh=reduced,rho=rho,water_area_m2=area,mixture_mass_per_m=float(rho@g['areas']))
def make_mesh(out):
    out=Path(out);assert out.is_dir() and not(out/'tank.msh').exists() and not(out/'geometry.npz').exists();g=geometry();faces={}
    for cid,q in enumerate(g['quads'],1):
        for a,b in zip(q,np.roll(q,-1)):
            key=tuple(sorted((int(a),int(b))))
            if key in faces:assert faces[key][3]==0;faces[key][3]=cid
            else:faces[key]=[int(a)+1,int(b)+1,cid,0]
    interior=[];walls=[];atmosphere=[]
    for face in faces.values():
        if face[3]:interior.append(face)
        elif np.all(abs(g['nodes'][np.array(face[:2])-1,1]-H)<1e-12):atmosphere.append(face)
        else:walls.append(face)
    assert len(atmosphere)==92 and len(walls)==324 and len(faces)==18352
    path=out/'tank.msh'
    with path.open('x',encoding='ascii',newline='\n') as f:
        f.write('(0 "Model35 nonuniform 9072-cell dam-break mesh")\n(2 2)\n')
        f.write(f'(10 (0 1 {len(g["nodes"]):x} 1 2))\n(10 (1 1 {len(g["nodes"]):x} 1 2)(\n')
        for xx,yy in g['nodes']:f.write(f'{xx:.17g} {yy:.17g}\n')
        f.write(f'))\n(12 (0 1 {NCELLS:x} 0))\n(12 (2 1 {NCELLS:x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
        first=1
        for zone,bc,rows in [(3,2,interior),(4,3,walls),(5,5,atmosphere)]:
            f.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
            for row in rows:f.write(' '.join(f'{v:x}' for v in row)+'\n')
            f.write('))\n');first+=len(rows)
        f.write('(45 (2 fluid fluid)())\n(45 (3 interior interior)())\n(45 (4 wall walls)())\n(45 (5 pressure-outlet atmosphere)())\n')
    np.savez_compressed(out/'geometry.npz',**g)
    return path
def native_field(path,components,n=NCELLS):
    text=Path(path).read_text(encoding='ascii');m=re.search(r'internalField\s+uniform\s+([^;]+);',text)
    if m:
        v=np.fromstring(m[1].replace('(',' ').replace(')',' '),sep=' ');assert len(v)==components;return np.tile(v,(n,1))
    m=re.search(r'internalField\s+nonuniform\s+List<\w+>\s+(\d+)\s*\(',text);assert m and int(m[1])==n,path
    end=re.search(r'\n\s*\)\s*;',text[m.end():]);assert end,path
    v=np.fromstring(text[m.end():m.end()+end.start()].replace('(',' ').replace(')',' '),sep=' ');assert v.size==n*components and np.isfinite(v).all(),path;return v.reshape(n,components)
def cell_values(group):
    result=np.full(NCELLS,np.nan)
    for ds in group.values():
        lo=int(np.asarray(ds.attrs['minId']).flat[0])-1;hi=int(np.asarray(ds.attrs['maxId']).flat[0]);assert 0<=lo<hi<=NCELLS and np.isnan(result[lo:hi]).all();result[lo:hi]=np.asarray(ds).ravel()
    assert np.isfinite(result).all();return result
