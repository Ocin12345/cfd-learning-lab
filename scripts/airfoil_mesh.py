from pathlib import Path

import math,json,argparse,os

parser=argparse.ArgumentParser()

parser.add_argument('--name',default='ansys-airfoil/alpha0')

parser.add_argument('--alpha',type=float,default=0)

parser.add_argument('--nt',type=int,default=384)

parser.add_argument('--nr',type=int,default=160)

args=parser.parse_args()

out=Path(os.environ['CFD_OUTPUT_ROOT']).resolve()/args.name

out.mkdir(parents=True,exist_ok=True)

nt,nr=args.nt,args.nr

nodes=[]

alpha=math.radians(args.alpha)

extent=1.008930411365

for j in range(nr+1):

    t=math.expm1(6*j/nr)/math.expm1(6)

    for i in range(nt):

        theta=2*math.pi*i/nt

        x=.5*(1+math.cos(theta))

        xx=x*extent

        y=.6*(.2969*math.sqrt(xx)-.1260*xx-.3516*xx**2+.2843*xx**3-.1015*xx**4)/extent

        y=math.copysign(max(0,y),math.sin(theta))

        # Clockwise rotation gives positive incidence relative to horizontal inflow.

        sx=.5+(x-.5)*math.cos(alpha)+y*math.sin(alpha)

        sy=-(x-.5)*math.sin(alpha)+y*math.cos(alpha)

        ox=.5+20*math.cos(theta-alpha)

        oy=20*math.sin(theta-alpha)

        nodes.append(((1-t)*sx+t*ox,(1-t)*sy+t*oy))

def nid(j,i): return j*nt+i%nt+1

faces={}

for j in range(nr):

    for i in range(nt):

        cid=j*nt+i+1

        # Counterclockwise cells, c0 on the left of each directed edge.

        v=[nid(j,i),nid(j+1,i),nid(j+1,i+1),nid(j,i+1)]

        for a,b in zip(v,v[1:]+v[:1]):

            key=tuple(sorted((a,b)))

            if key in faces: faces[key][3]=cid

            else: faces[key]=[a,b,cid,0]

groups={3:[],4:[],5:[],6:[]}

for a,b,c,d in faces.values():

    if d: zone=3

    elif a<=nt and b<=nt: zone=4

    else: zone=5 if nodes[a-1][0]+nodes[b-1][0]<1 else 6

    groups[zone].append((a,b,c,d))

with (out/'airfoil.msh').open('w') as f:

    f.write('(0 "Scaled closed NACA0012 O-grid: chord=1m, outer radius=20m")\n(2 2)\n')

    f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')

    for x,y in nodes:f.write(f'{x:.14g} {y:.14g}\n')

    f.write('))\n')

    f.write(f'(12 (0 1 {nt*nr:x} 0))\n(12 (2 1 {nt*nr:x} 1 3))\n')

    f.write(f'(13 (0 1 {len(faces):x} 0))\n')

    start=1

    for zone,bc in [(3,2),(4,3),(5,10),(6,5)]:

        ff=groups[zone]

        f.write(f'(13 ({zone:x} {start:x} {start+len(ff)-1:x} {bc:x} 2)(\n')

        for face in ff:f.write(' '.join(f'{v:x}' for v in face)+'\n')

        f.write('))\n');start+=len(ff)

    for zone,typ,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','airfoil'),(5,'velocity-inlet','inlet'),(6,'pressure-outlet','outlet')]:

        f.write(f'(45 ({zone} {typ} {name})())\n')

(out/'setup.json').write_text(json.dumps(dict(chord_m=1,angle_degrees=args.alpha,geometry="scaled closed NACA0012",outer_radius_m=20,density_kg_m3=1,viscosity_Pa_s=.01,inlet_velocity_m_s=1,reynolds_number=100,cells=nt*nr,radial_layers=nr,circumferential_cells=nt),indent=2))

print(out/'airfoil.msh')
