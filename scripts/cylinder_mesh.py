from pathlib import Path

import math,json,argparse,os

parser=argparse.ArgumentParser()

parser.add_argument('--name',default='ansys-cylinder')

parser.add_argument('--nt',type=int,default=256)

parser.add_argument('--nr',type=int,default=96)
parser.add_argument('--outer-radius',type=float,default=20.0)

args=parser.parse_args()

output_root=Path(os.environ['CFD_OUTPUT_ROOT']).resolve()
out=(output_root/args.name).resolve()
if not out.is_relative_to(output_root) or output_root.is_relative_to(Path(__file__).resolve().parents[1]):
    raise ValueError('Outputs must stay in an external run root')

out.mkdir(parents=True,exist_ok=False)

nt,nr=args.nt,args.nr
if nt < 8 or nt%4 or nr < 2 or not math.isfinite(args.outer_radius) or args.outer_radius < 20:
    raise ValueError('Use nt divisible by four, nr >= 2 and outer radius >= 20 m')
radii=[.5+19.5*math.expm1(5*j/nr)/math.expm1(5) for j in range(nr+1)]
core_nr=nr
if args.outer_radius > 20:
    growth=math.exp(5/nr)
    last_spacing=radii[-1]-radii[-2]
    distance=args.outer_radius-20
    count=math.ceil(math.log1p(distance*(growth-1)/(last_spacing*growth))/math.log(growth))
    steps=[last_spacing*growth**(i+1) for i in range(count)]
    scale=distance/sum(steps)
    for step in steps:radii.append(radii[-1]+scale*step)
    radii[-1]=args.outer_radius
    nr=len(radii)-1

nodes=[]

for j in range(nr+1):

    r=radii[j]

    for i in range(nt):

        a=2*math.pi*i/nt

        nodes.append((r*math.cos(a),r*math.sin(a)))

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

    else: zone=5 if nodes[a-1][0]+nodes[b-1][0]<0 else 6

    groups[zone].append((a,b,c,d))

with (out/'cylinder.msh').open('w') as f:

    f.write(f'(0 "Cylinder O-grid: D=1m, outer radius={args.outer_radius:g}m")\n(2 2)\n')

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

    for zone,typ,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','cylinder'),(5,'velocity-inlet','inlet'),(6,'pressure-outlet','outlet')]:

        f.write(f'(45 ({zone} {typ} {name})())\n')

(out/'setup.json').write_text(json.dumps(dict(diameter_m=1,outer_radius_m=args.outer_radius,density_kg_m3=1,viscosity_Pa_s=.05,inlet_velocity_m_s=1,reynolds_number=20,cells=nt*nr,radial_layers=nr,circumferential_cells=nt,preserved_radius20_radial_layers=core_nr,inner_nodes_preserved=True),indent=2))

print(out/'cylinder.msh')
