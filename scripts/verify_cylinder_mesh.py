"""Check an extended cylinder mesh against its preserved radius-20 source."""
from pathlib import Path
import argparse,hashlib,json,math

def node_lines(path):
    lines=path.read_text().splitlines()
    start=next(i for i,line in enumerate(lines) if line.startswith('(10 (1 '))+1
    end=next(i for i in range(start,len(lines)) if lines[i]=='))')
    return lines[start:end]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--mesh',type=Path,required=True)
    args=parser.parse_args()
    before=node_lines(args.source)
    after=node_lines(args.mesh)
    setup=json.loads((args.mesh.parent/'setup.json').read_text())
    nt=setup['circumferential_cells']
    assert len(before)==nt*(setup['preserved_radius20_radial_layers']+1)
    assert after[:len(before)]==before,'Original core coordinates changed'
    assert len(after)==nt*(setup['radial_layers']+1)
    xy=[tuple(map(float,line.split())) for line in after]
    assert max(abs(math.hypot(x,y)-setup['outer_radius_m']) for x,y in xy[-nt:])<1e-10
    min_area=math.inf
    for j in range(setup['radial_layers']):
        for i in range(nt):
            polygon=[xy[j*nt+i],xy[(j+1)*nt+i],xy[(j+1)*nt+(i+1)%nt],xy[j*nt+(i+1)%nt]]
            area=sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(polygon,polygon[1:]+polygon[:1]))/2
            min_area=min(min_area,area)
    assert min_area>0,'Non-positive oriented cell area'
    record=dict(original_nodes=len(before),extended_nodes=len(after),
                preserved_nodes_identical=True,
                preserved_coordinates_sha256=hashlib.sha256('\n'.join(before).encode()).hexdigest(),
                outer_radius_m=setup['outer_radius_m'],cells=setup['cells'],minimum_oriented_area_m2=min_area)
    (args.mesh.parent/'mesh-preservation-check.json').write_text(json.dumps(record,indent=2))
    print(json.dumps(record,indent=2))

if __name__=='__main__':main()
