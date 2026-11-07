"""Check matched physical strip mapping against actual accepted initial fields."""
from pathlib import Path
import json,hashlib,re
import h5py,numpy as np
from analyse_nozzle_fluent_v005 import common_mapping
from compare_coarse_fine_v001 import chunk,field
R=Path('E:/CFD-Learning-Lab');P=R/'local/sloshing-refinement-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    results={};sources=[]
    for solver in ['fluent','openfoam']:
        rows=[]
        for fine in [False,True]:
            nx,ny,dx=(160,80,.00125) if fine else (80,40,.0025);n=nx*ny
            if fine:run=R/'native/runs/sloshing-refinement-v001/qualification'/(solver+'-tilted');audit=P/(solver+'-tilted-audit.json')
            else:
                run=R/'native/runs/sloshing-production-v001'/('fluent' if solver=='fluent' else 'openfoam-v003');audit=R/'local/sloshing-production-v001'/('fluent-production-audit-v001.json' if solver=='fluent' else 'openfoam-production-audit-v003.json')
            a=json.loads(audit.read_text());assert a['numerical_checks_pass'] if fine else a['production_checks_pass']
            paths=[]
            if solver=='fluent':
                paths=[run/'tank-0000.dat.h5',run/'tank.cas.h5',run/'geometry.npz']
                with np.load(paths[2]) as g:
                    # Coarse audit did not hash this helper file: prove its exact
                    # row-major geometry against the declared grid and bound case.
                    expected_nodes=np.array([(i*dx,j*dx) for j in range(ny+1) for i in range(nx+1)])
                    expected_quads=np.array([(k,k+1,k+nx+2,k+nx+1) for j in range(ny) for i in range(nx) for k in [j*(nx+1)+i]])
                    assert np.array_equal(g['quads'],expected_quads)
                    assert np.max(abs(g['nodes']-expected_nodes))<1e-14
                    assert np.max(abs(g['centres']-expected_nodes[expected_quads].mean(axis=1)))<1e-14
                    mapping=common_mapping(paths[1],g)
                with h5py.File(paths[0]) as f:raw=chunk(f['results/1/phase-3/cells/SV_VOF'],n)
            else:
                paths=[run/'0/alpha.water',run/'0/C'];raw=field(paths[0],n)
                text=paths[1].read_text();m=re.search(r'internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\(',text);assert m and int(m[1])==n
                end=re.search(r'\n\s*\)\s*;',text[m.end():]);xy=np.fromstring(text[m.end():m.end()+end.start()].replace('(',' ').replace(')',' '),sep=' ').reshape(n,3)[:,:2]
                ij=np.rint(xy/dx-.5).astype(int);mapping=ij[:,0]+nx*ij[:,1];assert np.max(abs(xy-(ij+.5)*dx))<1e-10
            assert len(np.unique(mapping))==n
            for path in paths:
                if path.name=='geometry.npz' and not fine:continue # exact analytic/native structural proof above
                assert a['native_source_sha256'][str(path)]==sha(path)
            aligned=np.empty(n);aligned[mapping]=raw;h=aligned.reshape(ny,nx).sum(axis=0)*dx*1000
            if fine:h=h.reshape(80,2).mean(axis=1)
            rows.append(h);sources.extend(paths)
        difference=float(abs(rows[1]-rows[0]).max());assert difference<1e-9
        results[solver]=dict(actual_native_initial_columns_match=True,maximum_depth_difference_mm=difference,matched_strips=80,strip_width_mm=2.5,right_strip_bounds_mm=[197.5,200],right_strip_centre_mm=198.75)
    out=P/'matched-initial-strips-check-v001.json';assert not out.exists();out.write_text(json.dumps(dict(results=results,native_source_sha256={str(p):sha(p) for p in sources},moving_water_mesh_independence=False,physical_validation=False),indent=2)+'\n');print(json.dumps(results,indent=2))
if __name__=='__main__':main()
