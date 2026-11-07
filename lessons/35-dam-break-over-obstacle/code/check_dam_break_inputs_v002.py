"""Prove native mesh volumes/order and analytic initial fields before any solve."""
from pathlib import Path
import hashlib,json,re,numpy as np
from dam_break_common_v001 import *
P=Path('E:/CFD-Learning-Lab/local/dam-break-over-obstacle-v001');N=Path('E:/CFD-Learning-Lab/native/runs/dam-break-over-obstacle-v001')
def main():
    g=geometry();order=g['openfoam_to_common'];checks=[];sources={}
    for mode,rel in [('flat','qualification/openfoam-flat-v002'),('collapse','qualification/openfoam-collapse-v002'),('production','production/openfoam-v002')]:
        out=N/rel;mesh=(out/'mesh-check.log').read_text();assert 'Mesh OK.' in mesh and re.search(r'\bcells:\s+9072\b',mesh)
        c=native_field(out/'0/C',3);volume=native_field(out/'0/Vc',1).ravel()
        assert np.max(abs(c[:,:2]-g['centres'][order]))<1e-12 and np.max(abs(c[:,2]-DEPTH/2))<1e-12
        assert np.max(abs(volume/DEPTH-g['areas'][order]))<1e-15
        expected=initial_fields(g,mode);a=native_field(out/'0/alpha.water',1).ravel();pr=native_field(out/'0/p_rgh',1).ravel();u=native_field(out/'0/U',3)
        assert np.max(abs(a-expected['alpha'][order]))<1e-14 and np.max(abs(pr-expected['p_rgh'][order]))<1e-10 and np.max(abs(u))==0
        assert abs(a@volume/DEPTH-WATER_AREA)<1e-13
        checks.append(dict(mode=mode,native_mesh_and_topology_pass=True,native_centres_bijection_pass=True,native_volume_per_cell_pass=True,initial_fields_match=True,cells=len(a),water_area_m2=float(a@volume/DEPTH),native_fluid_area_m2=float(volume.sum()/DEPTH)))
        for name in ['0/C','0/Vc','0/alpha.water','0/p_rgh','0/U','system/blockMeshDict','mesh-check.log','declared-settings.json']:
            path=out/name;sources[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    record=dict(checks=checks,native_source_sha256=sources,flow_solvers_launched=False,physical_validation=False)
    out=P/'openfoam-native-input-check-v002.json';assert not out.exists();out.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(checks,indent=2))
if __name__=='__main__':main()