"""Diagnostic reader: continuity with zero scaling factor is unscaled, verified against native console residuals. Other zero denominators remain rejected."""
import numpy as np
import h5py
from plot_nozzle_native_v004 import require,indexed,TOL
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
            require(np.isfinite(iteration) and iteration>0 and iteration==int(iteration) and len(row)>=2 and np.isfinite(row[:2]).all() and row[1]>=0 and (row[1]>0 or name=='continuity'),'Invalid actual Fluent residual evidence')
            residuals[name]=float(row[0]/row[1]) if row[1]>0 else float(row[0]); terminals.append(int(iteration))
        require(len(residuals)>=3 and len(set(terminals))==1,'Incomplete actual Fluent solve evidence')
    result=dict(nodes=nodes,cells=cells,U=U,p=pressure,T=temp,layer='native planar2D quads',solve_evidence=dict(actual_terminal_iteration=terminals[0],terminal_scaled_residuals=residuals,operating_pressure_Pa=operating_pressure,initialized_only_data_refused=True))
    if rho is not None: result['rho']=rho
    return result
