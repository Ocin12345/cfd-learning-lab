"""Read-only rank-first audit of600actual thermal steps; preserves interrupted derived serial reconstruction honestly."""
from pathlib import Path
import argparse,csv,json,re,sys
import numpy as np
import prepare_openfoam_heated_channel_production_v005 as prep
import analyse_heated_channel_pilot_openfoam_v003 as thermal
from analyse_heated_channel_flow_openfoam_v002 import Evidence,field,cell_mapping
require=prep.require
sha=prep.sha
PATCHES=prep.PATCHES
PREPARATION_SHA='3b62ae10551d9cb3c62330469306215f82af706453400fbe2387a43f0ce9cd85'

def scalar_field(raw,count,components,dimensions):
    # Numerical nonfinite values remain values and false numerical flags, never parser substitutes.
    t=prep.clean(raw);require(re.search(r'\bformat\s+ascii\s*;',t),'Actual ASCII field required')
    dm=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',t);require(dm is not None and list(map(int,dm[1].split()))==dimensions,'Native field dimensions differ')
    return thermal.thermal_value(t,count,components,'internalField')

def runtime_flow(ur,pr,phr,patches,geometry):
    U=scalar_field(ur,prep.CELLS,3,[0,1,-1,0,0,0,0]);P=scalar_field(pr,prep.CELLS,1,[0,2,-2,0,0,0,0]);internal=scalar_field(phr,prep.INTERNAL_FACES,1,[0,3,-1,0,0,0,0])
    ub=thermal.thermal_blocks(ur);pb=thermal.thermal_blocks(pr);fb=thermal.thermal_blocks(phr);Phi=np.zeros(prep.TOTAL_FACES);Phi[:prep.INTERNAL_FACES]=internal;by={};us=[];ps=[];owner=np.asarray(geometry['owner'])
    for n in PATCHES:
        q=patches[n];own=owner[q['start']:q['start']+q['count']];count=q['count']
        uk='fixedValue' if n=='inlet' else 'zeroGradient' if n=='outlet' else 'noSlip';pk='fixedValue' if n=='outlet' else 'zeroGradient'
        require(re.search(r'\btype\s+'+uk+r'\s*;',ub[n]) and re.search(r'\btype\s+'+pk+r'\s*;',pb[n]),'Actual frozen flow boundary kinds differ')
        uv=thermal.thermal_value(ub[n],count,3) if re.search(r'\bvalue\s+',ub[n]) else np.zeros((count,3)) if n in ('top','bottom') else U[own]
        pv=thermal.thermal_value(pb[n],count,1) if re.search(r'\bvalue\s+',pb[n]) else P[own]
        v=thermal.thermal_value(fb[n],count,1);Phi[q['start']:q['start']+count]=v;by[n]=v;us.append(uv);ps.append(pv)
    require(re.search(r'\btype\s+empty\s*;',fb[prep.EMPTY]),'Actual empty phi constraint differs')
    return U,P,Phi,by,np.vstack(us),np.concatenate(ps)

class RankBoundaryReader:
    """Cache exact unchanged processor topology once; read every actual rank T/gradient thereafter."""
    def __init__(self,e,case,pilot,geometry,FG):
        self.e=e;self.FG=FG;self.maps=[];self.flowrefs={};self.face_maps={};seen=[];cells_seen=[];self.global_ids=np.concatenate([FG[n]['ids'] for n in PATCHES]);look={int(g):i for i,g in enumerate(self.global_ids)}
        gp=np.asarray(geometry['points']);gf=np.asarray(geometry['faces']);go=np.asarray(geometry['owner'])
        globalU,globalP,globalPhi,_,_,_=runtime_flow(e.text(case/'2000/U'),e.text(case/'2000/p'),e.text(case/'2000/phi'),prep.native_patches(prep.native_poly(case)),geometry)
        for rank in (0,1):
            poly=case/f'processor{rank}/0/polyMesh';old=pilot/f'processor{rank}/0/polyMesh'
            for n in ('points','faces','owner','neighbour','boundary','cellProcAddressing','faceProcAddressing','pointProcAddressing'):require(sha(e.touch(poly/n))==sha(e.touch(old/n)),'Inherited processor topology changed '+n)
            def labels(n):
                count,b=prep.counted(e.touch(poly/n));a=np.array(list(map(int,b.split())),dtype=np.int64);require(len(a)==count,'Processor label count differs');return a
            cells=labels('cellProcAddressing');addr=labels('faceProcAddressing');owner=labels('owner');pm=labels('pointProcAddressing');cells_seen.extend(cells.tolist())
            npnt,b=prep.counted(e.touch(poly/'points'));points=np.array([list(map(float,x.split())) for x in re.findall(r'\(([^()]*)\)',b)])
            require(points.shape==(npnt,3) and len(pm)==npnt and np.all((pm>=0)&(pm<prep.POINTS)) and np.max(abs(points-gp[pm]))<=1e-12,'Processor point coordinates/addressing differ')
            nf,b=prep.counted(e.touch(poly/'faces'));faces=np.array([list(map(int,x.split())) for n,x in re.findall(r'(\d+)\s*\(([^()]*)\)',b)])
            require(faces.shape==(nf,4) and len(addr)==len(owner)==nf,'Processor faces/address counts differ')
            nb,b=prep.counted(e.touch(poly/'boundary'));bp={n:dict(count=int(re.search(r'\bnFaces\s+(\d+)\s*;',q)[1]),start=int(re.search(r'\bstartFace\s+(\d+)\s*;',q)[1])) for n,q in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',b)}
            mappings=[]
            for n in PATCHES:
                q=bp[n];local=np.arange(q['start'],q['start']+q['count']);a=addr[local];ids=np.abs(a)-1
                require(np.all(a>0) and set(ids).issubset(set(FG[n]['ids'])),'Native physical face sign/patch differs');require(np.array_equal(go[ids],cells[owner[local]]),'Native physical owner/cell mapping differs')
                slots=np.array([look[int(fi)] for fi in ids],dtype=int);seen.extend(ids.tolist())
                for lf,fi in zip(local,ids):
                    require(set(pm[faces[lf]].tolist())==set(gf[fi].tolist()),'Processor/global face vertex set differs')
                    v=points[faces[lf]];sf=.5*np.sum(np.cross(v,np.roll(v,-1,axis=0)),axis=0);v=gp[gf[fi]];gsf=.5*np.sum(np.cross(v,np.roll(v,-1,axis=0)),axis=0)
                    require(np.max(abs(sf-gsf))<=1e-14,'Native physical outward orientation differs')
                mappings.append((n,len(local),slots))
            self.maps.append((rank,cells,mappings))
            ni=prep.counted(poly/'neighbour')[0];ref=self.rank_flow(case,rank,'2000',len(cells),nf,ni,bp)
            require(np.array_equal(ref[0],globalU[cells]) and np.array_equal(ref[1],globalP[cells]),'Actual decomposed frozen U/p do not map to global source')
            require(np.all(addr!=0) and np.all(abs(addr)<=prep.TOTAL_FACES) and np.max(abs(ref[2]-np.sign(addr)*globalPhi[abs(addr)-1]))<=1e-12,'Every native rank reference phi face/orientation differs from actual global frozen phi')
            self.flowrefs[rank]=(len(cells),nf,ni,bp,ref)
            self.face_maps[rank]=(addr,owner,bp,mappings);local_neighbour=labels('neighbour')
            require(len(np.unique(abs(addr)))==len(addr),'A processor maps a global face more than once')
            gn=np.asarray(geometry['neighbour'])
            for lf,a in enumerate(addr):
                fi=abs(int(a))-1;sign=1 if a>0 else -1
                require(set(pm[faces[lf]].tolist())==set(gf[fi].tolist()),'Every processor/global face vertex set must match')
                v=points[faces[lf]];sf=.5*np.sum(np.cross(v,np.roll(v,-1,axis=0)),axis=0);v=gp[gf[fi]];gsf=.5*np.sum(np.cross(v,np.roll(v,-1,axis=0)),axis=0)
                require(np.max(abs(sf-sign*gsf))<=1e-14,'Every native signed face orientation must match')
                expected_owner=go[fi] if sign>0 else gn[fi] if fi<len(gn) else -1
                require(cells[owner[lf]]==expected_owner,'Every native rank cell/face owner mapping must match')
                if lf<len(local_neighbour):require(cells[local_neighbour[lf]]==(gn[fi] if sign>0 else go[fi]),'Every native internal neighbour cell mapping must match')
        require(sorted(seen)==sorted(self.global_ids.tolist()) and len(set(seen))==240,'Each physical face must map exactly once');require(sorted(cells_seen)==list(range(prep.CELLS)),'Rank cells not bijective')
        self.report=dict(physical_faces=240,each_face_once=True,native_signed_face_addressing_positive=True,native_point_vertex_sets_owner_cells_and_outward_Sf_match=True,unchanged_pilot_production_processor_topology=True,rank_cell_counts=[len(x[1]) for x in self.maps])
    def assemble(self,case,label):
        T=np.full(prep.CELLS,np.nan);U=np.full((prep.CELLS,3),np.nan);P=np.full(prep.CELLS,np.nan);Phi=np.full(prep.TOTAL_FACES,np.nan);counts=np.zeros(prep.TOTAL_FACES,dtype=int)
        Tb=np.full(240,np.nan);Ub=np.full((240,3),np.nan);Pb=np.full(240,np.nan);duplicate_error=0.;rank_drift=np.zeros(3)
        for rank,cells,mappings in self.maps:
            nc,nf,ni,bp,ref=self.flowrefs[rank];addr,owner,_,_=self.face_maps[rank];folder=case/f'processor{rank}'/label
            rt=scalar_field(self.e.text(folder/'T'),nc,1,[0,0,0,1,0,0,0]);cur=self.rank_flow(case,rank,label,nc,nf,ni,bp)
            T[cells]=rt;U[cells]=cur[0];P[cells]=cur[1];rank_drift=np.maximum(rank_drift,[np.max(abs(x-y)) for x,y in zip(cur,ref)])
            ids=abs(addr)-1;oriented=np.sign(addr)*cur[2];existing=counts[ids]>0
            if np.any(existing):duplicate_error=max(duplicate_error,float(np.max(abs(Phi[ids[existing]]-oriented[existing]))))
            Phi[ids[~existing]]=oriented[~existing];counts[ids]+=1
            tb=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',prep.clean(self.e.text(folder/'T')).split('boundaryField',1)[1]));ub=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',prep.clean(self.e.text(folder/'U')).split('boundaryField',1)[1]));pb=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',prep.clean(self.e.text(folder/'p')).split('boundaryField',1)[1]))
            for n,count,slots in mappings:
                q=bp[n];own=owner[np.arange(q['start'],q['start']+q['count'])]
                tk='zeroGradient' if n=='outlet' else 'fixedValue';uk='fixedValue' if n=='inlet' else 'zeroGradient' if n=='outlet' else 'noSlip';pk='fixedValue' if n=='outlet' else 'zeroGradient'
                require(re.search(r'\btype\s+'+tk+r'\s*;',tb[n]) and re.search(r'\btype\s+'+uk+r'\s*;',ub[n]) and re.search(r'\btype\s+'+pk+r'\s*;',pb[n]),'Actual physical rank boundary kinds differ')
                tv=thermal.thermal_value(tb[n],count,1) if re.search(r'\bvalue\s+',tb[n]) else rt[own]
                if n!='outlet':require(np.all(tv==(300 if n=='inlet' else 350)),'Actual prescribed native rank temperature differs')
                uv=thermal.thermal_value(ub[n],count,3) if re.search(r'\bvalue\s+',ub[n]) else np.zeros((count,3)) if n in ('top','bottom') else cur[0][own]
                pv=thermal.thermal_value(pb[n],count,1) if re.search(r'\bvalue\s+',pb[n]) else cur[1][own]
                Tb[slots]=tv;Ub[slots]=uv;Pb[slots]=pv
        require(np.all((counts>=1)&(counts<=2)),'Actual native rank faces do not cover every global face')
        by={};tby={};offset=0
        for n in PATCHES:
            g=self.FG[n];count=len(g['ids']);by[n]=Phi[g['ids']];tby[n]=Tb[offset:offset+count];offset+=count
        return T,U,P,Phi,by,Ub,Pb,tby,rank_drift,duplicate_error

    def rank_flow(self,case,rank,label,ncells,nfaces,ninternal,bp):
        folder=case/f'processor{rank}'/label;U=scalar_field(self.e.text(folder/'U'),ncells,3,[0,1,-1,0,0,0,0]);P=scalar_field(self.e.text(folder/'p'),ncells,1,[0,2,-2,0,0,0,0]);raw=self.e.text(folder/'phi');Phi=np.zeros(nfaces);Phi[:ninternal]=scalar_field(raw,ninternal,1,[0,3,-1,0,0,0,0]);blocks=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',prep.clean(raw).split('boundaryField',1)[1]))
        for n,q in bp.items():
            block=blocks[n]
            if re.search(r'\btype\s+empty\s*;',block):continue
            Phi[q['start']:q['start']+q['count']]=thermal.thermal_value(block,q['count'],1)
        return U,P,Phi
    def frozen_drift(self,case,label):
        result=np.zeros(3)
        for rank,(nc,nf,ni,bp,ref) in self.flowrefs.items():
            cur=self.rank_flow(case,rank,label,nc,nf,ni,bp)
            result=np.maximum(result,np.array([np.max(abs(x-y)) for x,y in zip(cur,ref)]))
        return result
    def read(self,case,label,T):
        self.current_grad_internal_finite=True
        G=np.full((240,3),np.nan)
        for rank,cells,mappings in self.maps:
            folder=case/f'processor{rank}'/label;rt=scalar_field(self.e.text(folder/'T'),len(cells),1,[0,0,0,1,0,0,0]);require(np.array_equal(rt,T[cells],equal_nan=True),'Actual rank T differs from actual reconstructed T')
            raw=self.e.text(folder/'gradT');self.current_grad_internal_finite=self.current_grad_internal_finite and bool(np.isfinite(scalar_field(raw,len(cells),3,[0,-1,0,1,0,0,0])).all());blocks=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',prep.clean(raw).split('boundaryField',1)[1]))
            for n,count,slots in mappings:G[slots]=thermal.thermal_value(blocks[n],count,3)
        offset=0;out={}
        for n in PATCHES:c=len(self.FG[n]['ids']);out[n]=G[offset:offset+c];offset+=c
        return out

def time_metadata(e,folder,index):
    raw=e.text(folder/'uniform/time');d={n:float(re.search(r'\b'+n+r'\s+('+prep.FLOAT+r')\s*;',raw)[1]) for n in ('value','deltaT','index')}
    require(abs(d['value']-index*.005)<=1e-12 and abs(d['deltaT']-.005)<=1e-14 and d['index']==index,'Actual native time/dt/global index differs');return d

def residuals(e,case,indices):
    s=e.text(case/'solver.log');require(re.search(r'^End\s*$',s,re.M) and re.search(r'\bnProcs\s*:\s*2\b',s),'Completed actual two-rank solver log required');require(not re.search(r'FOAM FATAL|Solving for (?:Ux|Uy|Uz|p),',s),'Flow equation or fatal in thermal log')
    parts=re.split(r'^Time = ('+prep.FLOAT+r')s?\s*$',s,flags=re.M);require(len(parts)==1+2*len(indices),'Actual thermal step log count differs');rows=[];num='(?:'+prep.FLOAT+'|[-+]?(?:nan|inf))'
    for k,i in enumerate(indices):
        require(abs(float(parts[2*k+1])-i*.005)<=1e-12,'Actual thermal log physical time differs')
        rr=re.findall(r'Solving for T, Initial residual = ('+num+r'), Final residual = ('+num+r'), No Iterations (\d+)',parts[2*k+2],re.I);require(len(rr)==1,'Exactly one native T matrix solve each physical step required');rows.append([float(rr[0][0]),float(rr[0][1]),int(rr[0][2])])
    return np.asarray(rows)

def parsed_thermal_boundary(T,Tb_by,G_by,Phi_by,FG):
    arrays={n:[] for n in ('boundary_gradT','boundary_snGradT','boundary_T','reconstructed_boundary_snGradT','upwind_boundary_T','native_Qin_W_per_m','reconstructed_Qin_W_per_m','native_diffusive_Qin_W_per_m','native_convective_Qin_W_per_m')}
    for n in PATCHES:
        g=FG[n];Tb=Tb_by[n];G=G_by[n];Phi=Phi_by[n];own=g['owner'];sn=np.sum(G*g['normal'],axis=1);recon=np.zeros(len(Tb)) if n=='outlet' else (Tb-T[own])/g['distance'];up=T[own] if n=='outlet' else Tb
        qd=10*sn*g['area']/.01;qc=-1000*Phi*(up-300)/.01;qr=10*recon*g['area']/.01-1000*Phi*(np.where(Phi>=0,T[own],Tb)-300)/.01
        for key,value in [('boundary_gradT',G),('boundary_snGradT',sn),('boundary_T',Tb),('reconstructed_boundary_snGradT',recon),('upwind_boundary_T',up)]:arrays[key].append(value)
        for key,value in [('native_Qin_W_per_m',qd+qc),('reconstructed_Qin_W_per_m',qr),('native_diffusive_Qin_W_per_m',qd),('native_convective_Qin_W_per_m',qc)]:arrays[key].append(float(value.sum()))
    return {n:np.vstack(v) if n=='boundary_gradT' else np.concatenate(v) if n in ('boundary_snGradT','boundary_T','reconstructed_boundary_snGradT','upwind_boundary_T') else np.array(v) for n,v in arrays.items()}

def rank_precision_consistency(T,Tb,G,FG,C):
    ok=True
    for n,(axis,sign) in {'inlet':(0,-1),'outlet':(0,1),'top':(1,1),'bottom':(1,-1)}.items():
        g=FG[n];sn=G[n][:,axis]*sign;qn=thermal.significant_half_unit(G[n][:,axis],14)+8*np.spacing(np.maximum(abs(sn),np.finfo(float).tiny))
        if n=='outlet':lo=np.zeros(len(sn));hi=lo.copy()
        else:
            To=T[g['owner']];qT=thermal.significant_half_unit(To,14);Co=C[g['owner'],axis];Cf=g['centre'][:,axis];qd=thermal.significant_half_unit(Co,14)+8*abs(np.spacing(Cf))+8*abs(np.spacing(Co));d=g['distance'];require(np.all(d>qd),'Printed geometry interval includes zero distance')
            numerator=Tb[n]-To;v=np.stack([(numerator-qT)/(d-qd),(numerator-qT)/(d+qd),(numerator+qT)/(d-qd),(numerator+qT)/(d+qd)]);lo=v.min(axis=0);hi=v.max(axis=0)
        fp=8*np.finfo(float).eps*np.maximum(1,np.maximum(abs(lo),abs(hi)));ok=ok and bool(np.all((sn+qn>=lo-fp)&(sn-qn<=hi+fp)))
    return ok

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--authorization',type=Path,required=True);a=ap.parse_args();case=prep.checked_case(a.run);out=a.out.resolve();require(out.drive.upper()=='E:' and not out.exists() and case not in out.parents,'Fresh external E output required')
    require(sha(Path(prep.__file__))==PREPARATION_SHA,'Frozen production preparer changed');prep.verify_dependencies();e=Evidence();checks={}
    for p in (Path(__file__),Path(prep.__file__),*[Path(__file__).with_name(n) for n in prep.DEPENDENCIES]):e.touch(p)
    status=e.read_json(case/'run-status.json');require(status.get('stage') in ('reconstructing_actual_599_solved_states_preserving_rank_boundary_gradients','failed') and status.get('solver_started') is True,'Preserved interrupted derived-postprocess status required')
    if status.get('stage')=='failed':require(status.get('failed_stage')=='reconstructing_actual_599_solved_states_preserving_rank_boundary_gradients','Failure must concern derived reconstruction only')
    audit_auth=e.read_json(a.authorization);require(audit_auth.get('approved') is True and audit_auth.get('approved_by')=='root' and audit_auth.get('action')=='analyse_heated_channel_native_ranks_v006' and audit_auth.get('analysis_source_sha256')==sha(Path(__file__)) and prep.canonical_native_source_path(audit_auth['native_case'])==case and audit_auth.get('solver_completed') is True and audit_auth.get('postprocessing_stopped_by_root') is True,'Separate exact root rank-audit authorization required')
    receipt=prep.canonical_native_source_path(audit_auth['postprocess_receipt_path']);require(sha(e.touch(receipt))==audit_auth['postprocess_receipt_sha256'],'Root postprocess stop receipt differs');e.read_json(receipt)
    completion_proof=prep.private_root()/'root-native-rank-completion-proof-v001.json';require(sha(e.touch(completion_proof))=='4a81886e57c63275f2fbfc4574437b779bb06b50fd4edb9036990018b13acb36','Actual root native completion receipt changed');e.read_json(completion_proof)
    for rank in (0,1):
        physical=sorted(float(q.name) for q in (case/f'processor{rank}').iterdir() if q.is_dir() and re.fullmatch(prep.FLOAT,q.name) and 0<=float(q.name)<=3)
        require(len(physical)==601 and np.max(abs(np.asarray(physical)-np.arange(601)*.005))<=1e-12,'Every rank must retain exactly601actual physical states')
    dec=e.read_json(case/'declared-comparison.json');require(dec['thermal_start_time_s']==.005 and dec['thermal_end_time_s']==3 and dec['deltaT_s']==.005 and dec['additional_actual_steps_authorized']==599 and dec['total_actual_steps_authorized']==600,'Production cap/dt/continuation changed')
    for n,text in prep.production_dicts().items():require(e.text(case/n)==text,'Actual production dictionary differs '+n)
    for n,d in e.read_json(case/'provenance/prepared-input-sha256.json').items():require(sha(e.touch(case/n))==d,'Immutable prepared input changed '+n)
    require(sha(e.touch(case/'provenance/TARGETS.md'))==prep.TARGET_SHA==dec['target_document_sha256'],'Frozen targets changed')
    for folder in ('constant','system'):
        for n in ('fvModels','fvConstraints','fvOptions'):
            p=case/folder/n
            if p.exists():require(not re.sub(r'FoamFile\s*\{[^{}]*\}','',prep.clean(e.text(p))).strip(),'Unexpected external scalar models/sources')
    checks['actual_production_dictionaries_match_fixed_Euler_upwind_599_additional_targets']=True
    pilot=prep.run_root()/'openfoam-thermal-pilot-v003';pa=e.read_json(case/'provenance/accepted-pilot-checks.json');require(sha(case/'provenance/accepted-pilot-checks.json')==prep.PILOT_CHECKS_SHA and pa['all_checks_passed'] is True,'Accepted corrected pilot audit changed')
    for n,d in [('checks.json',prep.PILOT_CHECKS_SHA),('fields.npz',prep.PILOT_FIELDS_SHA),('source-manifest.json',prep.PILOT_MANIFEST_SHA)]:require(sha(e.touch(prep.private_root()/'openfoam-pilot-checks-v004'/n))==d,'Selected actual pilot audit changed')
    inherited=e.read_json(case/'provenance/inherited-pilot-inputs.json')
    for r in inherited['files']:require(sha(e.touch(prep.canonical_native_source_path(r['source'])))==r['sha256'],'Original inherited pilot input changed')
    for label in ('0','0.005'):require(sha(e.touch(case/label/'T'))==sha(e.touch(pilot/label/'T')),'Inherited actual pilot T changed')
    checks['accepted_actual_pilot_one_step_inputs_and_audit_preserved_byte_bound']=True
    refs=e.read_json(case/'provenance/frozen-flow-references.json');flow=prep.canonical_native_source_path(refs['actual_flow_case']);require(flow==prep.run_root()/'openfoam-flow-v002','Unexpected frozen flow source')
    require(sha(e.touch(case/'provenance/flow-checks.json'))==prep.FLOW_CHECKS_SHA and e.read_json(case/'provenance/flow-checks.json')['all_checks_passed'],'Original selected flow audit changed')
    for n,v in refs['fields'].items():require(sha(e.touch(prep.canonical_native_source_path(v['flow_path'])))==v['flow_sha256']==sha(e.touch(case/'2000'/n)),'Frozen flow actual source/copy changed')
    for n in ('U','p','phi'):require(sha(e.touch(case/'provenance/baseline-native-inputs'/n))==refs['fields'][n]['flow_sha256'],'Baseline native input changed')
    checks['actual_audited_flow2000_and_frozen_reference_copies_byte_bound_unchanged']=True
    poly=prep.native_poly(case);fp=prep.native_poly(flow)
    for n in ('points','faces','owner','neighbour','boundary'):require(sha(e.touch(poly/n))==sha(e.touch(fp/n)),'Production source mesh changed')
    patches=prep.native_patches(poly);geometry=prep.native_geometry(poly);C=field(e.text(case/'2000/C'),prep.CELLS,3,[0,1,0,0,0,0,0]);V=field(e.text(case/'2000/Vc'),prep.CELLS,1,[0,3,0,0,0,0,0]);require(np.all(V>0) and np.max(abs(V-np.asarray(geometry['volumes'])))<=1e-14,'Actual positive Vc/oriented volumes differ')
    mapping,mapping_record=cell_mapping(C,e.touch(case/'provenance/geometry.npz'),V,geometry['points']);FG=thermal.face_geometry(geometry,C,patches);ranks=RankBoundaryReader(e,case,pilot,geometry,FG)
    checks['actual_native_mesh_C_Vc_full_source_vertex_cell_mapping_positive']=True
    R=np.vstack((residuals(e,pilot,[1]),residuals(e,case,range(2,601))));checks['actual_two_rank_logs_600_physical_T_solves_without_flow_equations']=True
    U0,P0,Phi0,phiby0,Ub0,Pb0=runtime_flow(e.text(case/'2000/U'),e.text(case/'2000/p'),e.text(case/'2000/phi'),patches,geometry)
    histories={n:[] for n in ('T','U','p','phi','grad','sn','Tb','upwind','serialgrad','Q','Qr','Qd','Qc','freeze','bulk','precision','strict')};paths=[];finite=[];bounds=[];precision_reports=[];metadata=[]
    outlet_owner=FG['outlet']['owner'];weights=phiby0['outlet'];positive=weights>0
    require(np.any(positive),'Actual outward outlet convective flux needed')
    serial_checks=[];partial_serial=[];interface_errors=[];rank_meshes={}
    for kind,src in [('pilot',pilot),('production',case)]:
        rank_meshes[kind]=[{n:str(src/f'processor{r}/0/polyMesh'/n) for n in ('points','faces','owner','neighbour','boundary','cellProcAddressing','faceProcAddressing','pointProcAddressing')} for r in (0,1)]
    for i in range(601):
        label=f'{i*.005:.10g}';src=pilot if i<=1 else case;folder=src/label
        T,U,P,Phi,by,Ub,Pb,tby,rank_drift,interface_error=ranks.assemble(src,label);interface_errors.append(interface_error)
        if i:
            G=ranks.read(src,label,T);md=[]
            for r in (0,1):md.append(time_metadata(e,src/f'processor{r}'/label,i))
            metadata.append(md[0]);histories['precision'].append(rank_precision_consistency(T,tby,G,FG,C))
        else:
            raw=e.text(pilot/'0/gradT');scalar_field(raw,prep.CELLS,3,[0,-1,0,1,0,0,0]);blocks=thermal.thermal_blocks(raw);G={n:thermal.thermal_value(blocks[n],patches[n]['count'],3) for n in PATCHES}
        b=parsed_thermal_boundary(T,tby,G,by,FG);SG=np.full((240,3),np.nan)
        required=[folder/n for n in ('T','U','p','phi','gradT')]+([folder/'uniform/time'] if i else [])
        complete=all(q.is_file() for q in required)
        if complete and i>=2:complete=all(q.read_text(encoding='ascii').rstrip().endswith('// ************************************************************************* //') for q in required[:5])
        if complete:
            tr=e.text(folder/'T');ur=e.text(folder/'U');pr=e.text(folder/'p');phr=e.text(folder/'phi');gr=e.text(folder/'gradT');st=scalar_field(tr,prep.CELLS,1,[0,0,0,1,0,0,0]);su,sp,sphi,_,_,_=runtime_flow(ur,pr,phr,patches,geometry)
            serial_checks.append(dict(index=i,time_s=i*.005,T_path=str(folder/'T'),max_abs_T_K=float(np.max(abs(st-T))),max_abs_U_m_s=float(np.max(abs(su-U))),max_abs_p=float(np.max(abs(sp-P))),max_abs_phi_m3_s=float(np.max(abs(sphi-Phi)))))
            gb=thermal.thermal_blocks(gr);SG=np.vstack([thermal.thermal_value(gb[n],patches[n]['count'],3) for n in PATCHES])
        elif any(q.exists() for q in required):partial_serial.append(i)
        ok=bool(all(np.isfinite(q).all() for q in (T,U,P,Phi,Ub,Pb,rank_drift,*b.values())) and (not i or ranks.current_grad_internal_finite));finite.append(ok);bounds.append(bool(ok and np.all((T>=299.99)&(T<=350.01))))
        if i:histories['strict'].append(float(np.max(abs(b['boundary_snGradT']-b['reconstructed_boundary_snGradT'])/np.maximum(1,abs(b['boundary_snGradT'])))))
        for n,val in [('T',T),('U',U),('p',P),('phi',Phi),('grad',b['boundary_gradT']),('sn',b['boundary_snGradT']),('Tb',b['boundary_T']),('upwind',b['upwind_boundary_T']),('serialgrad',SG),('Q',b['native_Qin_W_per_m']),('Qr',b['reconstructed_Qin_W_per_m']),('Qd',b['native_diffusive_Qin_W_per_m']),('Qc',b['native_convective_Qin_W_per_m'])]:histories[n].append(val)
        histories['freeze'].append(np.maximum(rank_drift,[float(max(np.max(abs(U-U0)),np.max(abs(Ub-Ub0)))),float(max(np.max(abs(P-P0)),np.max(abs(Pb-Pb0)))),float(np.max(abs(Phi-Phi0)))]))
        w=by['outlet'];pos=w>0;histories['bulk'].append(float(np.sum(w[pos]*T[outlet_owner][pos])/w[pos].sum()) if np.any(pos) else float('nan'))
        paths.append(dict(index=i,time_s=i*.005,native_case=str(src),time_label=label,source_kind='native_rank_reassembly',rank_T=[str(src/f'processor{r}'/label/'T') for r in (0,1)],rank_U=[str(src/f'processor{r}'/label/'U') for r in (0,1)],rank_p=[str(src/f'processor{r}'/label/'p') for r in (0,1)],rank_phi=[str(src/f'processor{r}'/label/'phi') for r in (0,1)],rank_gradT=[str(src/f'processor{r}'/label/'gradT') for r in (0,1)] if i else [],rank_time_metadata=[str(src/f'processor{r}'/label/'uniform/time') for r in (0,1)] if i else [],rank_mesh=rank_meshes['pilot' if i<=1 else 'production'],baseline_gradT=str(pilot/'0/gradT') if not i else None))
        if i%100==0:print('Read actual rank thermal state '+str(i)+'/600',flush=True)
    H={n:np.asarray(v) for n,v in histories.items()};E=1000*np.sum(V[None,:]*(H['T']-300),axis=1)/.01;storage=np.diff(E)/.005;Q=H['Q'][1:];Qr=H['Qr'][1:];den=np.maximum.reduce((np.sum(abs(Q),axis=1),abs(storage),np.full(600,5000.)));defect=abs(storage-Q.sum(axis=1))/den;pde=abs(storage-Qr.sum(axis=1))/den
    checks['all601_actual_native_T_flow_gradient_values_finite']=bool(np.all(finite));checks['all601_T_bounds_300_350_plus0_01_K']=bool(np.all(bounds));checks['actual_initial_T300_and_first_solved_T_changes']=bool(np.all(H['T'][0]==300) and np.max(abs(H['T'][1]-H['T'][0]))>0)
    checks['all600_final_linear_T_residuals_le1e_8']=bool(np.isfinite(R).all() and np.all(R[:,1]<=1e-8));checks['all601_runtime_internal_boundary_U_p_allface_phi_frozen_le1e_12']=bool(np.isfinite(H['freeze']).all() and np.max(H['freeze'])<=1e-12)
    checks['actual_inlet_inward_outlet_outward_walls_zero_native_phi']=bool(np.all(phiby0['inlet']<0) and np.all(phiby0['outlet']>=0) and sum(abs(phiby0[n]).sum() for n in ('top','bottom'))<=1e-14)
    serial_summary=dict(indices=[r['index'] for r in serial_checks],count=len(serial_checks),**{n:max(r[n] for r in serial_checks) for n in ('max_abs_T_K','max_abs_U_m_s','max_abs_p','max_abs_phi_m3_s')})
    require(0 in serial_summary['indices'] and 1 in serial_summary['indices'],'Authoritative original pilot initial/first-step serial validation required')
    checks['all_available_closed_serial_states_equal_native_rank_reassembly']=bool(all(np.isfinite(serial_summary[n]) and serial_summary[n]<=1e-12 for n in ('max_abs_T_K','max_abs_U_m_s','max_abs_p','max_abs_phi_m3_s')))
    checks['all601_native_rank_interface_phi_values_agree_le1e_12']=bool(np.isfinite(interface_errors).all() and max(interface_errors)<=1e-12)
    checks['all600_actual_rank_boundary_gradT_face_owner_vertex_orientation_T_mapping_proved']=True
    checks['all600_global_native_Euler_relative_energy_defect_le1e_5']=bool(np.all(finite) and np.isfinite(defect).all() and np.all(defect<=1e-5));checks['all600_separate_reconstructed_PDE_relative_defect_le0_005']=bool(np.all(finite) and np.isfinite(pde).all() and np.all(pde<=.005))
    gauss=e.text(case/'provenance/installed-sources/src/finiteVolume/finiteVolume/gradSchemes/gaussGrad/gaussGrad.C');require('vsf.boundaryField()[patchi].snGrad()' in gauss and '(n & gGradbf[patchi])' in gauss,'Installed native boundary normal correction proof differs')
    scalar=e.text(case/'provenance/installed-sources/src/functionObjects/solvers/scalarTransport/scalarTransport.H');require(re.search(r'bool\s+executeAtStart\(\)\s+const\s*\{\s*return false;',scalar),'Native scalar must not duplicate startup T solve')
    for d,p in (q.split(None,1) for q in e.text(case/'provenance/installed-source-sha256.txt').splitlines() if q.strip()):require(sha(e.touch(case/'provenance/installed-sources'/p[len('/opt/openfoam14/'):]))==d,'Installed source digest differs')
    verification=e.text(case/'installed-source-verification.log');require(verification.count(': OK')==len(prep.INSTALLED_SOURCES)+1 and 'FAILED' not in verification,'Actual current native installed-source identity check missing')
    matrix=e.text(case/'provenance/installed-fvMatrixSolve.C');require('word(name + "Final")' in matrix and 'solutionControl::finalIteration' in matrix,'Native TFinal selection changed')
    checks['installed_native_functions_scalar_startup_grad_IO_sources_hash_bound']=True
    auth=e.read_json(case/'provenance/root-run-authorization.json');require(auth.get('approved') is True and auth.get('approved_by')=='root' and auth.get('action')=='run_heated_channel_openfoam_production_v005' and auth.get('allowed_mode')=='ThermalProduction','Separate root production admission missing')
    for n,k in [('prepare_openfoam_heated_channel_production_v005.py','preparation_sha256'),('run-openfoam-thermal-production-v005.sh','wrapper_sha256'),('Launch-OpenFOAM-ThermalProduction-v005.ps1','launcher_sha256')]:require(sha(e.touch(case/'provenance'/n))==auth[k],'Frozen reviewed execution source differs '+n)
    require(auth['preparation_sha256']==PREPARATION_SHA and sha(e.touch(case/'provenance/run-openfoam-executed.sh'))==auth['wrapper_sha256'],'Executed wrapper/preparer differs');checks['separate_root_599_step_production_admission_exact_source_hashes_bound']=True
    source=e.finish();checks['all_consumed_native_source_inputs_unchanged_during_read_only_analysis']=True;checks={n:bool(v) for n,v in checks.items()}
    mapping['source_cell_to_native_cell']=np.argsort(mapping['native_cell_to_source_cell']);mapping['source_ix']=mapping['ix'][mapping['source_cell_to_native_cell']];mapping['source_iy']=mapping['iy'][mapping['source_cell_to_native_cell']]
    record=dict(phase='600_step_thermal_production',native_case=str(case),native_state_paths=paths,physical_times_s=(np.arange(601)*.005).tolist(),actual_steps=600,inherited_actual_steps=1,additional_actual_steps=599,mpi_ranks=2,completed_600_step_transient=True,
      execution_completion=dict(actual_solver_completed=True,inherited_actual_steps=1,additional_actual_steps=599,total_actual_steps=600,solver_log_End=True,all_rank_states_present=True,raw_wrapper_stage=status['stage'],raw_wrapper_status=status,serial_reconstruction_complete=False,root_interrupted_postprocess=True,wrapper_finished_normally=False,root_postprocess_receipt_path=str(receipt),root_native_completion_proof_path=str(completion_proof)),
      native_rank_reassembly_proof=dict(all601_states_reassembled_from_actual_native_ranks=True,all600_both_rank_time_dt_global_indices_proved=True,global_cells=prep.CELLS,rank_cell_counts=ranks.report['rank_cell_counts'],full_native_cell_face_point_addressing_mapping_proved=True,available_closed_serial_state_crosschecks=serial_summary,original_pilot_initial_and_first_step_crosschecks_passed=checks['all_available_closed_serial_states_equal_native_rank_reassembly'],incomplete_partial_serial_indices=partial_serial,native_rank_phi_interface_max_abs_difference=float(max(interface_errors))),
      actual_settings=dec,actual_cell_mapping=mapping_record,actual_native_patches=patches,flow_source=refs,accepted_pilot_checks_sha256=prep.PILOT_CHECKS_SHA,retained_rank_native_gradient_mapping=ranks.report,
      final_linear_T_residual_max=float(np.max(R[:,1])),initial_linear_T_residual_max=float(np.max(R[:,0])),native_T_min_K=float(np.min(H['T'])),native_T_max_K=float(np.max(H['T'])),
      frozen_fields=dict(max_abs_change_U_m_s=float(np.max(H['freeze'][:,0])),max_abs_change_kinematic_p=float(np.max(H['freeze'][:,1])),max_abs_change_all_native_phi_m3_s=float(np.max(H['freeze'][:,2])),all601_actual_states_checked=True,internal_and_physical_boundary_values_checked=True,all_native_processor_face_phi_including_both_interface_sides_checked=True,rank_reference_phi_signed_face_addressing_matches_actual_global_source=True,empty_phi_zero_by_native_empty_constraint=True),
      energy=dict(boundary_order=list(PATCHES),reference_T_K=300,depth_m=.01,final_E_J_per_m=float(E[-1]),native_Euler_relative_defect_max=float(np.max(defect)),reconstructed_PDE_relative_defect_max=float(np.max(pde)),all_cell_energy_balance_claimed=False),
      outlet_outward_flux_weighted_bulk_T_final_K=float(H['bulk'][-1]),native_time_metadata_all600=metadata,
      diagnostic_comparisons=dict(original_strict_stencil_relative_target=1e-10,original_strict_stencil_comparison_passed=bool(np.isfinite(H['strict']).all() and np.all(H['strict']<=1e-10)),original_strict_stencil_comparison_all600_passed=bool(np.isfinite(H['strict']).all() and np.all(H['strict']<=1e-10)),original_strict_stencil_max_relative_error=float(np.max(H['strict'])),original_strict_failed_step_indices=(np.where(~(H['strict']<=1e-10))[0]+1).tolist(),
        print_precision_significant_digits=14,precision_interval_all600_consistent=bool(np.all(H['precision'])),precision_interval_inconsistent_step_indices=(np.where(~H['precision'].astype(bool))[0]+1).tolist(),classification='Auxiliary instrumentation diagnostics, distinct from unchanged physical/numerical targets; strict1e-10 preserved without waiver'),
      reconstructed_serial_gradT_caveat='extrapolatedCalculated field reload resets saved physical gradients to owner internal gradients; direct heat flux uses literal retained rank boundary gradients with actual native addressing',
      checks=checks,check_results=checks,failed_checks=[n for n,v in checks.items() if not v],all_checks_passed=all(checks.values()),all_declared_checks_passed=all(checks.values()),analysis_source_sha256=sha(Path(__file__)),preparation_source_sha256=PREPARATION_SHA,physical_validation=False,
      interpretation='Actual native rank fields reassembled through exact cell/signed face/point addressing. All600 physical steps solved; derived serial postprocessing interrupted by root. No synthetic global native files or transient momentum/physical-validation claim.')
    out.mkdir(parents=True)
    arrays=dict(time_s=np.arange(601)*.005,solved_time_s=np.arange(1,601)*.005,T_history=H['T'],T=H['T'][-1],initial_T=H['T'][0],centres=C,volumes=V,areas=V/.01,grid_shape=np.array([24,96]),**mapping,
      U_history=H['U'],p_history=H['p'],native_phi_history=H['phi'],frozen_reference_U=U0,frozen_reference_p=P0,frozen_reference_native_phi=Phi0,
      physical_boundary_face_indices=ranks.global_ids,boundary_names=np.array(PATCHES),boundary_offsets=np.array([0,24,48,144,240]),native_boundary_gradT_history=H['grad'],native_boundary_snGradT_history=H['sn'],native_boundary_T_history=H['Tb'],native_upwind_boundary_T_history=H['upwind'],reconstructed_serial_boundary_gradT_history=H['serialgrad'],serial_diagnostic_available=np.array([i in serial_summary['indices'] for i in range(601)]),
      energy_above300_J_per_m=E,storage_W_per_m=storage,native_Qin_W_per_m=Q,reconstructed_Qin_W_per_m=Qr,native_Euler_relative_defect=defect,reconstructed_PDE_relative_defect=pde,native_diffusive_Qin_W_per_m=H['Qd'][1:],native_convective_Qin_W_per_m=H['Qc'][1:],boundary_power_reconstructed_minus_native_W_per_m=Qr-Q,outlet_bulk_T_K=H['bulk'],
      initial_linear_T_residual=R[:,0],final_linear_T_residual=R[:,1],linear_T_iterations=R[:,2].astype(int),frozen_U_max_abs_change_m_s=H['freeze'][:,0],frozen_p_max_abs_change=H['freeze'][:,1],frozen_native_phi_max_abs_change_m3_s=H['freeze'][:,2],strict_stencil_max_relative_error=H['strict'],precision_interval_consistent=H['precision'].astype(bool),rho_kg_m3=np.array(1.),cp_J_kg_K=np.array(1000.),k_W_m_K=np.array(10.),depth_m=np.array(.01))
    np.savez_compressed(out/'fields.npz',**arrays);record['audited_fields_npz_sha256']=sha(out/'fields.npz')
    with (out/'thermal-step-history.csv').open('x',encoding='ascii',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','dt_s','E_J_per_m','storage_W_per_m','Qin_inlet_W_per_m','Qin_outlet_W_per_m','Qin_top_W_per_m','Qin_bottom_W_per_m','native_Euler_relative_defect','reconstructed_PDE_relative_defect','initial_linear_T_residual','final_linear_T_residual'])
        for j in range(600):w.writerow([(j+1)*.005,.005,E[j+1],storage[j],*Q[j],defect[j],pde[j],R[j,0],R[j,1]])
    record['thermal_history_csv_sha256']=sha(out/'thermal-step-history.csv');source['audited_outputs']={n:dict(sha256=sha(out/n),bytes=(out/n).stat().st_size) for n in ('fields.npz','thermal-step-history.csv')};prep.json_new(out/'checks.json',thermal.finite_json(record));prep.json_new(out/'source-manifest.json',source)
    for p in (Path(__file__),Path(prep.__file__),*[Path(__file__).with_name(n) for n in prep.DEPENDENCIES]):prep.copied(p,out/('executed-'+p.name))
    print(json.dumps(dict(output=str(out),all_checks_passed=all(checks.values()),failed_checks=record['failed_checks']),indent=2));return 0 if all(checks.values()) else 2
if __name__=='__main__':sys.exit(main())
