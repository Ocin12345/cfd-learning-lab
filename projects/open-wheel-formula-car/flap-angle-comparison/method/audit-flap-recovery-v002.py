"""Audit steady SST runs and render the solved vehicle surfaces and flow."""
import hashlib
import json
import re
import runpy
import struct
from pathlib import Path

import numpy as np

restart=runpy.run_path(str(Path(__file__).resolve().with_name('prepare-flap-restart-v002.py')))
configured=restart['configured']

NUMBER=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
PARTS=("chassis","frontWheels","rearWheels","frontMain","frontFlap","rearWing")
PATCHES=("inlet","outlet","road","farfield","centrePlane",*PARTS)
REF_FORCE=.5*1.225*30**2
REF_P=.5*30**2
EQS=("Ux","Uy","Uz","p","k","omega")


def audit(case,reader,iteration_cap,stages=()):
    read=runpy.run_path(str(reader))["foam_values"]
    hashes={}

    def raw(name, index=0):
        folder = (case,*stages)[index]
        data = (folder/name).read_bytes()
        key = f'{folder.name}/{name}' if index else name
        hashes[key] = hashlib.sha256(data).hexdigest()
        return data

    def native(name, index=0):
        return raw(name,index).decode("ascii")

    def field_stage(step):
        return int(step>600) + int(iteration_cap==1800 and step>1500)

    for manifest in ("pre-mesh-sha256.txt","input-sha256.txt"):
        for line in native(manifest).splitlines():
            sha,name=line.split("  ",1)
            raw(name)
            assert hashes[name]==sha,name
    for name in ("block-mesh.log","snappy.log","generator-sha256.txt","original-inputs-verified.log",
                 "inputs-verified.log","centres.log","volumes.log"):
        native(name)
    for part,expected in zip(PARTS,(1,2,2,1,1,1)):
        text=native(f"surface-{part}.log")
        assert "Surface is closed." in text
        assert re.search(rf"Number of unconnected parts\s*:\s*{expected}",text)
    stage_records=[]
    layout = (('blended',600,1200),) if iteration_cap==1200 else (('local',600,1500),('steady',1500,1800))
    for index,(phase,start,end) in enumerate(layout,1):
        if len(stages)<index:
            break
        folder=stages[index-1]
        receipt=json.loads(raw('restart-provenance.json',index))
        expected_source=(case,*stages)[index-1]
        assert (receipt['phase'],receipt['start'],receipt['end'])==(phase,start,end)
        assert receipt['parent']==case.name and receipt['checkpoint_case']==expected_source.name
        controls=('controlDict','fvSchemes','fvSolution')
        expected=configured({name:raw(f'system/{name}') for name in controls},phase)
        inputs=dict(line.split('  ',1)[::-1] for line in native('input-sha256.txt',index).splitlines())
        assert set(receipt['copied'])==set(inputs)
        for name,sha in receipt['copied'].items():
            from_checkpoint=name.startswith(f'{start}/')
            original=raw(name,index-1 if from_checkpoint else 0)
            assert hashlib.sha256(original).hexdigest()==sha,name
            updated=raw(name,index)
            assert hashlib.sha256(updated).hexdigest()==inputs[name],name
            changed=name.removeprefix('system/')
            assert updated==(expected[changed] if name.startswith('system/') and changed in expected else original),name
        assert 'Mesh OK.' in native('mesh-readback.log',index)
        assert native('inputs-verified.log',index).count(': OK')==2*len(inputs)
        stage_records.append({'phase':phase,'start':start,'end':end,'case':folder.name})
    solid_bounds={}
    for part in PARTS:
        stl=raw(f"constant/geometry/{part}.stl")
        count=struct.unpack_from("<I",stl,80)[0]
        records=list(struct.iter_unpack("<12fH",stl[84:]))
        assert len(records)==count and count>0
        xyz=np.array([row[3:12] for row in records]).reshape((-1,3))
        solid_bounds[part]={"triangles":count,"min_xyz_m":xyz.min(axis=0).tolist(),
                            "max_xyz_m":xyz.max(axis=0).tolist()}
    assert solid_bounds["frontMain"]["max_xyz_m"][0] < solid_bounds["frontFlap"]["min_xyz_m"][0]
    assert solid_bounds["frontFlap"]["max_xyz_m"][0] < solid_bounds["chassis"]["min_xyz_m"][0]
    assert solid_bounds["chassis"]["max_xyz_m"][1] < .79
    for part in PARTS:
        assert solid_bounds[part]["min_xyz_m"][2] > 0
    mesh=native("check-mesh.log")
    assert "Mesh OK." in mesh
    cells=int(re.search(r"^\s*cells:\s*(\d+)",mesh,re.M)[1])
    assert 50000<=cells<=130000
    residual=[]
    pressure_corrections=[]
    bounds=[(1,600),*( (record['start']+1,record['end']) for record in stage_records)]
    assert bounds[-1][1]==iteration_cap
    for index,(first,last) in enumerate(bounds):
        log=native('solver.log',index)
        assert log.rstrip().endswith('End') and 'kOmegaSST' in log and '14-7b05503f98a8' in log
        local = index>0 and stage_records[index-1]['phase']=='local'
        assert ('Using LTS' in log)==local
        times=list(re.finditer(r'^Time = (\d+)s?\s*$',log,re.M))
        assert [int(match[1]) for match in times]==list(range(first,last+1))
        for i,time in enumerate(times):
            stop=times[i+1].start() if i+1<len(times) else len(log)
            vals=re.findall(r'Solving for (Ux|Uy|Uz|p|k|omega), Initial residual = ('+NUMBER+r')',log[time.end():stop])
            series={key:[float(value) for field,value in vals if field==key] for key in EQS}
            count=2 if local else 1
            assert all(len(values)==(count if field=='p' else 1) for field,values in series.items()),(time[1],vals)
            residual.append([series[field][0] for field in EQS])
            if local:
                pressure_corrections.append(series['p'][1])
    residual=np.array(residual)
    assert np.isfinite(residual).all() and np.isfinite(pressure_corrections).all()

    def field(name,kind="scalar"):
        step=int(name.split('/',1)[0])
        a=read(native(name,field_stage(step)),"internalField",cells,kind)
        assert a.shape==((cells,3) if kind=="vector" else (cells,))
        assert np.isfinite(a).all(),name
        return a

    def body(text,name):
        start=re.search(r"\b"+name+r"\s*\{",text).end()
        depth=1
        for marker in re.finditer(r"[{}]",text[start:]):
            depth += 1 if marker[0]=="{" else -1
            if depth==0:
                return text[start:start+marker.start()]
        raise ValueError(f"Unclosed native patch {name}")

    centres,vol=field("0/C","vector"),field("0/Vc")
    assert (vol>0).all()
    faces=[np.fromstring(v,sep=" ",dtype=int) for v in re.findall(
        r"^\d+\(([0-9 ]+)\)$",native("constant/polyMesh/faces"),re.M)]
    owner=np.array(re.findall(r"^(\d+)$",native("constant/polyMesh/owner"),re.M),dtype=int)
    assert len(faces)==len(owner)-1==owner[0]
    owner=owner[1:]
    points=np.array(re.findall(r"^\(("+NUMBER+r") ("+NUMBER+r") ("+NUMBER+r")\)$",
                               native("constant/polyMesh/points"),re.M),dtype=float)
    boundary=native("constant/polyMesh/boundary")
    assert set(re.findall(r"^\s*(\w+)\s*\n\s*\{\s*type\s+",boundary,re.M))==set(PATCHES)
    patches={}
    for name in PATCHES:
        patch=body(boundary,name)
        count,start=(int(re.search(rf"\b{key}\s+(\d+);",patch)[1]) for key in ("nFaces","startFace"))
        assert count>0
        patches[name]=(start,count)
    assert points[:,1].min()>=-1e-12
    start,n=patches["centrePlane"]
    plane_vertices=np.unique(np.concatenate(faces[start:start+n]))
    assert np.max(abs(points[plane_vertices,1]))<1e-12
    assert re.search(r"type\s+symmetryPlane;",body(boundary,"centrePlane"))
    native_force={}
    for function in ("vehicleForce","chassisForce","frontWingForce","flapForce","rearWingForce"):
        merged=[]
        for index,(first,last) in enumerate(bounds):
            start=first-1
            reports=native(f'postProcessing/{function}/{start}/forces.dat',index)
            assert 'forces(pressure viscous)' in reports
            vals=np.array([re.findall(NUMBER,line) for line in reports.splitlines()
                           if line.strip() and not line.startswith('#')],dtype=float)
            assert vals.shape==(last-start+1,13) and np.array_equal(vals[:,0],np.arange(start,last+1))
            assert np.isfinite(vals).all()
            if index:
                assert np.allclose(merged[-1][-1,1:4],vals[0,1:4],rtol=1e-9,atol=1e-8),(function,start)
            merged.append(vals[1:] if index else vals)
        native_force[function]=np.vstack(merged)
    save_dirs=[]
    for index,(first,last) in enumerate(bounds):
        folder=(case,*stages)[index]
        native_dirs=sorted(int(p.name) for p in folder.iterdir() if p.is_dir() and p.name.isdigit())
        expected=list(range(first+49,last+1,50))
        assert native_dirs==([0,first-1,*expected] if index else [0,*expected])
        save_dirs.extend(expected)
    samples=[]
    for step in save_dirs:
        t=native(f"{step}/uniform/time",field_stage(step))
        assert int(re.search(r"^index\s+(\d+);",t,re.M)[1])==step
        assert float(re.search(r"^value\s+([^;]+);",t,re.M)[1])==step
        u,p=field(f"{step}/U","vector"),field(f"{step}/p")
        k,omega,nut=(field(f"{step}/{q}") for q in ("k","omega","nut"))
        assert (k>0).all() and (omega>0).all() and (nut>=0).all()
        flux=native(f"{step}/phi",field_stage(step))
        assert re.search(r"dimensions\s+\[0 3 -1 0 0 0 0\];",flux)
        rates={}
        wall=0.
        for name,(_,count) in patches.items():
            f=read(body(flux,name),"value",count,"scalar")
            rates[name]=float(f.sum())
            if name not in ("inlet","outlet"):
                wall=max(wall,float(np.max(abs(f))))
        assert rates["inlet"]<0<rates["outlet"]
        samples.append({"iteration":step,"minimum_pressure_m2_s2":float(p.min()),
                        "maximum_velocity_m_s":float(np.linalg.norm(u,axis=1).max()),
                        "minimum_tke_m2_s2":float(k.min()),"minimum_omega_s_inverse":float(omega.min()),
                        "maximum_nut_m2_s":float(nut.max()),"patch_volume_flux_m3_s":rates,
                        "relative_boundary_imbalance":abs(sum(rates.values()))/abs(rates["inlet"]),
                        "maximum_impermeable_face_flux_m3_s":wall})

    geom={}
    for name in PARTS:
        start,count=patches[name]
        poly=[points[ids] for ids in faces[start:start+count]]
        sf=np.array([.5*np.cross(v,np.roll(v,-1,axis=0)).sum(axis=0) for v in poly])
        norm=np.linalg.norm(sf,axis=1)
        c=np.array([v.mean(axis=0) for v in poly])
        own=owner[start:start+count]
        gap=abs(np.sum((c-centres[own])*sf,axis=1))/norm
        assert (norm>0).all() and (gap>0).all()
        geom[name]=(poly,c,sf,gap,p[own]/REF_P)

    u_text=native(f"{iteration_cap}/U",field_stage(iteration_cap))
    road=read(body(u_text,"road"),"value",patches["road"][1],"vector")
    assert np.max(abs(road-np.array([30.,0.,0.])))<=1e-12
    wheel_motion={}
    for name in ("frontWheels","rearWheels"):
        condition=body(u_text,name)
        assert "rotatingWallVelocity;" in condition
        _,cc,sf,_,_=geom[name]
        native_u=read(condition,"value",patches[name][1],"vector")
        normals=sf/np.linalg.norm(sf,axis=1)[:,None]
        lower,upper=cc[:,2]<.12,cc[:,2]>.64
        assert lower.sum()>=10 and upper.sum()>=10
        ux_lo,ux_hi=float(native_u[lower,0].mean()),float(native_u[upper,0].mean())
        normal_speed=float(np.max(abs(np.sum(native_u*normals,axis=1))))
        assert 18<ux_lo<36 and -36<ux_hi<-18 and normal_speed<1e-8
        wheel_motion[name]={"lower_mean_ux_m_s":ux_lo,"upper_mean_ux_m_s":ux_hi,
                            "maximum_normal_velocity_m_s":normal_speed}
    _,cc,sf,_,cp=geom["chassis"]
    area=np.linalg.norm(sf,axis=1)
    pressure_zones={}
    for name,xmin,xmax in (("flat_throat",2.6,3.17),("rear_floor",3.30,4.37)):
        mask=(cc[:,0]>=xmin)&(cc[:,0]<=xmax)&(abs(cc[:,1])<=.52)&(cc[:,2]<.47)&(sf[:,2]>.75*area)
        assert mask.sum()>=20
        pressure_zones[name]={"face_count":int(mask.sum()),"area_m2":float(area[mask].sum()),
                              "area_weighted_pressure_coefficient":float(np.average(cp[mask],weights=area[mask]))}
    patch_pressure={name:1.225*np.sum(a[2]*a[4][:,None]*REF_P,axis=0) for name,a in geom.items()}
    estimated_pressure=sum(patch_pressure.values())
    actual=native_force["vehicleForce"]
    half_total=actual[:,1:4]+actual[:,4:7]
    total=half_total*np.array([2.,0.,2.])
    chassis=native_force["chassisForce"]
    chassis_total=(chassis[:,1:4]+chassis[:,4:7])*np.array([2.,0.,2.])
    pressure_error=float(np.linalg.norm(estimated_pressure-actual[-1,1:4])
                         /max(np.linalg.norm(actual[-1,1:4]),1.))
    pressure_error_body=float(np.linalg.norm(patch_pressure["chassis"]-chassis[-1,1:4])
                              /max(np.linalg.norm(chassis[-1,1:4]),1.))
    assert pressure_error<1e-5 and pressure_error_body<1e-5,(pressure_error,pressure_error_body)
    groups={"chassisForce":("chassis",),"frontWingForce":("frontMain","frontFlap"),
            "flapForce":("frontFlap",),"rearWingForce":("rearWing",)}
    component={}
    for name,parts in groups.items():
        rows=native_force[name]
        p_force=sum(patch_pressure[part] for part in parts)
        discrepancy=float(np.linalg.norm(p_force-rows[-1,1:4])/max(1.,np.linalg.norm(p_force)))
        assert discrepancy<1e-5,(name,discrepancy)
        vector=rows[-50:,1:4]+rows[-50:,4:7]
        component[name]={"mean_drag_n":float(2*vector[:,0].mean()),
                         "mean_downforce_n":float(-2*vector[:,2].mean()),
                         "pressure_force_relative_readback":discrepancy}
    component["wheelsCombined"]={
        "mean_drag_n":float(total[-50:,0].mean()-sum(component[name]["mean_drag_n"]
                           for name in ("chassisForce","frontWingForce","rearWingForce"))),
        "mean_downforce_n":float(-total[-50:,2].mean()-sum(component[name]["mean_downforce_n"]
                                for name in ("chassisForce","frontWingForce","rearWingForce")))}
    tail=total[-50:,(0,2)]
    average=tail.mean(axis=0)
    span=np.ptp(tail,axis=0)
    force_tol=np.maximum(.03*abs(average),.0015*REF_FORCE)
    gates={
        "terminal_initial_residuals":bool(residual[-1].max()<=1e-3),
        "terminal_mass":samples[-1]["relative_boundary_imbalance"]<=1e-4,
        "impermeability":samples[-1]["maximum_impermeable_face_flux_m3_s"]<=1e-10,
        "positive_drag":bool(average[0]>0),
        "terminal_force_stability":bool((span<=force_tol).all()),
    }
    if iteration_cap in (1200,1800):
        gates['steady_final_50_all_initial_residuals']=bool(residual[-50:].max()<=1e-3)
    summary={"pass":all(gates.values()),"gates":gates,"case":case.name,"cells":cells,
             "iterations":iteration_cap,"saved_states":len(samples),"reference_force_n":REF_FORCE,
             "last_50_mean_drag_n":float(average[0]),"last_50_mean_downforce_n":float(-average[1]),
             "last_50_mean_drag_coefficient":float(average[0]/REF_FORCE),
             "last_50_mean_downforce_coefficient":float(-average[1]/REF_FORCE),
             "terminal_full_equivalent_force_n":total[-1].tolist(),
             "terminal_native_half_force_n":half_total[-1].tolist(),
             "zero_yaw_centre_plane_symmetry":True,
             "terminal_full_equivalent_chassis_force_n":chassis_total[-1].tolist(),
             "terminal_native_half_pressure_force_n":actual[-1,1:4].tolist(),
             "terminal_native_half_viscous_force_n":actual[-1,4:7].tolist(),
             "mean_pitch_moment_n_m_about_x2p17_z0p18":float(2*np.mean(actual[-50:,8]+actual[-50:,11])),
             "last_50_force_range_n":span.tolist(),"last_50_force_range_limits_n":force_tol.tolist(),
             "terminal_initial_residuals":dict(zip(EQS,residual[-1])),
             "rotating_wall_readback":wheel_motion,"underfloor_pressure_zones":pressure_zones,
             "component_forces":component,"source_geometry_bounds":solid_bounds,
             "surface_pressure_force_relative_readback":pressure_error,
             "surface_pressure_force_chassis_relative_readback":pressure_error_body,
             "chassis_cell_centre_wall_distance_quantiles_m":np.quantile(geom["chassis"][3],[0,.5,1]).tolist(),
             'final_50_maximum_initial_residuals':dict(zip(EQS,residual[-50:].max(axis=0))),
             'restart_stages':stage_records,
             "samples":samples,"source_sha256":hashes}
    return summary,geom,centres,u,p,total,residual,pressure_corrections
