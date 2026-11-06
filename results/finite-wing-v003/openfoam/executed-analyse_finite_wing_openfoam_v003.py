"""Read-only native checks for the saved 3D finite-wing Foundation14 run.

--run names the actual OpenFOAM case, --out a fresh external analysis folder.
Exit2 retains genuine missed numerical targets. No CFD or WSL is started.
"""
from pathlib import Path
import argparse
import csv
import json
import re
import sys
import numpy as np
from prepare_openfoam_finite_wing_v003 import (
    PATCHES, FLOAT, require, sha, clean, dictionaries, native_poly, native_patches,
    mesh_quality, json_new, copied,
)


class Evidence:
    def __init__(self):
        self.before = {}

    def text(self, path):
        path = Path(path).resolve()
        require(path.is_file(), 'Missing actual native evidence: ' + str(path))
        if path not in self.before:
            self.before[path] = sha(path)
        return path.read_text(encoding='utf-8-sig')

    def read_json(self, path):
        return json.loads(self.text(path))

    def finish(self):
        rows = [dict(path=str(path),bytes=path.stat().st_size,sha256_before=before,sha256_after=sha(path))
                for path,before in self.before.items()]
        require(all(r['sha256_before']==r['sha256_after'] for r in rows), 'Native source changed during read-only analysis')
        return dict(inputs_unchanged_during_analysis=True,sources=rows,solver_started_or_modified=False)


def values(raw, count, components):
    text = clean(raw)
    require(re.search(r'\bformat\s+ascii\s*;',text), 'Native ASCII field required')
    require(not re.search(r'(?i)(?<![A-Za-z0-9_])[-+]?(?:nan|inf(?:inity)?)(?![A-Za-z0-9_])',text), 'Nonfinite native field/boundary token')
    m = re.search(r'\binternalField\s+uniform\s+([^;]+);',text)
    if m:
        row = np.fromstring(m[1].replace('(',' ').replace(')',' '),sep=' ')
        require(row.size==components and np.isfinite(row).all(),'Invalid uniform native field')
        return np.tile(row,(count,1)) if components>1 else np.full(count,row[0])
    kind = 'vector' if components==3 else 'scalar'
    m = re.search(r'\binternalField\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
    require(m is not None and int(m[1])==count, 'Native internal-field count differs')
    data = np.fromstring(m[2].replace('(',' ').replace(')',' '),sep=' ')
    require(data.size==count*components and np.isfinite(data).all(),'Native field length/finite check failed')
    return data.reshape((count,components)) if components>1 else data


def field(raw,count,components,dimensions):
    d=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',clean(raw))
    require(d is not None and list(map(int,d[1].split()))==dimensions,'Native field dimensions differ')
    return values(raw,count,components)


def boundary_mass(raw,patches):
    text=clean(raw)
    require('boundaryField' in text,'Missing phi boundary fields')
    blocks=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',text.split('boundaryField',1)[1]))
    totals={}
    inward=0.0
    outward=0.0
    for name in PATCHES:
        require(name in blocks,'Missing native phi patch '+name)
        body=blocks[name]
        m=re.search(r'\bvalue\s+uniform\s+('+FLOAT+r')\s*;',body)
        if m:
            data=np.full(patches[name]['count'],float(m[1]))
        else:
            m=re.search(r'\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;',body,re.S)
            require(m is not None and int(m[1])==patches[name]['count'],'Native phi patch count differs '+name)
            data=np.fromstring(m[2],sep=' ')
        require(data.size==patches[name]['count'] and np.isfinite(data).all(),'Nonfinite/wrong phi boundary length')
        totals[name]=float(data.sum())
        inward-=float(data[data<0].sum())
        outward+=float(data[data>0].sum())
    require(inward>0,'Missing native inward boundary flux')
    net=sum(totals.values())
    return dict(signed_outward_volume_flux_m3_s=totals,total_inward_volume_flux_m3_s=inward,
        total_outward_volume_flux_m3_s=outward,net_outward_volume_flux_m3_s=net,
        relative_imbalance=abs(net)/inward,
        definition='abs(sum outward phi over every physical boundary face)/sum magnitude of inward face phi; includes farfield and wing')


def residual_history(text):
    require(not re.search(r'FOAM FATAL|Segmentation fault|Floating point exception\s*$',text,re.M),'Actual solver fatal error')
    sections=re.split(r'^Time = ('+FLOAT+r')s?\s*$',text,flags=re.M)
    rows=[]
    for i in range(1,len(sections),2):
        iteration=float(sections[i])
        samples={name:[] for name in ('p','Ux','Uy','Uz')}
        for name,initial in re.findall(r'Solving for (p|Ux|Uy|Uz), Initial residual = ('+FLOAT+r')',sections[i+1]):
            samples[name].append(float(initial))
        require(all(samples.values()),'Missing actual p/Ux/Uy/Uz residual at iteration '+str(iteration))
        row=dict(iteration=int(iteration),**{name:max(v) for name,v in samples.items()})
        require(iteration==row['iteration'] and all(np.isfinite(v) for v in row.values()),'Invalid steady residual samples')
        rows.append(row)
    require([r['iteration'] for r in rows]==list(range(1,1201)),'Need actual fixed-cap iterations1..1200')
    return rows


def forces(text):
    header=re.search(r'^#\s*Time\s+Cm\s+Cd\s+Cl\s+Cl\(f\)\s+Cl\(r\)\s*$',text,re.M)
    require(header is not None,'Native forceCoeffs columns differ from Foundation14 official source')
    rows=[]
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        data=list(map(float,line.split()))
        require(len(data)==6 and all(np.isfinite(data)),'Invalid native force row')
        iteration=data[0]
        require(iteration.is_integer(),'Steady force index must be an iteration, not physical time')
        if iteration>0:
            rows.append(dict(iteration=int(iteration),Cm=data[1],Cd=data[2],Cl=data[3]))
    indices=[r['iteration'] for r in rows]
    require(indices==sorted(set(indices)) and all(1<=i<=1200 for i in indices),
            'Native force iteration indices must be unique, increasing, within1..1200')
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    case,out=args.run.resolve(),args.out.resolve()
    require(not out.exists(),'Analysis output must be fresh')
    require(out.drive.upper()=='E:','Analysis outputs belong on E:')
    e=Evidence()
    for source in (Path(__file__),Path(__file__).with_name('prepare_openfoam_finite_wing_v003.py')):
        e.text(source)
    declaration=e.read_json(case/'declared-comparison.json')
    status=e.read_json(case/'run-status.json')
    require(status.get('stage')=='finished_execution' and status.get('final_iteration')==1200,'Native fixed-cap execution is incomplete')
    checks={}
    for name,expected in dictionaries().items():
        require(e.text(case/name)==expected,'Actual native dictionary differs: '+name)
    checks['actual_native_dictionaries_match_exact_recorded_settings']=True
    prepared=e.read_json(case/'provenance/prepared-input-sha256.json')
    for name,digest in prepared.items():
        e.text(case/name)
        require(sha(case/name)==digest,'Prepared source changed: '+name)
    checks['prepared_input_hashes_match']=True
    mesh=case/'input-wing.msh'
    require(sha(mesh)==declaration['mesh_source_sha256'],'Shared source mesh hash differs')
    e.text(mesh)
    checks['shared_source_mesh_hash_matches']=True
    preflight=e.read_json(case/'provenance/preflight-validation.json')
    require(preflight.get('passed') is True,'Native preflight did not pass')
    poly=native_poly(case)
    for name,digest in preflight['native_mesh_sha256'].items():
        e.text(poly/name)
        require(sha(poly/name)==digest,'Native static mesh changed')
    patches=native_patches(poly)
    quality=mesh_quality(e.text(case/'mesh-check.log'))
    require(quality['cells']==declaration['intended_cells'],'Native cell count differs')
    checks['native_positive_volume_nonorth_skew_and_full_mesh_check_passed']=True
    cells=quality['cells']
    internal_faces=min(p['start'] for p in patches.values())
    stats={}
    for iteration in (300,600,900,1200):
        u=field(e.text(case/str(iteration)/'U'),cells,3,[0,1,-1,0,0,0,0])
        p=field(e.text(case/str(iteration)/'p'),cells,1,[0,2,-2,0,0,0,0])
        phi=field(e.text(case/str(iteration)/'phi'),internal_faces,1,[0,3,-1,0,0,0,0])
        stats[str(iteration)]=dict(U_min_m_s=u.min(axis=0).tolist(),U_max_m_s=u.max(axis=0).tolist(),
            speed_min_m_s=float(np.linalg.norm(u,axis=1).min()),speed_max_m_s=float(np.linalg.norm(u,axis=1).max()),
            gauge_pressure_Pa_min=float(p.min()),gauge_pressure_Pa_max=float(p.max()),internal_phi_finite=True)
    centres=field(e.text(case/'1200/C'),cells,3,[0,1,0,0,0,0,0])
    require(all(np.ptp(centres[:,j])>1 for j in range(3)),'Native cell centres do not span a full 3D domain')
    checks['all_saved_native_fields_finite_with_correct_counts_and_dimensions']=True
    checks['native_solution_is_three_dimensional_domain']=True
    residual_rows=residual_history(e.text(case/'solver.log'))
    terminal={name:residual_rows[-1][name] for name in ('p','Ux','Uy','Uz')}
    checks['final_iteration_max_initial_p_Ux_Uy_Uz_residual_le_1e_5']=max(terminal.values())<=1e-5
    mass=boundary_mass(e.text(case/'1200/phi'),patches)
    checks['native_upstream_inflow_and_outlet_outflow']=mass['signed_outward_volume_flux_m3_s']['inlet']<0<mass['signed_outward_volume_flux_m3_s']['outlet']
    checks['relative_all_physical_boundary_mass_imbalance_le_1e_4']=mass['relative_imbalance']<=1e-4
    raw_forces=sorted((case/'postProcessing/forceCoefficients').glob('*/forceCoeffs.dat'))
    require(len(raw_forces)==1,'Fresh fixed-cap run must have one actual force history')
    force_rows=forces(e.text(raw_forces[0]))
    stability={}
    for name in ('Cd','Cl'):
        window=force_rows[-100:]
        data=np.array([row[name] for row in window])
        complete=[row['iteration'] for row in window]==list(range(1101,1201))
        mean=float(data.mean()) if len(data) else None
        span=float(np.ptp(data)) if len(data) else None
        relative=span/max(abs(mean),1e-3) if len(data) else None
        stability[name]=dict(actual_samples=len(data),final100_iteration_window_complete=complete,
            first_iteration=window[0]['iteration'] if window else None,
            last_iteration=window[-1]['iteration'] if window else None,mean=mean,
            minimum=float(data.min()) if len(data) else None,maximum=float(data.max()) if len(data) else None,
            absolute_span=span,relative_span=relative,denominator='max(abs(mean),1e-3)')
        checks['last100_'+name+'_relative_span_le_1e_3']=complete and relative<=1e-3
    checks['actual_force_samples_complete']=[row['iteration'] for row in force_rows]==list(range(1,1201))
    auth=e.read_json(case/'provenance/root-authorization.json')
    require(auth.get('approved') is True and auth.get('approved_by')=='root' and auth.get('action')=='run_finite_wing_openfoam_v003','Missing run authorization evidence')
    executed=case/'provenance/run-openfoam-executed.sh'
    wrapper=case/'provenance/run-openfoam.sh'
    e.text(executed)
    require(sha(executed)==sha(wrapper),'Executed wrapper differs from prepared reviewed wrapper')
    checks['executed_wrapper_matches_prepared_source']=True
    source=e.finish()
    checks['inputs_unchanged_during_read_only_analysis']=True
    record=dict(solver='OpenFOAM Foundation14 incompressibleFluid',steady_iteration_indices_are_not_physical_time=True,
        mesh_cells=cells,actual_settings=declaration,mesh_quality=quality,actual_native_patches=patches,
        native_fields=stats,boundary_mass_balance=mass,final_iteration=1200,
        residuals=dict(final_maximum_printed_initial_residual_per_variable=terminal,
            last100_maximum_printed_initial_residual_per_variable={name:max(row[name] for row in residual_rows[-100:]) for name in terminal},
            definition='maximum printed initial linear-system residual per variable in one steady SIMPLE iteration',
            norm_equivalent_to_fluent=False),
        final100_force_stability=stability,final_forces=force_rows[-1] if force_rows else None,force_samples=len(force_rows),
        check_results=checks,failed_checks=[name for name,value in checks.items() if not value],
        all_declared_checks_passed=all(checks.values()),physical_validation=False,mesh_independence_established=False,
        steady_physical_wake_established=False,benchmark_validation=False,
        interpretation='Bounded numerical checks on one illustrative steady 3D mesh; missed targets retained without automatic extension')
    out.mkdir(parents=True)
    json_new(out/'checks.json',record)
    json_new(out/'source-manifest.json',source)
    for name,rows in [('residual-history.csv',residual_rows),('force-history.csv',force_rows)]:
        with (out/name).open('x',encoding='ascii',newline='') as stream:
            columns=list(rows[0]) if rows else ['iteration','Cm','Cd','Cl']
            writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    for path in (Path(__file__),Path(__file__).with_name('prepare_openfoam_finite_wing_v003.py')):
        copied(path,out/('executed-'+path.name))
    print(json.dumps(dict(output=str(out),all_declared_checks_passed=record['all_declared_checks_passed'],failed_checks=record['failed_checks']),indent=2))
    return 0 if record['all_declared_checks_passed'] else 2


if __name__=='__main__':
    sys.exit(main())
