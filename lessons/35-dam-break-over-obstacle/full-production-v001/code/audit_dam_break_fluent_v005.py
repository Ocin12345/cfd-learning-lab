"""Independently audit native Fluent dam-break states and open-boundary mass."""
from pathlib import Path
import argparse,csv,hashlib,json,re
import h5py,numpy as np
from ansys.fluent.core.offline import _lispy
from analyse_nozzle_fluent_v005 import common_mapping
from dam_break_common_v002 import *
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def water_group(f):
    groups=[phase['cells/SV_VOF'] for phase in f['results/1'].values() if 'cells/SV_VOF' in phase]
    if len(groups)==1:return groups[0]
    groups=[g for g in groups if '/phase-3/' in g.name];assert len(groups)==1;return groups[0]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);ap.add_argument('out',type=Path);ap.add_argument('--mode',choices=['flat','collapse','production'],required=True);args=ap.parse_args();run=args.run;assert not args.out.exists()
    steps=10000 if args.mode=='production' else 400;interval=100 if args.mode=='production' else 20
    a=json.loads((run/'setup.json').read_text());assert a['mode']==args.mode and a['cells']==NCELLS and a['dt_s']==DT and a['steps']==steps and a['save_interval']==interval
    assert json.loads((run/'status.json').read_text())['state']=='closed'
    controls=a['solution_controls']['p_v_controls'];assert controls['flow_courant_number']==20 and controls['explicit_pressure_under_relaxation']==.5 and controls['explicit_momentum_under_relaxation']==.5
    for name,(rho,mu) in MATERIALS.items():assert a['materials'][name]['density']['value']==rho and a['materials'][name]['viscosity']['value']==mu
    mp=a['model'];assert mp['model']=='vof' and mp['vof_parameters']['vof_formulation']=='explicit' and mp['advanced_formulation']['implicit_body_force'] is True
    force=mp['phase_interaction']['forces'];assert force['surface_tension_model'] is True and force['surface_tension_model_type']=='Continuum Surface Force' and force['wall_adhesion'] is False
    assert len(force['surface_tension'])==1 and force['surface_tension'][0]['ppt_methods']=='constant' and force['surface_tension'][0]['ppt_constant']==SIGMA
    op=a['operating_conditions'];assert op['gravity']['components']==[0,-G] and op['operating_density']==dict(method='user-input',value=0.)
    out=a['outlet']['phase'];momentum=out['mixture']['momentum'];assert momentum['gauge_pressure']['value']==0 and momentum['backflow_pressure_spec']=='Static Pressure' and momentum['backflow_dir_spec_method']=='Normal to Boundary'
    assert out['phase-2']['multiphase']['backflow_volume_fraction']['value']==0
    wall=a['wall']['phase']['mixture']['momentum'];assert wall['wall_motion']=='Stationary Wall' and wall['shear_condition']=='No Slip'
    assert a['patch_options']['patch_reconstructed_interface'] is False and a['patch_options']['use_volumetric_smoothing'] is False
    methods=a['numerical_methods'];assert methods['p_v_coupling']['flow_scheme']=='Coupled' and methods['transient_formulation']=='unsteady-1st-order'
    scheme=methods['spatial_discretization']['discretization_scheme'];assert scheme['mom']=='second-order-upwind' and scheme['pressure']=='body-force-weighted' and scheme['mp']=='geo-reconstruct'
    targets={'continuity':1e-12,'x-velocity':1e-10,'y-velocity':1e-10} if args.mode=='flat' else {k:1e-4 for k in ['continuity','x-velocity','y-velocity']}
    reporting=a['residual_reporting'];assert reporting['targets']==targets and reporting['options']['normalize'] is False and reporting['options']['residual_values']['scale_residuals']==(args.mode!='flat')
    for key,value in targets.items():assert reporting['equations'][key]['absolute_criteria']==value and reporting['equations'][key]['check_convergence'] is True
    definition=a['report_definitions']['flux']['mdot_atmosphere'];assert definition['report_type']=='flux-massflow' and definition['phase']=='mixture' and definition['boundaries']==['atmosphere'] and definition['average_over']==1
    monitor=a['report_files']['atmosphere_mass'];assert monitor['frequency_of']=='time-step' and monitor['frequency']==1 and monitor['write_instantaneous_values'] is True and monitor['active'] is True
    text=(run/'solver-transcript.txt').read_text(errors='replace');quality=float(re.search(r'Minimum Orthogonal Quality\s*=\s*('+NUMBER+')',text)[1]);assert quality>.2
    assert float(re.search(r'minimum volume \(m3\):\s*('+NUMBER+')',text)[1])>0
    boundaries=list(re.finditer(r'Flow time = ('+NUMBER+r')s, time step = (\d+)',text));assert len(boundaries)==steps
    pattern=re.compile(r'^\s*(\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+\d+:\d+:\d+\s+\d+\s*$',re.M);history=[];start=0;previous=0
    for k,b in enumerate(boundaries,1):
        assert int(b[2])==k and abs(float(b[1])-k*DT)<1e-10
        rows=list(pattern.finditer(text,start,b.start()));assert rows;last=rows[-1];iteration=int(last[1]);values=[float(last[i]) for i in [2,3,4]];assert np.isfinite(values).all() and 0<iteration-previous<=400
        history.append(dict(step=k,time_s=k*DT,residuals=dict(zip(targets,values)),inner_iterations=iteration-previous,passed=all(v<=limit for v,limit in zip(values,targets.values()))));previous=iteration;start=b.end()
    co=[float(v) for v in re.findall(r'Global Courant Number \[Explicit VOF Criteria\]\s*:\s*('+NUMBER+')',text)];assert len(co)==steps and np.isfinite(co).all()
    fluxpath=run/'massflow-history.out';ft=fluxpath.read_text();assert '"Time Step"' in ft and '"mdot_atmosphere"' in ft
    rows=[re.fullmatch(r'\s*(\d+)\s+('+NUMBER+r')\s*',line) for line in ft.splitlines()];rows=[m for m in rows if m]
    assert len(rows)==steps and [int(m[1]) for m in rows]==list(range(1,steps+1));qout=-np.array([float(m[2]) for m in rows]);assert np.isfinite(qout).all()
    integrated=np.r_[0,np.cumsum(qout)*DT]
    g=geometry();mapping=common_mapping(run/'dam.cas.h5',g);expected=initial_fields(g,args.mode);M0=expected['mixture_mass_per_m']
    with (run/'frames.csv').open(newline='') as stream:ledger=list(csv.DictReader(stream))
    assert [int(r['step']) for r in ledger]==list(range(0,steps+1,interval))
    records=[];sources=[run/name for name in ['dam.cas.h5','dam-final.cas.h5','dam-final.dat.h5','geometry.npz','setup.json','frames.csv','solver-transcript.txt','massflow-history.out','initial-field-check.json','executed-driver.py']]
    for row in ledger:
        step=int(row['step']);path=run/row['file'];fields=[]
        with h5py.File(path) as f:
            clock={v[0]:v[1] for v in _lispy.parse(f['settings/Data Variables'][0].decode())[1]};assert int(clock['time-step'])==step and abs(float(clock['flow-time'])-step*DT)<1e-10 and abs(float(clock['physical-time-step'])-DT)<1e-14
            assert abs(float(row['time_s'])-float(clock['flow-time']))<1e-10
            for group in [water_group(f),f['results/1/phase-1/cells/SV_U'],f['results/1/phase-1/cells/SV_V'],f['results/1/phase-1/cells/SV_P']]:
                raw=cell_values(group);aligned=np.empty(NCELLS);aligned[mapping]=raw;fields.append(aligned)
        alpha,u,v,pressure=fields;area=float(alpha@g['areas']);mass=float((1+999*alpha)@g['areas']);speed=float(np.hypot(u,v).max());change=float(abs(alpha-expected['alpha']).max())
        if step==0:assert change<1e-10 and np.max(abs(pressure-expected['pressure']))<1e-8 and abs(area-WATER_AREA)<1e-12 and speed==0
        assert abs(float(row['water_area_m2'])-area)<1e-12 and abs(float(row['mixture_mass_per_m'])-mass)<1e-10
        closure=(mass-M0+integrated[step])/M0;bounds=alpha.min()>=-1e-6 and alpha.max()<=1+1e-6;quiet=args.mode!='flat' or (speed<=1e-4 and change<=.01)
        records.append(dict(step=step,time_s=step*DT,water_area_m2=area,mixture_mass_per_m=mass,integrated_outward_mass_per_m=float(integrated[step]),relative_mixture_mass_closure=float(closure),max_speed_m_s=speed,maximum_alpha_change=change,alpha_min=float(alpha.min()),alpha_max=float(alpha.max()),pressure_min_Pa=float(pressure.min()),pressure_max_Pa=float(pressure.max()),fields_pass=bool(bounds and quiet and abs(closure)<=.005)));sources.append(path)
    proof={}
    for name in ['solver-inputs-v004.json','fluent-inputs-v007.json']:
        path=P/name;receipt=json.loads(path.read_text())
        for source,h in receipt['source_sha256'].items():assert sha(Path(source))==h;proof[source]=h
        sources.append(path)
    sources.append(Path(__file__))
    for source in sources:proof[str(source)]=sha(source)
    passed=all(r['passed'] for r in history) and all(r['fields_pass'] for r in records) and max(co)<=.25
    result=dict(solver='ANSYS Fluent2026R1',mode=args.mode,cells=NCELLS,steps=steps,dt_s=DT,saved_states=len(records),final_time_s=steps*DT,numerical_checks_pass=bool(passed),reported_interface_courant_count=len(co),expected_courant_reports=steps,peak_interface_courant=max(co),native_flux_rows=len(rows),maximum_relative_mixture_mass_closure=max(abs(r['relative_mixture_mass_closure']) for r in records),maximum_saved_speed_m_s=max(r['max_speed_m_s'] for r in records),minimum_orthogonal_quality=quality,mass_flux_sign='Fluent native inward-positive; outward=-reported, planar unit depth',mass_integration='Euler: end-of-step native mixture mass flux, every actual step',records=records,step_records=history,native_source_sha256=proof,physical_validation=False,production_complete=args.mode=='production' and bool(passed))
    args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['records','step_records','native_source_sha256']},indent=2))
    if not passed:raise SystemExit(2)
if __name__=='__main__':main()