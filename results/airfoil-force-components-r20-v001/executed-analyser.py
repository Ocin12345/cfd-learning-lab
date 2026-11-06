"""Aggregate hash-bound native Fluent/OpenFOAM pressure and viscous forces.

Standard library only. Native evidence is read-only; all outputs must be fresh.
Failed main convergence/symmetry checks remain failed in the resulting record.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil

ROOT = Path('E:/CFD-Learning-Lab' if os.name == 'nt' else '/mnt/e/CFD-Learning-Lab')
NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
COEFFICIENT_EQUIVALENCE_TOLERANCE = 1e-8


def require(value, message):
    if not value:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def local_path(value):
    if os.name != 'nt' and re.match(r'^E:[/\\]', str(value), re.I):
        return Path('/mnt/e/' + str(value)[3:].replace('\\', '/'))
    return Path(value)


def key(path):
    return str(path.resolve()).casefold()


def parse_native_force(raw):
    rows = re.findall(r'^\s*Net\s+((?:'+NUMBER+r'\s+){5}'+NUMBER+r')\s*$', raw, re.M)
    require(len(rows) == 1, 'Expected one native scalar pressure/viscous/total force row')
    values = list(map(float, rows[0].split()))
    require(all(math.isfinite(v) for v in values), 'Nonfinite native force')
    return dict(zip(['pressure_force_N','viscous_force_N','total_force_N',
                     'pressure_coefficient','viscous_coefficient','total_coefficient'], values))


def coefficient_row(raw, time):
    rows = [list(map(float, line.split())) for line in raw.splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    selected = [row for row in rows if row[0] == time]
    require(len(selected) == 1 and len(selected[0]) == 6
            and all(math.isfinite(v) for v in selected[0]), 'Missing/invalid native coefficient row')
    return dict(zip(['Cm','Cd','Cl'], selected[0][1:4]))


def force_vectors(raw, time):
    rows = [line for line in raw.splitlines() if line.strip() and not line.lstrip().startswith('#')]
    require(len(rows) == 1, 'Expected one saved-time force decomposition row')
    values = list(map(float, re.findall(NUMBER, rows[0])))
    require(len(values) == 13 and values[0] == time and all(math.isfinite(v) for v in values),
            'Expected pressure/viscous force and moment vectors without porosity')
    return values[1:4], values[4:7], values[7:10], values[10:13]


def check_diagnostic_flags(record):
    checks = record['check_results']
    require(isinstance(checks,dict) and checks and all(type(v) is bool for v in checks.values()),
            'Main diagnostic checks must be named booleans')
    failed = [name for name,value in checks.items() if not value]
    require(record['failed_checks'] == failed and type(record['numerical_checks_passed']) is bool
            and record['numerical_checks_passed'] == (not failed), 'Main failed flags are inconsistent')
    convergence = [value for name,value in checks.items() if not name.endswith('zero_angle_lift_diagnostic')]
    symmetry = [value for name,value in checks.items() if name.endswith('zero_angle_lift_diagnostic')]
    require(record['numerical_convergence_passed'] == all(convergence)
            and record['symmetry_diagnostic_passed'] == all(symmetry),
            'Main convergence/symmetry classification changed')


def signed_share(component, total):
    return None if total == 0 else 100*component/total


def report_close(actual, reported):
    # Native Fluent reports print fewer digits than saved fields; check rounding.
    return math.isclose(actual, reported, rel_tol=1e-7, abs_tol=1e-12)


def analyse(fluent_run, component_path, comparison_dir):
    inputs = {}
    def bind(path, expected=None):
        path = local_path(path).resolve()
        require(path.is_relative_to(ROOT.resolve()) and path.is_file(), 'Missing/outside-project evidence '+str(path))
        digest = sha256(path)
        require(expected is None or digest == expected, 'Evidence SHA256 mismatch: '+str(path))
        previous = inputs.get(key(path))
        require(previous is None or previous['sha256'] == digest, 'Evidence changed during analysis')
        inputs[key(path)] = dict(path=str(path),sha256=digest)
        return path
    def read(path):
        return json.loads(bind(path).read_text())
    main = read(comparison_dir/'comparison.json')
    check_diagnostic_flags(main)
    main_manifest = read(comparison_dir/'source-manifest.json')
    require(main_manifest['numerical_checks_passed'] == main['numerical_checks_passed']
            and main_manifest['failed_checks'] == main['failed_checks'], 'Main record/manifest diagnostic state differs')
    main_inputs = {}
    for item in main_manifest['inputs']:
        path = bind(item['path'],item['sha256'])
        main_inputs[key(path)] = item['sha256']
    def main_bound(path):
        path = local_path(path).resolve()
        require(key(path) in main_inputs, 'Selected native evidence is not bound to the main analysis: '+str(path))
        return bind(path,main_inputs[key(path)])
    of = read(component_path)
    scratch = component_path.parent.resolve()
    preparation_path = bind(scratch/'input-manifest.json',of['input_manifest_sha256'])
    preparation = json.loads(preparation_path.read_text())
    case = local_path(preparation['source_case']).resolve()
    time = preparation['selected_iteration']
    require(of['selected_iteration'] == time == main['openfoam']['final_iteration'], 'Saved OpenFOAM times differ')
    status = read(case/'run-status.json')
    require(status['stage'] == 'finished_execution' and status['final_iteration'] == time, 'OpenFOAM source not final/reconstructed')
    require(of['total_equivalence_passed'] is True and of['flow_iterations_executed'] == 0
            and of['physical_validation'] is False
            and of['absolute_coefficient_equivalence_tolerance'] == COEFFICIENT_EQUIVALENCE_TOLERANCE
            and preparation['absolute_coefficient_equivalence_tolerance'] == COEFFICIENT_EQUIVALENCE_TOLERANCE,
            'Postprocessing equivalence/zero-flow contract differs')
    for item in preparation['sources']:
        path = main_bound(item['source'])
        require(sha256(path) == item['sha256'], 'Scratch source differs from main evidence')
        copied = (scratch/item['copy']).resolve()
        require(copied.is_relative_to(scratch), 'Copied evidence escapes the selected scratch')
        bind(copied,item['sha256'])
    for item in preparation['generated']:
        generated = (scratch/item['copy']).resolve()
        require(generated.is_relative_to(scratch), 'Generated input escapes the selected scratch')
        bind(generated,item['sha256'])
    for item in of['outputs']:
        path = (scratch/item['path']).resolve()
        require(path.is_relative_to(scratch), 'Postprocessing output escapes scratch')
        bind(path,item['sha256'])
    raw_forces = list((scratch/'postProcessing/forcesDiagnostic').glob('*/forces.dat'))
    raw_replay = list((scratch/'postProcessing/coefficientsReplay').glob('*/forceCoeffs.dat'))
    require(len(raw_forces) == len(raw_replay) == 1, 'Expected one force/replay output')
    expected_outputs = {key(raw_forces[0]),key(raw_replay[0]),key(scratch/'postprocess.log')}
    require({key(scratch/item['path']) for item in of['outputs']} == expected_outputs,
            'Result manifest does not bind the actual force/replay/log files')
    log = (scratch/'postprocess.log').read_text()
    require(re.findall(r'^Time =\s*(\S+)',log,re.M) in [[str(time)],[str(time)+'s']]
            and re.search(r'^nProcs\s*:\s*1\s*$',log,re.M)
            and re.search(r'^End\s*$',log,re.M) and 'FOAM FATAL' not in log
            and not re.search(r'Solving for|SIMPLE solution converged|PIMPLE: iteration',log),
            'Unexpected postprocessing time/process/equation-solving evidence')
    controls = read(case/'declared-comparison.json')
    require(sha256(bind(comparison_dir/'declared-comparison.json')) == sha256(case/'declared-comparison.json'),
            'Public/native comparison declarations differ')
    refs = preparation['force_references']
    for name,target in [('rhoInf','density_kg_m3'),('magUInf','velocity_m_s'),('Aref','reference_area_m2'),('lRef','chord_m')]:
        require(refs[name] == controls[target] and refs[name] > 0, 'Actual force references differ')
    require(preparation['pRef_Pa'] == controls['pressure_reference_Pa_gauge'], 'Force pressure reference differs')
    denominator = .5*refs['rhoInf']*refs['magUInf']**2*refs['Aref']
    pf,vf,pm,vm = force_vectors(raw_forces[0].read_text(),time)
    openfoam = {}
    openfoam_moments = {}
    for name,force,moment in [('pressure',pf,pm),('viscous',vf,vm),
                              ('total',[a+b for a,b in zip(pf,vf)],[a+b for a,b in zip(pm,vm)])]:
        openfoam[name] = dict(Cd=force[0]/denominator,Cl=force[1]/denominator,
                              drag_force_N=force[0],lift_force_N=force[1])
        openfoam_moments[name] = moment[2]/(denominator*refs['lRef'])
        saved = next((row for row in of['components'] if row['component'] == name),None)
        require(saved is not None and all(saved[k] == openfoam[name][k] for k in ['Cd','Cl'])
                and saved['Fx_N'] == force[0] and saved['Fy_N'] == force[1]
                and saved['Cm'] == openfoam_moments[name], 'Raw/saved OpenFOAM decomposition differs')
    expected = coefficient_row((scratch/'evidence/source-forceCoeffs.dat').read_text(),time)
    replay = coefficient_row(raw_replay[0].read_text(),time)
    require(expected == preparation['expected_coefficients'] == of['expected_saved_forceCoeffs']
            and replay == of['replayed_forceCoeffs'], 'Native/source/replayed coefficient records differ')
    for coefficient in ['Cd','Cl']:
        require(abs(openfoam['total'][coefficient]-expected[coefficient]) <= COEFFICIENT_EQUIVALENCE_TOLERANCE
                and abs(replay[coefficient]-expected[coefficient]) <= COEFFICIENT_EQUIVALENCE_TOLERANCE
                and abs(main['openfoam'][coefficient]-expected[coefficient]) <= 1e-12,
                'OpenFOAM totals disagree with saved main result')
    require(abs(openfoam_moments['total']-expected['Cm']) <= COEFFICIENT_EQUIVALENCE_TOLERANCE
            and abs(replay['Cm']-expected['Cm']) <= COEFFICIENT_EQUIVALENCE_TOLERANCE,
            'Postprocessed moment normalization/reference differs')
    geometry = read(main_bound(fluent_run/'geometry-setup.json'))
    require(json.loads(main_bound(case/'source-setup.json').read_text()) == geometry
            and geometry['cells'] == main['cells'] and geometry['angle_degrees'] == main['angle_degrees']
            and sha256(main_bound(case/'input-airfoil.msh')) == main['source_mesh_sha256'], 'Matched mesh/setup differs')
    fluent_status = read(fluent_run/'status.json')
    require(fluent_status['stage'] == 'completed', 'Fluent source is not completed')
    settings = read(main_bound(fluent_run/'loaded-settings.json'))
    require(settings['general']['solver']['type'] == 'pressure-based'
            and settings['general']['solver']['time'] == 'steady'
            and settings['general']['solver']['two_dim_space'] == 'planar'
            and not settings['general']['operating_conditions']['gravity']['enable']
            and settings['viscous']['model'] == 'laminar'
            and settings['methods'] == main['fluent']['numerical_methods'], 'Fluent model/methods differ')
    for name,value in [('density',refs['rhoInf']),('velocity',refs['magUInf']),('area',refs['Aref']),
                       ('length',refs['lRef']),('depth',controls['span_m']),('pressure',preparation['pRef_Pa'])]:
        require(settings['references'][name] == value, 'Fluent coefficient reference differs: '+name)
    require(settings['material']['density']['value'] == refs['rhoInf']
            and settings['material']['viscosity']['value'] == controls['dynamic_viscosity_Pa_s']
            and settings['material']['density']['option'] == settings['material']['viscosity']['option'] == 'value'
            and settings['inlet']['momentum']['velocity_specification_method'] == 'Components'
            and all(v['option'] == 'value' for v in settings['inlet']['momentum']['velocity_components'])
            and [v['value'] for v in settings['inlet']['momentum']['velocity_components']] == [refs['magUInf'],0]
            and settings['outlet']['momentum']['gauge_pressure']['value'] == preparation['pRef_Pa'],
            'Fluent physical conditions differ')
    for name in ['airfoil-Re100-final.cas.h5','airfoil-Re100-final.dat.h5']:
        main_bound(fluent_run/name)
    require(main['fluent']['saved_data_cell_ranges_and_material_verified'] is True
            and main['fluent']['native_mesh_connectivity_verified'] is True,
            'Main native HDF field/connectivity verification missing')
    criteria = read(main_bound(fluent_run/'declared-targets.json'))
    summary = read(main_bound(fluent_run/'summary.json'))
    require(criteria == main['fluent']['acceptance_thresholds'], 'Fluent acceptance criteria differ')
    with main_bound(fluent_run/'force-history.csv').open(newline='') as stream:
        history = [{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]
    require([row['additional_iteration'] for row in history] == list(range(criteria['additional_iterations']+1))
            and all(math.isfinite(v) for row in history for v in row.values()), 'Incomplete/nonfinite Fluent history')
    terminal = {}
    for row in history:
        for direction in ['drag','lift']:
            parsed = parse_native_force(main_bound(fluent_run/f"{direction}-{int(row['additional_iteration']):03d}.txt").read_text())
            require(parsed['total_force_N'] == row[direction+'_force_N']
                    and parsed['total_coefficient'] == row[direction+'_coefficient'], 'History/native scalar force differs')
            if direction == 'drag':
                require(parsed['pressure_force_N'] == row['pressure_drag_N']
                        and parsed['viscous_force_N'] == row['viscous_drag_N'], 'History/native force components differ')
            for component in ['pressure','viscous','total']:
                require(report_close(parsed[component+'_force_N']/denominator,parsed[component+'_coefficient']),
                        'Native force/coefficient normalization differs')
            require(report_close(parsed['pressure_force_N']+parsed['viscous_force_N'],parsed['total_force_N'])
                    and report_close(parsed['pressure_coefficient']+parsed['viscous_coefficient'],parsed['total_coefficient']),
                    'Native component sums exceed report rounding')
            terminal[direction] = parsed
    require(history[-1] == summary['final'], 'Stale Fluent summary')
    fluent = {name:dict(Cd=terminal['drag'][name+'_coefficient'],Cl=terminal['lift'][name+'_coefficient'],
                        drag_force_N=terminal['drag'][name+'_force_N'],lift_force_N=terminal['lift'][name+'_force_N'])
              for name in ['pressure','viscous','total']}
    require(all(fluent['total'][name] == main['fluent'][name] for name in ['Cd','Cl']), 'Native Fluent/main endpoint differs')
    window = history[-criteria['last_window_samples']:]
    require(len(window) == criteria['last_window_samples'] and window, 'Missing final Fluent force window')
    cd,cl = [[row[name+'_coefficient'] for row in window] for name in ['drag','lift']]
    cd_range = (max(cd)-min(cd))/abs(math.fsum(cd)/len(cd));cl_range=max(cl)-min(cl)
    text = main_bound(fluent_run/'solver-transcript.txt').read_text()
    residuals = re.findall(r'^\s*(\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')',text,re.M)
    require(residuals and all(math.isfinite(float(v)) for v in residuals[-1]), 'Invalid native residual evidence')
    final_residuals = list(map(float,residuals[-1][1:]))
    require(all(v >= 0 for v in final_residuals), 'Negative native residual evidence')
    quality = float(re.search(r'Minimum Orthogonal Quality\s*=\s*('+NUMBER+r')',text)[1])
    require(quality >= criteria['minimum_orthogonal_quality']
            and quality == main['fluent']['minimum_orthogonal_quality'], 'Native mesh admission differs')
    mass = main_bound(fluent_run/'mass-balance.txt').read_text()
    net=float(re.search(r'^\s*Net\s+('+NUMBER+r')',mass,re.M)[1]);inlet=float(re.search(r'^\s*inlet\s+('+NUMBER+r')',mass,re.M)[1])
    require(inlet != 0, 'Invalid native inlet mass flux')
    recomputed = {'fluent.residuals':max(final_residuals)<=criteria['residual_target'],
                  'fluent.Cd_stability':cd_range<=criteria['drag_relative_range_max'],
                  'fluent.Cl_stability':cl_range<=criteria['lift_absolute_range_max'],
                  'fluent.zero_angle_lift_diagnostic':max(map(abs,cl))<=criteria['zero_angle_abs_lift_diagnostic_target'],
                  'fluent.mass_balance':abs(net/inlet)<=criteria['relative_mass_imbalance_max']}
    require(all(main['check_results'][name] == value for name,value in recomputed.items())
            and int(residuals[-1][0]) == main['fluent']['final_iteration'], 'Main/native Fluent diagnostics differ')
    differences = []
    total_gap = {name:fluent['total'][name]-openfoam['total'][name] for name in ['Cd','Cl']}
    cross_percent = 100*abs(total_gap['Cd'])/abs(fluent['total']['Cd'])
    require(abs(cross_percent-main['cross_code_Cd_difference_percent']) <= 1e-10
            and main['cross_code_learning_target_percent'] == controls['cross_code_cd_difference_learning_target_percent']
            and main['cross_code_target_pass'] == (cross_percent <= main['cross_code_learning_target_percent']),
            'Main cross-code difference/criterion flags differ from native endpoints')
    for component in ['pressure','viscous','total']:
        delta = {name:fluent[component][name]-openfoam[component][name] for name in ['Cd','Cl']}
        differences.append(dict(component=component,Cd_Fluent_minus_OpenFOAM=delta['Cd'],Cl_Fluent_minus_OpenFOAM=delta['Cl'],
                                Cd_signed_share_of_total_gap_percent=signed_share(delta['Cd'],total_gap['Cd']),
                                Cl_signed_share_of_total_gap_percent=signed_share(delta['Cl'],total_gap['Cl'])))
    flag_names = ['check_results','failed_checks','numerical_checks_passed','numerical_convergence_passed',
                  'symmetry_diagnostic_passed','cross_code_target_pass']
    record = dict(kind='Saved force-component accounting; no new CFD',
                  parameters=dict(Re=main['reynolds_number'],angle_degrees=main['angle_degrees'],cells=main['cells'],
                                  outer_radius_m=geometry['outer_radius_m'],span_m=main['span_m']),
                  fluent=fluent,openfoam=openfoam,differences=differences,
                  **{name:main[name] for name in flag_names},
                  evidence_consistency_checks_passed=True,
                  main_diagnostic_values=dict(openfoam=main['openfoam'],fluent={key:main['fluent'][key] for key in
                      ['final_iteration','final_scaled_residuals','Cd_relative_range','Cl_absolute_range',
                       'maximum_absolute_Cl_in_window','relative_mass_imbalance']}),
                  main_acceptance_thresholds=main['acceptance_thresholds'],
                  native_Fluent_checks_recomputed=recomputed,
                  selected_openfoam_iteration=time,selected_fluent_iteration=int(residuals[-1][0]),
                  cross_code_Cd_difference_percent=main['cross_code_Cd_difference_percent'],
                  cross_code_learning_target_percent=main['cross_code_learning_target_percent'],
                  native_HDF_evidence_hash_bound_to_main_checks=True,flow_iterations_executed=0,
                  numerical_cause_identified=False,physical_validation=False,
                  component_closure_residual={name:sum(fluent[c][name] for c in ['pressure','viscous'])-fluent['total'][name]
                                              for name in ['Cd','Cl']},
                  note='Signed component shares locate the observed coefficient gap, not its numerical cause. Shares can be negative or exceed 100% when terms cancel; an exactly zero total gap has no defined share. Native Fluent report rounding limits closure and percentage precision. Main failed convergence/symmetry and cross-code flags are preserved.')
    for item in inputs.values():
        require(sha256(Path(item['path'])) == item['sha256'], 'Input changed before output')
    return record,list(inputs.values())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fluent-run',type=Path,required=True)
    parser.add_argument('--openfoam-components',type=Path,required=True,help='Verified scratch force-components.json')
    parser.add_argument('--comparison-dir',type=Path,required=True,help='Main matched-pair comparison/manifest directory')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();fluent,components,comparison,out=[p.resolve() for p in
        [args.fluent_run,args.openfoam_components,args.comparison_dir,args.out]]
    allowed_output = out.is_relative_to((ROOT/'repo/cfd-learning-lab/results').resolve()) or out.is_relative_to((ROOT/'local').resolve())
    require(not out.exists() and allowed_output
            and not out.is_relative_to(fluent) and not out.is_relative_to(components.parent), 'Use a fresh project analysis directory')
    record,inputs=analyse(fluent,components,comparison)
    out.mkdir(parents=True)
    with (out/'comparison.json').open('x',encoding='utf-8') as stream:json.dump(record,stream,indent=2,allow_nan=False);stream.write('\n')
    rows=[dict(solver=solver,component=component,**values) for solver in ['fluent','openfoam'] for component,values in record[solver].items()]
    for name,data in [('force-components.csv',rows),('force-differences.csv',record['differences'])]:
        with (out/name).open('x',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    manifest=dict(script_sha256=sha256(Path(__file__)),evidence_consistency_checks_passed=True,
                  numerical_checks_passed=record['numerical_checks_passed'],failed_checks=record['failed_checks'],inputs=inputs)
    with (out/'source-manifest.json').open('x',encoding='utf-8') as stream:json.dump(manifest,stream,indent=2);stream.write('\n')
    shutil.copy2(Path(__file__),out/'executed-analyser.py')
    print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__':main()
