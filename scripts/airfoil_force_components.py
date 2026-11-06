"""Prepare or verify an isolated Foundation-14 saved-field force diagnostic.

Preparation never runs OpenFOAM. The generated shell runs foamPostProcess only.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shlex
import shutil

ROOT = Path('E:/CFD-Learning-Lab' if os.name == 'nt' else '/mnt/e/CFD-Learning-Lab')
NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
ABS_COEFFICIENT_TOLERANCE = 1e-8


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
    if os.name != 'nt' and re.match(r'^E:[/\\]', value, re.I):
        return Path('/mnt/e/' + value[3:].replace('\\', '/'))
    return Path(value)


def linux_path(path):
    value = str(path).replace('\\', '/')
    if value.lower().startswith('e:/'):
        return '/mnt/e/' + value[3:]
    require(value.startswith('/mnt/e/'), 'Only project E paths are supported')
    return value


def clean(path):
    raw = re.sub(r'/\*.*?\*/|//[^\n]*', '', path.read_text(), flags=re.S)
    require(not re.search(r'^\s*#|\bcoded\b|\bcodeStream\b', raw, re.M),
            'External include/code dependency is outside this bounded diagnostic: ' + str(path))
    return raw


def block(raw, name):
    match = re.search(r'\b' + re.escape(name) + r'\s*\{', raw)
    require(match, 'Missing dictionary block ' + name)
    first, depth = match.end(), 1
    for i in range(first, len(raw)):
        depth += (raw[i] == '{') - (raw[i] == '}')
        if depth == 0:
            return raw[first:i]
    raise ValueError('Unclosed block ' + name)


def entry(raw, key):
    match = re.search(r'\b' + re.escape(key) + r'\s+([^;]+);', raw)
    require(match, 'Missing entry ' + key)
    return match[1].strip()


def parse_coefficients(raw, time):
    rows = [list(map(float, line.split())) for line in raw.splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    selected = [row for row in rows if row[0] == time]
    require(len(selected) == 1 and len(selected[0]) == 6
            and all(math.isfinite(v) for v in selected[0]), 'Expected one finite coefficient row at selected time')
    return dict(zip(['Cm', 'Cd', 'Cl'], selected[0][1:4]))


def parse_forces(raw, time):
    rows = [line for line in raw.splitlines() if line.strip() and not line.lstrip().startswith('#')]
    require(len(rows) == 1, 'Expected one saved-time force decomposition row')
    numbers = list(map(float, re.findall(NUMBER, rows[0])))
    require(len(numbers) == 13 and numbers[0] == time and all(math.isfinite(v) for v in numbers),
            'Expected pressure/viscous force and moment vectors only (no porous contribution)')
    return dict(pressure_force=numbers[1:4], viscous_force=numbers[4:7],
                pressure_moment=numbers[7:10], viscous_moment=numbers[10:13])


def dump(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def prepare(case, time, out):
    case, out = case.resolve(), out.resolve()
    require(case.is_relative_to(ROOT.resolve()) and out.is_relative_to(ROOT.resolve())
            and not out.is_relative_to(case) and not out.exists(), 'Use a fresh external project scratch directory')
    require(time > 0, 'Select an explicit positive saved iteration')
    status = json.loads((case/'run-status.json').read_text())
    require(status['stage'] == 'finished_execution' and status['final_iteration'] == time,
            'Selected case must have completed execution/reconstruction at the selected iteration')
    targets = json.loads((case/'declared-comparison.json').read_text())
    control = clean(case/'system/controlDict')
    coefficients = block(block(control, 'functions'), 'forceCoefficients')
    for key, value in [('type','forceCoeffs'), ('rho','rhoInf'), ('patches','(airfoil)'),
                       ('p','p'), ('U','U'), ('dragDir','(1 0 0)'), ('liftDir','(0 1 0)'),
                       ('pitchAxis','(0 0 1)')]:
        require(entry(coefficients,key) == value, 'Unsupported source force reference ' + key)
    refs = {key:float(entry(coefficients,key)) for key in ['rhoInf','magUInf','Aref','lRef']}
    for key, target in [('rhoInf','density_kg_m3'), ('magUInf','velocity_m_s'),
                        ('Aref','reference_area_m2'), ('lRef','chord_m')]:
        require(refs[key] == targets[target] and refs[key] > 0, 'Declared/actual force reference differs: ' + key)
    p_ref = float(entry(coefficients,'pRef')) if re.search(r'\bpRef\s',coefficients) else 0.
    require(p_ref == targets['pressure_reference_Pa_gauge'], 'Pressure reference differs')
    require(entry(clean(case/'constant/momentumTransport'),'simulationType') == 'laminar'
            and float(entry(clean(case/'constant/physicalProperties'),'nu')) == targets['kinematic_viscosity_m2_s'],
            'Unsupported momentum/viscosity model')
    expected_candidates = []
    for path in sorted((case/'postProcessing/forceCoefficients').glob('*/forceCoeffs.dat')):
        raw = path.read_text()
        if any(line.strip() and not line.lstrip().startswith('#') and float(line.split()[0]) == time
               for line in raw.splitlines()):
            expected_candidates.append((path,parse_coefficients(raw,time)))
    require(expected_candidates, 'No source forceCoeffs row at selected saved time')
    expected_path, expected = expected_candidates[-1]
    require(all(all(abs(value[k]-expected[k]) <= 1e-12 for k in expected)
                for _,value in expected_candidates), 'Inconsistent source checkpoint force coefficients')
    pairs = [(case/str(time)/name,Path(str(time))/name) for name in ['U','p','phi']]
    mesh_root = case/'0/polyMesh'
    require(mesh_root.is_dir() and not (case/str(time)/'polyMesh').exists(), 'Expected fixed native mesh in 0/polyMesh')
    pairs += [(mesh_root/name,Path('constant/polyMesh')/name)
              for name in ['points','faces','owner','neighbour','boundary']]
    pairs += [(case/'constant'/name,Path('constant')/name) for name in ['physicalProperties','momentumTransport']]
    pairs += [(case/'system'/name,Path('system')/name) for name in ['fvSchemes','fvSolution']]
    pairs += [(case/name,Path('evidence')/name) for name in ['declared-comparison.json','run-status.json']]
    pairs += [(case/'system/controlDict',Path('evidence/source-controlDict')),
              (expected_path,Path('evidence/source-forceCoeffs.dat'))]
    require(set(p.name for p in (case/'constant').iterdir()) == {'physicalProperties','momentumTransport'},
            'Extra constant model dependencies require explicit review')
    for path,destination in pairs:
        require(path.is_file() and not path.is_symlink(), 'Missing or linked input ' + str(path))
        if destination.parts[0] != 'evidence':
            clean(path)
    out.mkdir(parents=True)
    records = []
    for source,destination in pairs:
        target = out/destination
        target.parent.mkdir(parents=True,exist_ok=True)
        initial_hash = sha256(source)
        shutil.copy2(source,target)
        require(sha256(target) == initial_hash == sha256(source), 'Input changed during copy')
        records.append(dict(source=str(source),copy=destination.as_posix(),sha256=initial_hash))
    header = lambda name: 'FoamFile\n{ version 2.0; format ascii; class dictionary; object '+name+'; }\n\n'
    (out/'system/controlDict').write_text(header('controlDict')+f'''solver incompressibleFluid;
startFrom startTime;
startTime {time};
stopAt endTime;
endTime {time};
deltaT 1;
writeControl timeStep;
writeInterval 1;
writeFormat ascii;
writePrecision 14;
writeCompression off;
runTimeModifiable false;
''',encoding='ascii')
    common = f'''libs ("libforces.so");
patches (airfoil);
rho rhoInf;
rhoInf {refs['rhoInf']:.17g};
p p;
U U;
pRef {p_ref:.17g};
CofR {entry(coefficients,'CofR')};
log true;
writeControl timeStep;
writeInterval 1;
writePrecision 14;
'''
    replay = f'''dragDir (1 0 0); liftDir (0 1 0); pitchAxis (0 0 1);
magUInf {refs['magUInf']:.17g}; Aref {refs['Aref']:.17g}; lRef {refs['lRef']:.17g};
'''
    (out/'system/functions').write_text(header('functions')+'forcesDiagnostic\n{\ntype forces;\n'+common+
                                       '}\ncoefficientsReplay\n{\ntype forceCoeffs;\n'+common+replay+'}\n',encoding='ascii')
    shutil.copy2(Path(__file__),out/'executed-helper.py')
    scratch = shlex.quote(linux_path(out))
    shell = f'''#!/usr/bin/env bash
set -eo pipefail
scratch={scratch}
cd "$scratch"
test ! -e "$scratch/postprocess.log"
export TMPDIR="$scratch/tmp"
export XDG_CACHE_HOME="$scratch/cache"
export WM_PROJECT_USER_DIR=/mnt/e/CFD-Learning-Lab/environments/openfoam-user
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$WM_PROJECT_USER_DIR"
source /opt/openfoam14/etc/bashrc
python3 "$scratch/executed-helper.py" verify --scratch "$scratch" --inputs-only
foamVersion > "$scratch/openfoam-version.txt"
sha256sum /opt/openfoam14/applications/utilities/postProcessing/foamPostProcess/foamPostProcess.C /opt/openfoam14/src/functionObjects/forces/forcesBase/forcesBase.C /opt/openfoam14/applications/modules/incompressibleFluid/incompressibleFluid.C > "$scratch/runtime-source-hashes.txt"
foamPostProcess -case "$scratch" -solver incompressibleFluid -dict "$scratch/system/functions" -time {time} > "$scratch/postprocess.log" 2>&1
python3 "$scratch/executed-helper.py" verify --scratch "$scratch"
'''
    (out/'run-postprocess.sh').write_text(shell,encoding='ascii',newline='\n')
    generated = ['system/controlDict','system/functions','executed-helper.py','run-postprocess.sh']
    manifest = dict(kind='Saved-field pressure/viscous force diagnostic; no flow iterations',source_case=str(case),
                    selected_iteration=time,force_references=refs,pRef_Pa=p_ref,expected_coefficients=expected,
                    absolute_coefficient_equivalence_tolerance=ABS_COEFFICIENT_TOLERANCE,
                    resource_contract=dict(host_available_ram_before_min_bytes=4*1024**3,
                                           project_disk_free_before_min_bytes=10*1024**3,
                                           assumed_peak_extra_ram_bytes=2*1024**3,
                                           host_available_ram_after_min_bytes=2*1024**3,
                                           simultaneous_postprocessors=1),
                    model='Foundation14 incompressibleFluid laminar',sources=records,
                    generated=[dict(copy=name,sha256=sha256(out/name)) for name in generated],
                    official_method_sources=['https://cpp.openfoam.org/v14/foamPostProcess_8C.html',
                    'https://doc.cfd.direct/openfoam/user-guide-v14/post-processing-functionality'],
                    scope='Decomposition localizes a force discrepancy; it does not establish its cause or physical validation')
    dump(out/'input-manifest.json',manifest)
    verify(out,True)
    print('Prepared only; not executed: wsl -d CFD-Foundation14 -- bash '+shlex.quote(linux_path(out/'run-postprocess.sh')))


def verify(scratch, inputs_only=False):
    scratch = scratch.resolve()
    require(scratch.is_relative_to(ROOT.resolve()), 'Scratch must stay under project root')
    manifest = json.loads((scratch/'input-manifest.json').read_text())
    for row in manifest['sources']:
        require(sha256(local_path(row['source'])) == row['sha256'] == sha256(scratch/row['copy']),
                'Original or copied native input changed: '+row['copy'])
    for row in manifest['generated']:
        require(sha256(scratch/row['copy']) == row['sha256'],'Generated diagnostic input changed')
    if inputs_only:
        print('All original/copy/generated hashes match; postprocessing not executed by this check')
        return
    time = manifest['selected_iteration']
    raw_forces = list((scratch/'postProcessing/forcesDiagnostic').glob('*/forces.dat'))
    raw_coefficients = list((scratch/'postProcessing/coefficientsReplay').glob('*/forceCoeffs.dat'))
    require(len(raw_forces) == len(raw_coefficients) == 1,'Expected one postprocessing output per function')
    log = (scratch/'postprocess.log').read_text()
    require('FOAM FATAL' not in log and re.search(r'^End\s*$',log,re.M)
            and not re.search(r'Solving for|SIMPLE solution converged|PIMPLE: iteration',log),
            'Postprocessing failed or unexpected equation-solving evidence')
    selected_times = re.findall(r'^Time =\s*(\S+)',log,re.M)
    require(len(selected_times) == 1 and float(selected_times[0].rstrip('s')) == time,'Unexpected processing times')
    force = parse_forces(raw_forces[0].read_text(),time)
    replay = parse_coefficients(raw_coefficients[0].read_text(),time)
    refs = manifest['force_references'];denominator=.5*refs['rhoInf']*refs['magUInf']**2*refs['Aref']
    rows = []
    for component in ['pressure','viscous']:
        f,m = force[component+'_force'],force[component+'_moment']
        rows.append(dict(component=component,Fx_N=f[0],Fy_N=f[1],Fz_N=f[2],
                         Cd=f[0]/denominator,Cl=f[1]/denominator,Cm=m[2]/(denominator*refs['lRef'])))
    total = dict(component='total',**{key:sum(row[key] for row in rows) for key in rows[0] if key!='component'})
    rows.append(total)
    tolerance = manifest['absolute_coefficient_equivalence_tolerance']
    deltas = {key:total[key]-manifest['expected_coefficients'][key] for key in ['Cd','Cl','Cm']}
    replay_deltas = {key:replay[key]-manifest['expected_coefficients'][key] for key in ['Cd','Cl','Cm']}
    passed = all(abs(v)<=tolerance for v in [*deltas.values(),*replay_deltas.values()])
    result = dict(kind=manifest['kind'],selected_iteration=time,components=rows,
                  expected_saved_forceCoeffs=manifest['expected_coefficients'],replayed_forceCoeffs=replay,
                  component_total_minus_saved=deltas,replay_minus_saved=replay_deltas,
                  absolute_coefficient_equivalence_tolerance=tolerance,total_equivalence_passed=passed,
                  input_hashes_still_match=True,flow_iterations_executed=0,physical_validation=False,
                  input_manifest_sha256=sha256(scratch/'input-manifest.json'),
                  outputs=[dict(path=str(p.relative_to(scratch)),sha256=sha256(p))
                           for p in [*raw_forces,*raw_coefficients,scratch/'postprocess.log']])
    dump(scratch/'force-components.json',result)
    with (scratch/'force-components.csv').open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(result,indent=2))
    require(passed,'Postprocessed totals differ from saved forceCoeffs; inspect retained diagnostic without changing criteria')


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action',required=True)
    prep=sub.add_parser('prepare');prep.add_argument('--case',type=Path,required=True);prep.add_argument('--time',type=int,required=True);prep.add_argument('--out',type=Path,required=True)
    check=sub.add_parser('verify');check.add_argument('--scratch',type=Path,required=True);check.add_argument('--inputs-only',action='store_true')
    args=parser.parse_args()
    if args.action=='prepare':prepare(args.case,args.time,args.out)
    else:verify(args.scratch,args.inputs_only)


if __name__=='__main__':main()
