"""Prepare exact localEuler or original-steady restart inputs for the formula car."""
import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

STAGES = {'local': (600, 1500), 'steady': (1500, 1800), 'blended': (600, 1200)}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def once(data, before, after):
    assert data.count(before) == 1, before
    return data.replace(before, after)


def configured(original, phase):
    start, end = STAGES[phase]
    control = once(original['controlDict'], b'startTime 0;', f'startTime {start};'.encode())
    control = once(control, b'endTime 600;', f'endTime {end};'.encode())
    result = {**original, 'controlDict': control}
    if phase == 'local':
        result['fvSchemes'] = once(original['fvSchemes'], b'default         steadyState;', b'default         localEuler;')
        solution = once(original['fvSolution'], b'    Phi\n', b'    pFinal\n    {\n        $p;\n    }\n\n    Phi\n')
        for field in ('U', 'k', 'omega'):
            pattern = rb'(?m)^    '+field.encode()+rb'\n    \{[^}]+\}'
            block, = re.findall(pattern, solution)
            solution = once(solution, block, block+f'\n\n    {field}Final\n    {{\n        ${field};\n    }}'.encode())
        old = b'SIMPLE\n{\n    nNonOrthogonalCorrectors 0;\n    consistent yes;\n}'
        new = b'''PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors 2;
    nNonOrthogonalCorrectors 0;
    maxCo 5;
    maxDeltaT 0.1;
    rDeltaTSmoothingCoeff 0.1;
}'''
        result['fvSolution'] = once(solution, old, new)
    if phase == 'blended':
        result['fvSchemes'] = once(
            original['fvSchemes'], b'bounded Gauss linearUpwindV grad(U);',
            b'bounded Gauss fixedBlended 0.25 linearUpwindV grad(U) upwind;'
        )
    return result


def complete(case, end, first):
    log = (case/'solver.log').read_text(encoding='ascii')
    assert log.rstrip().endswith('End')
    assert [int(i) for i in re.findall(r'^Time = (\d+)s?\s*$',log,re.M)] == list(range(first,end+1))
    clock = (case/str(end)/'uniform/time').read_text(encoding='ascii')
    assert re.search(rf'^index\s+{end};',clock,re.M)
    for name in ('U', 'p', 'phi', 'k', 'omega', 'nut'):
        assert (case/str(end)/name).stat().st_size > 0


def setup(parent, source, output, phase, original_audit):
    start, end = STAGES[phase]
    complete(parent, 600, 1)
    audit = json.loads(original_audit.read_bytes())
    assert audit['case'] == parent.name and audit['iterations'] == 600
    required = {name: sha for name,sha in audit['source_sha256'].items() if name.startswith(('0/','600/','constant/','system/'))}
    assert required and all(digest(parent/name)==sha for name,sha in required.items())
    if phase in ('local','blended'):
        assert source == parent
    else:
        complete(source, start, 601)
        assert (source/'input-sha256.txt').is_file()
        for line in (source/'input-sha256.txt').read_text().splitlines():
            sha,name=line.split('  ',1)
            assert digest(source/name)==sha,name
        for name in required:
            if name.startswith(('0/','constant/')):
                assert digest(source/name)==required[name],name
    output.mkdir()
    for name in ('0','constant','system'):
        shutil.copytree(parent/name,output/name)
    shutil.copytree(source/str(start),output/str(start))
    origin = {}
    for item in sorted(output.rglob('*')):
        if item.is_file():
            rel = item.relative_to(output).as_posix()
            src = (source if rel.startswith(f'{start}/') else parent)/rel
            origin[rel] = digest(src)
            assert digest(item)==origin[rel],rel
    controls = ('controlDict','fvSchemes','fvSolution')
    original = {name:(parent/'system'/name).read_bytes() for name in controls}
    for name,data in configured(original,phase).items():
        (output/'system'/name).write_bytes(data)
    receipt = {'phase':phase,'start':start,'end':end,'parent':parent.name,
               'checkpoint_case':source.name,'original_audit_sha256':digest(original_audit),'copied':origin}
    with (output/'restart-provenance.json').open('x',encoding='utf-8') as stream:
        json.dump(receipt,stream,indent=2)
    with (output/'input-sha256.txt').open('x',encoding='ascii') as stream:
        for name in sorted(origin):
            stream.write(f'{digest(output/name)}  {name}\n')


if __name__ == '__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('phase',choices=STAGES)
    for name in ('parent','source','output','original_audit'):
        cli.add_argument(name,type=lambda text:Path(text).resolve())
    args=cli.parse_args()
    setup(args.parent,args.source,args.output,args.phase,args.original_audit)
