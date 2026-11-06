"""Prepare or check the fresh 3D finite-wing case; never starts a solver.

Foundation 14: incompressibleFluid + steadyState + SIMPLE. Uniform external
velocity is applied by patch name, so no Fluent/native node ordering is used.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import shutil

RUN = '2026-11-01-finite-wing-v003'
PATCHES = ('inlet', 'outlet', 'farfield', 'wing')
FLOAT = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
SOURCES = {
    'steady_module_controls': 'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/tutorials/incompressibleFluid/pitzDailySteady/system/fvSolution',
    'force_output_and_normalisation': 'https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/src/functionObjects/forces/forceCoeffs/forceCoeffs.C',
    'release': 'https://openfoam.org/version/14/',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def lab_root():
    return Path(__file__).resolve().parents[3]


def checked_case(path):
    path = Path(path).resolve()
    require(path == (lab_root() / 'native/runs' / RUN / 'openfoam').resolve(),
            'Only the dedicated new finite-wing-v003/openfoam case is allowed')
    return path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def write_new(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='ascii', newline='\n') as stream:
        stream.write(text.rstrip() + '\n')


def json_new(path, value):
    write_new(path, json.dumps(value, indent=2, allow_nan=False))


def copied(source, target):
    before = sha(source)
    shutil.copy2(source, target)
    require(before == sha(source) == sha(target), 'Copy hash differs: ' + str(target))


def header(name, kind='dictionary'):
    return ('FoamFile\n{\n    version 2.0;\n    format ascii;\n'
            f'    class {kind};\n    object {name};\n}}\n\n')


def dictionaries():
    control = '''solver incompressibleFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 1200;
deltaT 1;
writeControl timeStep;
writeInterval 300;
purgeWrite 0;
writeFormat ascii;
writePrecision 14;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable false;
functions
{
    forceCoefficients
    {
        type forceCoeffs;
        libs ("libforces.so");
        patches (wing);
        rho rhoInf;
        rhoInf 1;
        p p;
        U U;
        CofR (0.25 0 0);
        liftDir (0 1 0);
        dragDir (1 0 0);
        pitchAxis (0 0 1);
        magUInf 1;
        lRef 1;
        Aref 4;
        log false;
        writeControl timeStep;
        writeInterval 1;
    }
'''
    for patch in PATCHES:
        control += f'''    flux_{patch}
    {{
        type surfaceFieldValue;
        libs ("libfieldFunctionObjects.so");
        patch {patch};
        operation sum;
        fields (phi);
        writeFields false;
        log false;
        writeControl timeStep;
        writeInterval 1;
    }}
'''
    control += '}\n'
    content = {
        'system/controlDict': control,
        'system/fvSchemes': '''ddtSchemes { default steadyState; }
gradSchemes { default leastSquares; }
divSchemes
{
    default none;
    div(phi,U) bounded Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''',
        'system/fvSolution': '''solvers
{
    p { solver GAMG; tolerance 1e-10; relTol 0.01; smoother GaussSeidel; }
    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent yes;
}
relaxationFactors { equations { U 0.7; } }
''',
        'system/decomposeParDict': '''numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
''',
        'constant/physicalProperties': 'viscosityModel constant;\nnu 0.0025;\n',
        'constant/momentumTransport': 'simulationType laminar;\n',
        '0/U': '''dimensions [0 1 -1 0 0 0 0];
internalField uniform (1 0 0);
boundaryField
{
    inlet { type fixedValue; value uniform (1 0 0); }
    farfield { type fixedValue; value uniform (1 0 0); }
    outlet { type zeroGradient; }
    wing { type noSlip; }
}
''',
        '0/p': '''dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet { type zeroGradient; }
    farfield { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    wing { type zeroGradient; }
}
''',
    }
    return {name: header(Path(name).name, 'volVectorField' if name == '0/U' else
                        'volScalarField' if name == '0/p' else 'dictionary') + body
            for name, body in content.items()}


def clean(text):
    return re.sub(r'/\*.*?\*/|//[^\n]*', '', text, flags=re.S)


def counted(path):
    text = clean(path.read_text(encoding='ascii'))
    h = re.search(r'FoamFile\s*\{([^{}]*)\}', text)
    require(h is not None and re.search(r'\bformat\s+ascii\s*;', h[1]), 'Require native ASCII list: ' + str(path))
    m = re.fullmatch(r'\s*(\d+)\s*\((.*)\)\s*', text[h.end():], re.S)
    require(m is not None, 'Invalid counted native list: ' + str(path))
    return int(m[1]), m[2]


def native_poly(case):
    candidates = [case / p for p in ('0/polyMesh', 'constant/polyMesh') if (case / p / 'points').is_file()]
    require(len(candidates) == 1, 'Require exactly one initial native polyMesh')
    return candidates[0]


def native_patches(poly):
    count, body = counted(poly / 'boundary')
    rows = {}
    for name, block in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}', body):
        name = name.strip('"')
        rows[name] = dict(kind=re.search(r'\btype\s+(\w+)\s*;', block)[1],
                          count=int(re.search(r'\bnFaces\s+(\d+)\s*;', block)[1]),
                          start=int(re.search(r'\bstartFace\s+(\d+)\s*;', block)[1]))
    require(len(rows) == count and set(rows) == set(PATCHES), 'Converted 3D patches differ from exact source names')
    require(rows['wing']['kind'] == 'wall' and all(v['count'] > 0 and v['kind'] != 'empty' for v in rows.values()),
            'Wing must be wall, all patches nonempty and none an empty 2D plane')
    ranges = sorted((v['start'], v['start'] + v['count']) for v in rows.values())
    require(all(a[1] == b[0] for a, b in zip(ranges, ranges[1:])), 'Native boundary face ranges are not contiguous')
    return rows


def mesh_quality(log):
    require('checkMesh -allGeometry -allTopology' in log and 'Mesh OK.' in log and
            not re.search(r'Failed\s+\d+\s+mesh checks|FOAM FATAL', log), 'Full native checkMesh did not pass')
    names = {'cells': r'\bcells:\s+(\d+)', 'minimum_volume_m3': r'Min volume\s*=\s*(' + FLOAT + ')',
             'maximum_nonorthogonality_degrees': r'Mesh non-orthogonality Max:\s*(' + FLOAT + ')',
             'maximum_skewness': r'Max skewness\s*=\s*(' + FLOAT + ')'}
    result = {}
    for name, pattern in names.items():
        m = re.search(pattern, log)
        require(m is not None, 'Missing native mesh quality value: ' + name)
        result[name] = int(m[1]) if name == 'cells' else float(m[1])
    require(result['minimum_volume_m3'] > 0 and result['maximum_nonorthogonality_degrees'] <= 75
            and result['maximum_skewness'] <= 4, 'Native quality gate missed: positive volumes, nonorth<=75, skew<=4')
    return result


def prepare(mesh, case):
    mesh, case = mesh.resolve(), checked_case(case)
    require(not case.exists(), 'Fresh native case required; existing projects cannot be overwritten')
    require(mesh.parent == case.parent / 'mesh-export-v002' and mesh.name == 'wing.msh', 'Require explicitly selected repaired mesh-export-v002/wing.msh')
    raw = mesh.read_bytes()
    require(re.search(rb'\(2\s+3\)', raw), 'Require 3D Fluent ASCII mesh, never -2D extrusion')
    cell_match = re.search(rb'^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)', raw, re.M)
    require(cell_match is not None, 'Missing shared Fluent cell count')
    cells = int(cell_match[1], 16)
    require(0 < cells <= 250000, '3D first-model cap is 250000 cells')
    for name in ('fluid', 'interior', *PATCHES):
        require(re.search(rb'\(45\s+\([^)]*\s' + name.encode() + rb'\)\(\)\)', raw), 'Missing named Fluent zone ' + name)
    require(shutil.disk_usage(case.anchor).free >= 10 * 1024**3, 'Require 10 GiB free output drive space')
    setup = json.loads((mesh.parent / 'setup.json').read_text(encoding='utf-8-sig'))
    require(setup['cells']==cells and setup['mesh_sha256']==sha(mesh), 'Shared setup and native mesh hashes/counts differ')
    require(setup['chord_m']==1 and setup['span_m']==4 and setup['reference_area_m2']==4
            and setup['incidence_degrees_clockwise']==8 and setup['Re_c']==400
            and setup['density_kg_m3']==1 and setup['kinematic_viscosity_m2_s']==.0025,
            'Shared geometry/physical setup differs from the approved model')
    patch_counts = {row['name']:row['faces'] for row in setup['boundary_zones'].values() if row['name'] in PATCHES}
    require(set(patch_counts)==set(PATCHES), 'Shared setup is missing physical patch counts')
    case.mkdir(parents=True)
    for name, text in dictionaries().items():
        write_new(case / name, text)
    copied(mesh, case / 'input-wing.msh')
    provenance = case / 'provenance'
    provenance.mkdir()
    for source in (Path(__file__), lab_root() / 'local/finite-wing-v003/run-openfoam.sh',
                   lab_root() / 'local/finite-wing-v003/DESIGN.md'):
        copied(source, provenance / source.name)
    copied(mesh.parent / 'setup.json', provenance / 'mesh-setup.json')
    declared = dict(model='Full-span 3D rectangular finite wing, illustrative laminar steady calculation',
        x_bounds_m=[-3, 9], y_bounds_m=[-4, 4], z_bounds_m=[-6, 6],
        chord_m=1, full_span_m=4, planform_area_m2=4, aspect_ratio=4, incidence_degrees=8,
        inlet_U_m_s=[1, 0, 0], density_kg_m3=1, dynamic_viscosity_Pa_s=.0025,
        kinematic_viscosity_m2_s=.0025, reynolds_chord=400,
        wing='Stationary no-slip exact source wing patch; shape/pivot in mesh setup.json',
        farfield='Uniform fixed U(1,0,0) on upstream and four lateral faces; zeroGradient p',
        outlet='zeroGradient U and fixed kinematic p=0', initial_U_m_s=[1,0,0], initial_p=0,
        solver='OpenFOAM Foundation14 incompressibleFluid', algorithm='steady SIMPLE consistent',
        max_iterations=1200, saved_iterations=[300,600,900,1200], mpi_processes=2,
        U_equation_relaxation=.7, nonorthogonal_correctors=1,
        force_reference_area_m2=4, force_reference_length_m=1, force_reference_U_m_s=1,
        force_reference_rho_kg_m3=1, drag_direction=[1,0,0], lift_direction=[0,1,0],
        moment_reference_m=[.25,0,0], force_sample_interval_iterations=1,
        residual_target=1e-5, relative_mass_imbalance_target=1e-4,
        last100_relative_Cd_Cl_range_target=1e-3, relative_force_denominator='max(abs(last100 mean),1e-3)',
        mesh_quality_targets=dict(positive_volume=True,max_nonorthogonality_degrees=75,max_skewness=4,full_checkMesh_pass=True),
        intended_cells=cells, source_boundary_face_counts=patch_counts,
        mesh_source_path=str(mesh), mesh_source_sha256=sha(mesh),
        physical_validation=False, mesh_independence_established=False,
        residual_norms_equivalent_to_fluent=False, pressure_note='Kinematic p; gauge Pa=rho*p',
        limitation='A fixed steady solve may fail force stability at this Reynolds number; preserve misses; do not infer a steady physical wake',
        official_sources=SOURCES, solver_started=False)
    json_new(case / 'declared-comparison.json', declared)
    json_new(provenance / 'prepared-input-sha256.json',
             {p.relative_to(case).as_posix(): sha(p) for p in sorted(case.rglob('*')) if p.is_file()})
    json_new(case / 'run-status.json', dict(stage='prepared_only', solver_started=False))
    return dict(case=str(case), cells=cells, mesh_sha256=sha(mesh), solver_started=False)


def validate_native(case):
    case = checked_case(case)
    require(not (case / 'solver.log').exists() and not list(case.glob('processor*')), 'Native preflight must precede any solver/decomposition')
    manifest = json.loads((case / 'provenance/prepared-input-sha256.json').read_text())
    require(all(sha(case / name) == digest for name,digest in manifest.items()), 'Prepared inputs changed before conversion/solve')
    poly = native_poly(case)
    patches = native_patches(poly)
    quality = mesh_quality((case / 'mesh-check.log').read_text(encoding='utf-8'))
    declared = json.loads((case / 'declared-comparison.json').read_text())
    require(quality['cells'] == declared['intended_cells'], 'Converted native cell count differs')
    require({name:row['count'] for name,row in patches.items()}==declared['source_boundary_face_counts'],
            'Converted native patch face counts differ from named shared source zones')
    count, body = counted(poly / 'points')
    points = [tuple(map(float, row.split())) for row in re.findall(r'\(([^()]*)\)', body)]
    require(len(points) == count and all(len(p)==3 and all(math.isfinite(v) for v in p) for p in points), 'Invalid native 3D points')
    bounds = [[min(p[j] for p in points), max(p[j] for p in points)] for j in range(3)]
    require(all(abs(a-b)<=1e-8 for row,target in zip(bounds, ([-3,9],[-4,4],[-6,6])) for a,b in zip(row,target)), 'Native box bounds differ')
    faces, body = counted(poly / 'faces')
    sizes = re.findall(r'(\d+)\s*\([^()]*\)', body)
    require(len(sizes)==faces and all(int(v)==3 for v in sizes), 'Expected linear tetrahedral native triangular faces')
    receipt = dict(passed=True, actual_native_patches=patches, actual_mesh_quality=quality,
        native_bounds_m=bounds, native_points=count, native_faces=faces, uniform_patch_values_require_no_node_ordering=True,
        native_mesh_sha256={name:sha(poly/name) for name in ('points','faces','owner','neighbour','boundary')},
        mesh_location=poly.relative_to(case).as_posix(), prepared_inputs_unchanged=True, solver_started=False)
    json_new(case / 'provenance/preflight-validation.json', receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--mesh', type=Path)
    group.add_argument('--validate-native', action='store_true')
    args = parser.parse_args()
    print(json.dumps(prepare(args.mesh,args.case) if args.mesh else validate_native(args.case), indent=2))


if __name__ == '__main__':
    main()
