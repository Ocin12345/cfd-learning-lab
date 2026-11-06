"""Prepare/check the shared Re=100 transient cylinder case for Foundation 14.

Preparation never runs OpenFOAM. --apply-inlet is used only after native mesh
conversion; it reads inlet faces in polyMesh order, never Fluent node order.
All persistent/native output belongs to the external CFD-Learning-Lab tree.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import re
import shutil

RUN = "2026-11-01-vortex-shedding-v001"
HEIGHT = .41


def lab_root():
    return Path(__file__).resolve().parents[3]


def checked_case(case):
    case = Path(case).resolve()
    run = lab_root() / "native" / "runs" / RUN
    allowed = {(run / name).resolve() for name in ("openfoam", "openfoam-extended")}
    if case not in allowed:
        raise ValueError("Use only the dedicated vortex-shedding-v001/openfoam or openfoam-extended case")
    return case


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="ascii", newline="\n") as stream:
        stream.write(text.rstrip() + "\n")


def json_new(path, data):
    write_new(path, json.dumps(data, indent=2))


def foam_header(name, kind="dictionary"):
    return ("FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            f"    class {kind};\n    object {name};\n}}\n\n")


def velocity_field(inlet):
    return foam_header("U", "volVectorField") + """
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    inlet
    {
        type fixedValue;
        VALUE
    }
    outlet { type zeroGradient; }
    walls { type noSlip; }
    cylinder { type noSlip; }
    "(frontAndBack.*|defaultFaces)" { type empty; }
}
""".replace("VALUE", inlet).strip() + "\n"


def control_body(end=2, latest=False):
    body = f"""
solver incompressibleFluid;
startFrom {'latestTime' if latest else 'startTime'};
startTime 0;
stopAt endTime;
endTime {end};
deltaT 0.005;
adjustTimeStep false;
writeControl timeStep;
writeInterval 5;
purgeWrite 0;
writeFormat ascii;
writePrecision 14;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable false;
functions
{{
    forceCoefficients
    {{
        type forceCoeffs;
        libs ("libforces.so");
        patches (cylinder);
        rho rhoInf;
        rhoInf 1;
        p p;
        U U;
        CofR (0.2 0.2 0);
        liftDir (0 1 0);
        dragDir (1 0 0);
        pitchAxis (0 0 1);
        magUInf 1;
        lRef 0.1;
        Aref 0.1;
        log false;
        writeControl timeStep;
        writeInterval 1;
    }}
"""
    for patch in ("inlet", "outlet", "walls", "cylinder"):
        body += f"""
    flux_{patch}
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
"""
    return body + "}\n"


def prepare(mesh, case):
    mesh, case = Path(mesh).resolve(), checked_case(case)
    if case.exists():
        raise FileExistsError("Refusing to overwrite any existing native case")
    content = mesh.read_bytes()
    if mesh.name != "vortex.msh" or not re.search(rb"\(2\s+2\)", content):
        raise ValueError("Expected the new shared 2D Fluent ASCII vortex.msh")
    cells = re.search(rb"^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)", content, re.M)
    if cells is None or int(cells[1], 16) <= 0:
        raise ValueError("Missing Fluent cell count")
    for zone in (b"fluid", b"inlet", b"outlet", b"walls", b"cylinder"):
        if not re.search(rb"\(45\s+\([^)]*\s" + zone + rb"\)\(\)\)", content):
            raise ValueError("Missing expected Fluent zone: " + zone.decode())
    if shutil.disk_usage(case.parent).free < 10 * 1024**3:
        raise RuntimeError("Require 10 GiB free on the output drive")
    local = lab_root() / "local" / "vortex-shedding-v001"
    wrapper = local / "run-openfoam.sh"
    if not wrapper.is_file():
        raise FileNotFoundError(wrapper)
    case.mkdir(parents=True)

    def write(relative, body, kind="dictionary"):
        path = case / relative
        write_new(path, foam_header(path.name, kind) + body.strip())

    write("system/controlDict", control_body())
    write("system/fvSchemes", """
ddtSchemes { default backward; }
gradSchemes { default leastSquares; }
divSchemes
{
    default none;
    div(phi,U) Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    write("system/fvSolution", """
solvers
{
    p
    {
        solver GAMG;
        tolerance 1e-9;
        relTol 0.01;
        smoother GaussSeidel;
    }
    pFinal { $p; relTol 0; }
    U
    {
        solver smoothSolver;
        smoother symGaussSeidel;
        tolerance 1e-10;
        relTol 0.01;
    }
    UFinal { $U; relTol 0; }
}
PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 3;
    nCorrectors 2;
    nNonOrthogonalCorrectors 1;
}
relaxationFactors { equations { ".*" 1; } }
""")
    write("system/decomposeParDict", """
numberOfSubdomains 2;
method simple;
simpleCoeffs { n (2 1 1); delta 0.001; }
""")
    write("constant/physicalProperties", "viscosityModel constant;\nnu 0.001;")
    write("constant/momentumTransport", "simulationType laminar;")
    write_new(case / "0/U", velocity_field("value uniform (0 0 0);"))
    write("0/p", """
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    walls { type zeroGradient; }
    cylinder { type zeroGradient; }
    "(frontAndBack.*|defaultFaces)" { type empty; }
}
""", "volScalarField")
    shutil.copy2(mesh, case / "input-vortex.msh")
    provenance = case / "provenance"
    provenance.mkdir()
    for source in (Path(__file__), wrapper, local / "DESIGN.md"):
        if source.is_file():
            shutil.copy2(source, provenance / source.name)
    for name in ("setup.json", "inlet-face-profile.json", "parabolic-inlet.prof"):
        source = mesh.parent / name
        if source.is_file():
            shutil.copy2(source, provenance / ("mesh-" + name))
    targets = dict(
        model="2D transient incompressible laminar confined cylinder; DFG-inspired",
        benchmark_validation=False, physical_validation=False,
        x_bounds_m=[0, 2.2], y_bounds_m=[0, HEIGHT],
        cylinder_centre_m=[.2, .2], cylinder_radius_m=.05,
        extrusion_depth_m=1, density_kg_m3=1,
        kinematic_viscosity_m2_s=.001, dynamic_viscosity_Pa_s=.001,
        reynolds_diameter_mean_velocity=100,
        inlet_ux_expression="6*y*(0.41-y)/0.41^2", inlet_uy=0,
        inlet_continuum_mean_m_s=1, inlet_continuum_max_m_s=1.5,
        inlet_assignment="Native converted polyMesh inlet-face centre and patch order",
        inlet_discrete_values_rescaled=False,
        initial_U_m_s=[0, 0, 0], initial_kinematic_p=0,
        startup="Abrupt parabolic inlet imposed on initially stationary fluid",
        outlet_gauge_pressure_Pa=0,
        mesh_source_path=str(mesh), mesh_source_sha256=sha(mesh),
        intended_cells=int(cells[1], 16), mpi_processes=2,
        solver="OpenFOAM Foundation 14 incompressibleFluid",
        time_scheme="backward", momentum_convection="Gauss linearUpwind grad(U)",
        delta_t_s=.005, initial_pilot_end_time_s=2,
        continuation_end_time_s=8, continuation_requires_root_authorization=True,
        saved_field_interval_s=.025, force_write_interval_s=.005,
        pimple_outer_correctors=3, pressure_correctors=2,
        force_reference_U_m_s=1, force_reference_area_m2=.1,
        force_reference_length_m=.1, force_CofR_m=[.2, .2, 0],
        residual_norms_equivalent_to_fluent=False,
        pressure_note="Kinematic p [m2/s2]; Pa gauge = rho*p",
        relative_mass_imbalance_target=1e-5,
        inlet_profile_applied=False, solver_started=False)
    json_new(case / "declared-comparison.json", targets)
    hashes = {p.relative_to(case).as_posix(): sha(p)
              for p in sorted(case.rglob("*")) if p.is_file()}
    json_new(provenance / "prepared-input-sha256.json", hashes)
    shutil.copytree(case / "system", provenance / "initial-system")
    shutil.copytree(case / "0", provenance / "initial-0-before-inlet")
    json_new(case / "run-status.json", dict(stage="prepared", solver_started=False,
                                           inlet_profile_applied=False))
    print(json.dumps(dict(case=str(case), intended_cells=int(cells[1], 16),
                          mesh_source_sha256=sha(mesh), solver_started=False), indent=2))


def foam_list(path):
    text = Path(path).read_text(encoding="ascii")
    text = re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)
    header = re.search(r"FoamFile\s*\{([^{}]*)\}", text)
    if not header or not re.search(r"\bformat\s+ascii\s*;", header[1]):
        raise ValueError("Require native ASCII polyMesh: " + str(path))
    match = re.fullmatch(r"\s*(\d+)\s*\((.*)\)\s*", text[header.end():], re.S)
    if not match:
        raise ValueError("Invalid counted OpenFOAM list: " + str(path))
    return int(match[1]), match[2]


def native_mesh(case):
    candidates = [case / "0/polyMesh", case / "constant/polyMesh"]
    present = [p for p in candidates if (p / "points").is_file()]
    if len(present) != 1:
        raise ValueError("Require exactly one native initial mesh in 0 or constant")
    return present[0]


def native_inlet_profile(poly):
    count, body = foam_list(poly / "points")
    points = [tuple(map(float, row.split())) for row in re.findall(r"\(([^()]*)\)", body)]
    if len(points) != count or any(len(p) != 3 for p in points):
        raise ValueError("Invalid native point count")
    count, body = foam_list(poly / "faces")
    faces = []
    for size, labels in re.findall(r"(\d+)\s*\(([^()]*)\)", body):
        face = list(map(int, labels.split()))
        if len(face) != int(size) or any(i < 0 or i >= len(points) for i in face):
            raise ValueError("Invalid native face")
        faces.append(face)
    if len(faces) != count:
        raise ValueError("Invalid native face count")
    count, body = foam_list(poly / "boundary")
    patches = {name.strip('"'): data for name, data in
               re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}', body)}
    if len(patches) != count:
        raise ValueError("Invalid patch count")
    for name in ("inlet", "outlet", "walls", "cylinder"):
        if name not in patches:
            raise ValueError("Missing native patch " + name)
    data = patches["inlet"]
    start = int(re.search(r"\bstartFace\s+(\d+)\s*;", data)[1])
    nfaces = int(re.search(r"\bnFaces\s+(\d+)\s*;", data)[1])
    if nfaces <= 0 or start + nfaces > len(faces):
        raise ValueError("Invalid inlet face range")
    rows = []
    for index in range(start, start + nfaces):
        vertices = [points[i] for i in faces[index]]
        # fluentMeshToFoam -2D extrudes inlet edges into planar rectangles.
        # A rectangle's vertex average is precisely OpenFOAM face::centre.
        if len(vertices) != 4 or any(abs(v[0]) > 1e-10 for v in vertices):
            raise ValueError("Expected rectangular extruded native inlet faces on x=0")
        ys, zs = sorted({v[1] for v in vertices}), sorted({v[2] for v in vertices})
        if len(ys) != 2 or len(zs) != 2 or set((v[1], v[2]) for v in vertices) != {(y,z) for y in ys for z in zs}:
            raise ValueError("Native inlet face is not a rectangle")
        centre = [sum(v[j] for v in vertices) / 4 for j in range(3)]
        y = centre[1]
        if not 0 < y < HEIGHT or abs((zs[1]-zs[0])-1) > 1e-10:
            raise ValueError("Unexpected native inlet height or span")
        area = (ys[1]-ys[0]) * (zs[1]-zs[0])
        ux = 6*y*(HEIGHT-y)/HEIGHT**2
        rows.append(dict(patch_face_index=index-start, native_face_index=index,
                         centre_m=centre, area_m2=area, U_m_s=[ux, 0, 0]))
    return rows


def apply_inlet(case):
    case = checked_case(case)
    receipt = case / "provenance/native-inlet-profile.json"
    if receipt.exists() or (case / "solver.log").exists() or list(case.glob("processor*")):
        raise FileExistsError("Inlet application is fresh-case only and cannot be repeated")
    prepared = json.loads((case / "provenance/prepared-input-sha256.json").read_text())
    if sha(case / "0/U") != prepared["0/U"]:
        raise ValueError("Initial U changed before native inlet application")
    poly = native_mesh(case)
    rows = native_inlet_profile(poly)
    values = "\n".join("            (%.16g 0 0)" % row["U_m_s"][0] for row in rows)
    value = f"value nonuniform List<vector>\n        {len(rows)}\n        (\n{values}\n        );"
    (case / "0/U").write_text(velocity_field(value), encoding="ascii", newline="\n")
    area = sum(r["area_m2"] for r in rows)
    flux = sum(r["U_m_s"][0]*r["area_m2"] for r in rows)
    if abs(area-HEIGHT) > 1e-10:
        raise ValueError("Native inlet area does not equal 0.41 m2")
    data = dict(order="Native polyMesh inlet startFace + patch-local index",
                face_centres="Exact rectangle centres in converted extruded mesh",
                ux_expression="6*y*(0.41-y)/0.41^2", rescaled=False,
                n_faces=len(rows), inlet_area_m2=area,
                discrete_volume_flow_m3_s=flux, continuum_volume_flow_m3_s=HEIGHT,
                discrete_flux_error_percent=100*(flux/HEIGHT-1),
                applied_U_sha256=sha(case / "0/U"),
                native_mesh_location=poly.relative_to(case).as_posix(),
                applicator_sha256=sha(Path(__file__)),
                native_mesh_sha256={name: sha(poly / name)
                                    for name in ("points", "faces", "boundary")},
                faces=rows)
    json_new(receipt, data)
    shutil.copy2(Path(__file__), case / "provenance/native-profile-applicator.py")
    shutil.copy2(case / "0/U", case / "provenance/initial-U-with-native-inlet")
    print(json.dumps({key: val for key, val in data.items() if key != "faces"}, indent=2))


def continuation(case, authorization):
    case = checked_case(case)
    auth = json.loads(Path(authorization).read_text(encoding="utf-8"))
    if auth.get("approved") is not True or auth.get("action") != "continue_openfoam_vortex_to_8s" or auth.get("approved_by") != "root":
        raise ValueError("Require explicit root authorization for continuation to 8 seconds")
    if not (case / "2/U").is_file() or not (case / "2/p").is_file():
        raise ValueError("Require reconstructed completed 2-second pilot")
    target = case / "provenance/continuation-to-8s"
    if target.exists():
        raise FileExistsError("Continuation evidence already exists")
    original = case / "provenance/initial-system/controlDict"
    if sha(case / "system/controlDict") != sha(original):
        raise ValueError("Pilot controlDict changed; preserve and inspect it before continuation")
    target.mkdir()
    shutil.copy2(authorization, target / "root-authorization.json")
    shutil.copy2(case / "system/controlDict", target / "pilot-controlDict")
    body = foam_header("controlDict") + control_body(8, latest=True).strip() + "\n"
    write_new(target / "continuation-controlDict", body)
    (case / "system/controlDict").write_text(body, encoding="ascii", newline="\n")


def write_force_history(case, end):
    files = sorted((case / "postProcessing/forceCoefficients").glob("*/forceCoeffs.dat"))
    samples = {}
    for path in files:
        header = None
        for line in path.read_text().splitlines():
            words = line.lstrip("# ").split()
            if line.lstrip().startswith("#"):
                if "Cd" in words and "Cl" in words and "Time" in words:
                    header = words
                continue
            if not words:
                continue
            if header is None:
                raise ValueError("Missing native Cd/Cl column header: " + str(path))
            values = [float(word) for word in words]
            if len(values) != len(header):
                raise ValueError("Unexpected native force coefficient column count")
            row = dict(zip(header, values))
            time = row["Time"]
            if 0 <= time <= end + 1e-8:
                if not all(math.isfinite(row[name]) for name in ("Time", "Cd", "Cl")):
                    raise ValueError("Non-finite native force coefficient")
                samples[time] = (row["Cd"], row["Cl"])
    if not samples or abs(max(samples)-end) > 1e-8:
        raise ValueError("Missing final native force coefficient sample")
    output = case / "force-history.csv"
    if output.exists():
        if end != 8:
            raise FileExistsError("Force history already exists")
        shutil.copy2(output, case / "provenance/force-history-before-8s.csv")
    with output.open("w", encoding="ascii", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("time_s", "Cd", "Cl"))
        for time in sorted(samples):
            writer.writerow((format(time, ".16g"), *(format(v, ".16g") for v in samples[time])))
    return dict(csv="force-history.csv", rows=len(samples),
                source_files=[p.relative_to(case).as_posix() for p in files],
                source_sha256={p.relative_to(case).as_posix(): sha(p) for p in files})

def summarize(case, end):
    case = checked_case(case)
    fluxes = {}
    for patch in ("inlet", "outlet", "walls", "cylinder"):
        samples = []
        for path in (case / "postProcessing" / ("flux_" + patch)).glob("*/*.dat"):
            for line in path.read_text().splitlines():
                if line.strip() and not line.lstrip().startswith("#"):
                    values = [float(x) for x in line.split()]
                    if len(values) == 2 and abs(values[0]-end) < 1e-8:
                        samples.append(values[1])
        if not samples or not all(math.isfinite(v) for v in samples):
            raise ValueError("Missing finite final native flux for " + patch)
        fluxes[patch] = samples[-1]
    closure = abs(sum(fluxes.values())) / max(abs(fluxes["inlet"]), 1e-30)
    force_history = write_force_history(case, end)
    summary = dict(end_time_s=end, density_kg_m3=1,
                   signed_patch_volume_flow_m3_s=fluxes,
                   relative_mass_imbalance=closure,
                   relative_mass_imbalance_target=1e-5,
                   mass_target_passed=closure <= 1e-5,
                   execution_finished=True, physical_validation=False,
                   paired_comparison_pending=True, force_history=force_history)
    json_new(case / ("numerical-summary-pilot.json" if end == 2 else "numerical-summary-8s.json"), summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, type=Path)
    commands = parser.add_mutually_exclusive_group(required=True)
    commands.add_argument("--mesh", type=Path)
    commands.add_argument("--apply-inlet", action="store_true")
    commands.add_argument("--prepare-continuation", action="store_true")
    commands.add_argument("--summarize", type=float, choices=(2, 8))
    parser.add_argument("--authorization", type=Path)
    args = parser.parse_args()
    if args.mesh:
        prepare(args.mesh, args.case)
    elif args.apply_inlet:
        apply_inlet(args.case)
    elif args.prepare_continuation:
        if args.authorization is None:
            parser.error("--prepare-continuation requires --authorization")
        continuation(args.case, args.authorization)
    else:
        summarize(args.case, args.summarize)
