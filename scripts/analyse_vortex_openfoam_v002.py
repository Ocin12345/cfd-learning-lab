"""Read-only numerical checks for the actual Re=100 OpenFOAM 2/4-second transient.

--run is the native OpenFOAM folder, not its parent. --out is a fresh directory.
Only existing solver evidence is analysed; this program never starts a solver.
Terminal initial residual means the last printed linear-system initial residual
for each variable in a physical time step, after the PIMPLE outer corrections.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import re
import sys
import zipfile

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
PHYSICAL_PATCHES = ("inlet", "outlet", "walls", "cylinder")
RESIDUAL_VARIABLES = ("Ux", "Uy", "p")
RESIDUAL_TARGET = 1e-5
MASS_TARGET = 1e-4


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean(text):
    return re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)


def near(a, b, absolute=1e-10, relative=1e-10):
    return abs(a-b) <= absolute + relative * max(abs(a), abs(b))


class Sources:
    def __init__(self, run):
        self.run = run
        self.before = {}

    def text(self, path):
        path = Path(path).resolve()
        require(path.is_file(), "Missing actual evidence: " + str(path))
        if path not in self.before:
            self.before[path] = sha(path)
        return path.read_text(encoding="utf-8")

    def local(self, relative):
        return self.text(self.run / relative)

    def json(self, relative):
        return json.loads(self.local(relative))

    def finish(self):
        sources = []
        unchanged = True
        for path, before in self.before.items():
            after = sha(path) if path.is_file() else None
            unchanged = unchanged and before == after
            sources.append(dict(path=str(path), bytes=path.stat().st_size if path.is_file() else None,
                                sha256_before=before, sha256_after=after))
        return dict(inputs_unchanged_during_analysis=unchanged,
                    native_fields_or_logs_copied_to_repository=False,
                    sources=sources)


def scalar_setting(text, name):
    match = re.search(r"\b" + re.escape(name) + r"\s+(" + FLOAT + r")\s*;", text)
    require(match is not None, "Missing actual scalar setting " + name)
    value = float(match[1])
    require(math.isfinite(value), "Non-finite actual setting " + name)
    return value


def flat_block(text, name):
    match = re.search(r"\b" + re.escape(name) + r"\s*\{([^{}]*)\}", text, re.S)
    require(match is not None, "Missing native dictionary block " + name)
    return match[1]


def numeric_values(raw, components):
    tokens = raw.replace("(", " ").replace(")", " ").split()
    values = [float(token) for token in tokens]
    require(all(math.isfinite(v) for v in values), "Non-finite native field value")
    require(len(values) % components == 0, "Native field component count differs")
    return values


def native_field(raw, count, components, dimensions):
    text = clean(raw)
    require(re.search(r"\bformat\s+ascii\s*;", text), "Only native ASCII field evidence is supported")
    dim = re.search(r"\bdimensions\s*\[([^]]+)\]\s*;", text)
    require(dim and list(map(int, dim[1].split())) == dimensions, "Native field dimensions differ")
    # Catch non-finite tokens anywhere, including boundary values, before parsing.
    require(not re.search(r"(?i)(?<![A-Za-z0-9_])[-+]?(?:nan|inf(?:inity)?)(?![A-Za-z0-9_])", text),
            "Non-finite token in native field or boundary")
    value = re.search(r"\binternalField\s+uniform\s+([^;]+);", text)
    if value:
        values = numeric_values(value[1], components)
        require(len(values) == components, "Invalid uniform native field")
        minimum, maximum = values, values
        actual_count = count
        max_abs_z = abs(values[2]) if components == 3 else None
    else:
        kind = "vector" if components == 3 else "scalar"
        value = re.search(r"\binternalField\s+nonuniform\s+List<" + kind +
                          r">\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
        require(value is not None and int(value[1]) == count, "Native internal field element count differs")
        values = numeric_values(value[2], components)
        require(len(values) == count*components, "Native internal field data length differs")
        minimum = [min(values[i::components]) for i in range(components)]
        maximum = [max(values[i::components]) for i in range(components)]
        actual_count = len(values)//components
        max_abs_z = max(abs(v) for v in values[2::3]) if components == 3 else None
    return dict(finite=True, elements=actual_count, components=components,
                minimum=minimum, maximum=maximum, maximum_absolute_z=max_abs_z)


def native_boundary(raw):
    text = clean(raw)
    patches = {}
    for name, data in re.findall(r"\b(\w+)\s*\{([^{}]*)\}", text):
        count = re.search(r"\bnFaces\s+(\d+)\s*;", data)
        if count:
            patches[name] = dict(count=int(count[1]),
                                 start=int(re.search(r"\bstartFace\s+(\d+)\s*;", data)[1]),
                                 kind=re.search(r"\btype\s+(\w+)\s*;", data)[1])
    require(set(patches) == set(PHYSICAL_PATCHES) | {"frontAndBackPlanes"}, "Unexpected native patch set")
    require(patches["frontAndBackPlanes"]["kind"] == "empty", "Front/back must be empty")
    require(patches["walls"]["kind"] == patches["cylinder"]["kind"] == "wall", "Expected native wall types")
    return patches


def boundary_mass(raw_phi, patches):
    text = clean(raw_phi)
    require("boundaryField" in text, "Missing native phi boundaryField")
    blocks = dict(re.findall(r"\b(\w+)\s*\{([^{}]*)\}", text.split("boundaryField",1)[1]))
    flow = {}
    for name in PHYSICAL_PATCHES:
        require(name in blocks, "Missing phi boundary " + name)
        data = blocks[name]
        uniform = re.search(r"\bvalue\s+uniform\s+([^;]+);", data)
        if uniform:
            value = float(uniform[1])
            require(math.isfinite(value), "Non-finite uniform boundary phi")
            flow[name] = value * patches[name]["count"]
        else:
            match = re.search(r"\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",data,re.S)
            require(match is not None and int(match[1]) == patches[name]["count"], "Boundary phi count differs: " + name)
            values = numeric_values(match[2],1)
            require(len(values) == patches[name]["count"], "Boundary phi data differs: " + name)
            flow[name] = math.fsum(values)
    denominator = abs(flow["inlet"])
    relative = abs(math.fsum(flow.values()))/denominator if denominator else None
    inlet_outlet = abs(flow["inlet"]+flow["outlet"])/denominator if denominator else None
    return dict(signed_outward_volume_flux_m3_s=flow,
                net_outward_volume_flux_m3_s=math.fsum(flow.values()),
                relative_all_physical_boundary_imbalance=relative,
                relative_inlet_outlet_imbalance=inlet_outlet,
                inlet_inflow_outlet_outflow=flow["inlet"] < 0 < flow["outlet"],
                density_kg_m3=1, mass_flux_kg_s_equals_volume_flux_numerically=True,
                empty_patches_excluded=["frontAndBackPlanes"],
                definition="abs(sum signed outward patch phi)/abs(inlet phi); inlet/outlet-only value also reported")


def residual_steps(log_text, source_name):
    sections = re.split(r"^Time =\s*("+FLOAT+r")\s*s?\s*$",log_text,flags=re.M)
    rows = []
    for i in range(1,len(sections),2):
        time = float(sections[i])
        require(math.isfinite(time), "Non-finite native solver time")
        solves = []
        for variable, initial, final, iterations in re.findall(
                r"Solving for\s+(\w+),\s*Initial residual =\s*([^,\s]+),\s*Final residual =\s*([^,\s]+),\s*No Iterations\s+(\d+)",sections[i+1]):
            if variable in RESIDUAL_VARIABLES or variable == "Uz":
                solves.append((variable,float(initial),float(final),int(iterations)))
        row = dict(time_s=time, source_log=source_name, solve_count=len(solves),
                   all_printed_residuals_finite=all(math.isfinite(v) for s in solves for v in s[1:3]),
                   variables={})
        for variable in RESIDUAL_VARIABLES:
            values=[s for s in solves if s[0]==variable]
            if values and all(math.isfinite(v) for s in values for v in s[1:3]):
                row["variables"][variable]=dict(first_initial=values[0][1],maximum_initial=max(s[1] for s in values),
                    terminal_initial=values[-1][1],terminal_final=values[-1][2],
                    maximum_final=max(s[2] for s in values),linear_solve_count=len(values),
                    total_linear_iterations=sum(s[3] for s in values))
        row["variables_complete"]=set(row["variables"])==set(RESIDUAL_VARIABLES)
        row["terminal_initial_target_passed"]=row["variables_complete"] and all(
            row["variables"][v]["terminal_initial"]<=RESIDUAL_TARGET for v in RESIDUAL_VARIABLES)
        row["terminal_final_target_passed"]=row["variables_complete"] and all(
            row["variables"][v]["terminal_final"]<=RESIDUAL_TARGET for v in RESIDUAL_VARIABLES)
        continuity=re.findall(r"time step continuity errors\s*:\s*sum local =\s*([^,]+),\s*global =\s*([^,]+),\s*cumulative =\s*([^\s]+)",sections[i+1])
        if continuity:
            values=[float(x) for x in continuity[-1]]
            row["terminal_continuity_error"]=values if all(map(math.isfinite,values)) else None
        courants=re.findall(r"Courant Number mean:\s*("+FLOAT+r")\s+max:\s*("+FLOAT+r")",sections[i+1])
        if courants:
            row["maximum_logged_Courant"]=max(float(v[1]) for v in courants)
        rows.append(row)
    return rows



def fatal_solver_error(text):
    # Foundation's startup notice enables trapping; it is not a solver crash.
    notice=r"^\s*sigFpe\s*:\s*Enabling floating point exception trapping \(FOAM_SIGFPE\)\.\s*$"
    errors=re.sub(notice,"",text,flags=re.M|re.I)
    return bool(re.search(r"FOAM FATAL|Floating point exception|Segmentation fault|SIGSEGV|MPI_ABORT",errors,re.I))


def force_data(text):
    rows=[]
    header=None
    for line in text.splitlines():
        words=line.lstrip("# ").split()
        if line.lstrip().startswith("#"):
            if all(name in words for name in ("Time","Cd","Cl")):
                header=words
            continue
        if not words:
            continue
        require(header is not None,"Missing native force coefficient header")
        values=list(map(float,words))
        require(len(values)==len(header) and all(map(math.isfinite,values)),"Invalid native force coefficient sample")
        row=dict(zip(header,values)); rows.append((row["Time"],row["Cd"],row["Cl"]))
    require(rows,"No actual native force samples")
    return rows


def analyse(run, sources):
    parent=run.parent
    setup=json.loads(sources.text(parent/"mesh/setup.json"))
    declared=sources.json("declared-comparison.json")
    status=sources.json("run-status.json")
    preflight=sources.json("provenance/preflight-validation.json")
    profile=sources.json("provenance/native-inlet-profile.json")
    freeze=sources.json("provenance/final-preparation-source-sha256-v002.json")
    sources.text(Path(__file__))
    mesh=parent/"mesh/vortex.msh"; sources.text(mesh)
    sources.local("input-vortex.msh")
    lab=Path(__file__).resolve().parents[3]
    driver=lab/"repo/cfd-learning-lab/scripts/prepare_openfoam_vortex.py"
    wrapper=lab/"local/vortex-shedding-v001/run-openfoam.sh"
    sources.text(driver);sources.text(wrapper)
    cells=int(setup["cells"])
    checks={}
    checks["shared_source_mesh_hash_matches"]=sha(mesh)==sha(run/"input-vortex.msh")==setup["mesh_sha256"]==declared["mesh_source_sha256"]
    checks["frozen_preparation_source_hashes_match"]=sha(driver)==freeze["prepare_openfoam_vortex.py"] and sha(wrapper)==freeze["run-openfoam.sh"]
    checks["preflight_profile_and_quality_passed"]=all(preflight.get(key) is True for key in ("native_face_order_verified","controlDict_matches_preserved_pilot","no_overwrite_guards_verified","mesh_check_passed"))
    checks["initial_native_inlet_field_hash_matches"]=sha(run/"0/U")==profile["applied_U_sha256"]==preflight["active_inlet_field_sha256"]
    control=clean(sources.local("system/controlDict"))
    schemes=clean(sources.local("system/fvSchemes"))
    solution=clean(sources.local("system/fvSolution"))
    physical=clean(sources.local("constant/physicalProperties"))
    transport=clean(sources.local("constant/momentumTransport"))
    decompose=clean(sources.local("system/decomposeParDict"))
    sources.local("0/U");sources.local("0/p")
    pilot_control=clean(sources.local("provenance/initial-system/controlDict"))
    dt=scalar_setting(control,"deltaT");end=scalar_setting(control,"endTime")
    require(end in (2,4),"Only the bounded 2-second pilot or 4-second continuation is supported")
    # Executed new versions are separate sources; original frozen checks remain intact.
    pilot_wrapper=lab/"local/vortex-shedding-v001/run-openfoam-v002.sh"
    sources.text(pilot_wrapper)
    expected_steps=round(end/dt)
    actual_settings=dict(delta_t_s=dt,end_time_s=end,saved_field_interval_s=dt*scalar_setting(control,"writeInterval"),
        viscosity_m2_s=scalar_setting(physical,"nu"),mpi_processes=scalar_setting(decompose,"numberOfSubdomains"),
        pimple_outer_correctors=scalar_setting(flat_block(solution,"PIMPLE"),"nOuterCorrectors"),
        pressure_correctors=scalar_setting(flat_block(solution,"PIMPLE"),"nCorrectors"))
    force=flat_block(control,"forceCoefficients")
    force_refs={name:scalar_setting(force,name) for name in ("rhoInf","magUInf","lRef","Aref")}
    force_refs["CofR_m"]=list(map(float,re.search(r"\bCofR\s*\(([^)]*)\)",force)[1].split()))
    force_refs["sample_interval_s"]=scalar_setting(force,"writeInterval")*dt
    checks["actual_properties_and_transient_settings_match"]=all((near(dt,.005),near(actual_settings["saved_field_interval_s"],.025 if end==2 else .05),near(actual_settings["viscosity_m2_s"],.001),actual_settings["mpi_processes"]==2,actual_settings["pimple_outer_correctors"]==(3 if end==2 else 8),actual_settings["pressure_correctors"]==2,"simulationType laminar;" in transport,"viscosityModel constant;" in physical,"default backward;" in schemes,"Gauss linearUpwind grad(U);" in schemes,"solver incompressibleFluid;" in control,"adjustTimeStep false;" in control,not re.search(r"\bresidualControl\b",solution)))
    checks["force_reference_and_sample_settings_match"]=all((near(force_refs["rhoInf"],1),near(force_refs["magUInf"],1),near(force_refs["lRef"],.1),near(force_refs["Aref"],.1),force_refs["CofR_m"]==[.2,.2,0],near(force_refs["sample_interval_s"],.005),"patches (cylinder);" in force,"rho rhoInf;" in force))
    checks["declared_geometry_properties_match_shared_setup"]=all((cells==declared["intended_cells"],setup["Re_D"]==100,near(setup["diameter_m"],.1),setup["cylinder_centre_m"]==[.2,.2],near(setup["channel_height_m"],.41),near(setup["channel_length_m"],2.2),near(setup["density_kg_m3"],1),near(setup["reference_depth_m"],1),near(setup["kinematic_viscosity_m2_s"],.001),near(setup["inlet_mean_velocity_m_s"],1),near(setup["inlet_maximum_velocity_m_s"],1.5)))
    if end==2:
        checks["pilot_control_matches_preserved_initial_dictionary"]=control==pilot_control
    else:
        continuation="provenance/continuation-to-4s"
        authorization=sources.json(continuation+"/root-authorization.json")
        checks["continuation_has_explicit_root_authorization"]=authorization.get("approved") is True and authorization.get("approved_by")=="root" and authorization.get("action")=="continue_openfoam_vortex_to_4s"
        checks["preserved_pilot_control_retains_2s_cap"]=near(scalar_setting(pilot_control,"endTime"),2) and near(scalar_setting(pilot_control,"deltaT"),.005)
        checks["preserved_pilot_sampling_interval_matches"]=near(scalar_setting(pilot_control,"deltaT")*scalar_setting(pilot_control,"writeInterval"),.025)
        pilot_solution=clean(sources.local(continuation+"/pilot-fvSolution"))
        continued_solution=clean(sources.local(continuation+"/continuation-fvSolution"))
        pilot_pimple=flat_block(pilot_solution,"PIMPLE")
        checks["preserved_pilot_pimple_settings_match"]=scalar_setting(pilot_pimple,"nOuterCorrectors")==3 and scalar_setting(pilot_pimple,"nCorrectors")==2 and not re.search(r"\bresidualControl\b",pilot_solution)
        checks["active_continuation_solution_matches_archived_dictionary"]=solution==continued_solution
        clone_driver=lab/"repo/cfd-learning-lab/scripts/prepare_openfoam_vortex_v003.py"
        clone_wrapper=lab/"local/vortex-shedding-v001/run-openfoam-extended-to4s-v003.sh"
        executed_wrapper=run/"provenance/run-openfoam-continuation-to-4s.sh"
        sources.text(clone_driver);sources.text(clone_wrapper);sources.text(executed_wrapper)
        checks["executed_continuation_wrapper_matches_versioned_source"]=sha(executed_wrapper)==sha(clone_wrapper)
        binding=sources.json(continuation+"/source-sha256.json")
        checks["continuation_source_and_dictionary_hashes_match"]=all((
            binding.get("driver_sha256")==sha(clone_driver),
            binding.get("wrapper_sha256")==sha(clone_wrapper),
            binding.get("pilot_fvSolution_sha256")==sha(run/continuation/"pilot-fvSolution"),
            binding.get("continuation_fvSolution_sha256")==sha(run/continuation/"continuation-fvSolution")))
        checks["continuation_provenance_settings_match"]=all(
            near(float(binding.get(name,-1)),value) for name,value in dict(
                pilot_outer_correctors=3,continuation_outer_correctors=8,pressure_correctors=2,
                delta_t_s=.005,authorized_end_time_s=4,saved_field_interval_s=.05,
                pilot_saved_field_interval_s=.025).items())
    meshlog=sources.local("mesh-check.log")
    quality=dict(cells=int(re.search(r"\bcells:\s+(\d+)",meshlog)[1]),
        minimum_volume_m3=float(re.search(r"Min volume =\s*("+FLOAT+r")",meshlog)[1]),
        maximum_nonorthogonality_degrees=float(re.search(r"Mesh non-orthogonality Max:\s*("+FLOAT+r")",meshlog)[1]),
        maximum_skewness=float(re.search(r"Max skewness =\s*("+FLOAT+r")",meshlog)[1]),
        full_geometry_and_topology="checkMesh -allGeometry -allTopology" in meshlog,
        mesh_ok="Mesh OK." in meshlog)
    checks["full_native_mesh_check_passed"]=quality["mesh_ok"] and quality["full_geometry_and_topology"] and quality["cells"]==cells and quality["minimum_volume_m3"]>0
    poly=run/profile["native_mesh_location"]
    for name in ("points","faces","boundary","owner","neighbour"):
        sources.text(poly/name)
    checks["native_initial_mesh_hash_matches_preflight"]=all(sha(poly/name)==digest for name,digest in profile["native_mesh_sha256"].items())
    patches=native_boundary(sources.text(poly/"boundary"))
    checks["native_boundary_counts_match_shared_mesh"]=all(patches[n]["count"]==v for n,v in dict(cylinder=128,inlet=41,outlet=41,walls=570,frontAndBackPlanes=2*cells).items())
    step_rows=[];log_summaries=[]
    for relative in ["solver.log"]+(["solver-continuation-to-4s.log"] if end==4 else []):
        text=sources.local(relative)
        rows=residual_steps(text,relative);step_rows.extend(rows)
        log_summaries.append(dict(path=relative,has_End=bool(re.search(r"^End\s*$",text,re.M)),steps=len(rows),last_time_s=rows[-1]["time_s"] if rows else None,fatal_error_detected=fatal_solver_error(text)))
    # Restart logs may print an identical saved boundary time; prefer the latest log.
    by_time={row["time_s"]:row for row in step_rows};step_rows=[by_time[t] for t in sorted(by_time) if t>0]
    require(step_rows,"Actual native solver log has no transient steps")
    final=step_rows[-1]
    checks["execution_reached_declared_end_time"]=near(final["time_s"],end) and all(v["has_End"] and not v["fatal_error_detected"] for v in log_summaries)
    checks["wrapper_reports_finished_execution"]=status.get("stage")=="finished_execution" and near(float(status.get("final_time_s",-1)),end)
    checks["expected_fixed_time_steps_present"]=len(step_rows)==expected_steps and all(near(row["time_s"],(i+1)*dt) for i,row in enumerate(step_rows))
    checks["all_printed_residuals_finite"]=all(r["all_printed_residuals_finite"] for r in step_rows)
    checks["all_step_residual_variables_present"]=all(r["variables_complete"] for r in step_rows)
    checks["final_step_terminal_initial_p_U_residuals_le_1e_5"]=final["terminal_initial_target_passed"]
    checks["final_step_terminal_final_p_U_residuals_le_1e_5"]=final["terminal_final_target_passed"]
    residual_record=dict(target=RESIDUAL_TARGET,terminal_initial_definition="Last printed initial linear-system residual for each p/Ux/Uy within the physical step, after final outer correction",residual_norms_equivalent_to_fluent=False,final_step=final,physical_steps=len(step_rows),diagnostic_terminal_initial_target_misses=sum(not r["terminal_initial_target_passed"] for r in step_rows),diagnostic_terminal_final_target_misses=sum(not r["terminal_final_target_passed"] for r in step_rows),diagnostic_missing_variable_steps=sum(not r["variables_complete"] for r in step_rows),diagnostic_nonfinite_steps=sum(not r["all_printed_residuals_finite"] for r in step_rows),earlier_step_threshold_misses_are_final_acceptance_gates=False,logs=log_summaries)
    time_name=format(end,"g")
    fields={}
    for name,components,dimensions,count in (("U",3,[0,1,-1,0,0,0,0],cells),("p",1,[0,2,-2,0,0,0,0],cells),("C",3,[0,1,0,0,0,0,0],cells),("phi",1,[0,3,-1,0,0,0,0],min(p["start"] for p in patches.values()))):
        raw=sources.local(time_name+"/"+name)
        fields[name]=native_field(raw,count,components,dimensions)
        if name=="phi":
            flux=boundary_mass(raw,patches)
    checks["final_native_fields_finite_and_expected_counts"]=all(v["finite"] for v in fields.values())
    checks["two_dimensional_native_solution"]=fields["U"]["maximum_absolute_z"]<=1e-12 and fields["C"]["maximum_absolute_z"]<=1e-10
    checks["inlet_inflow_and_outlet_outflow"]=flux["inlet_inflow_outlet_outflow"]
    checks["relative_inlet_outlet_mass_imbalance_le_1e_4"]=flux["relative_inlet_outlet_imbalance"] is not None and flux["relative_inlet_outlet_imbalance"]<=MASS_TARGET
    checks["relative_all_physical_boundary_mass_imbalance_le_1e_4"]=flux["relative_all_physical_boundary_imbalance"] is not None and flux["relative_all_physical_boundary_imbalance"]<=MASS_TARGET
    summary_name="numerical-summary-pilot.json" if end==2 else "numerical-summary-4s.json"
    receipt=sources.json(summary_name)
    checks["native_phi_agrees_with_final_flux_receipt"]=near(float(receipt["end_time_s"]),end) and all(near(flux["signed_outward_volume_flux_m3_s"][n],float(receipt["signed_patch_volume_flow_m3_s"][n]),absolute=1e-10,relative=1e-9) for n in PHYSICAL_PATCHES)
    sampled_times=[]
    for path in run.iterdir():
        if path.is_dir():
            try:value=float(path.name)
            except ValueError:continue
            if 0<value<=end+1e-8 and all((path/n).is_file() for n in ("U","p","phi")):
                sampled_times.append(value)
    sampled_times.sort()
    expected_snapshot_times=[i*.025 for i in range(1,81)]
    if end==4:
        expected_snapshot_times.extend(2+i*.05 for i in range(1,41))
    checks["saved_native_snapshot_times_complete"]=len(sampled_times)==len(expected_snapshot_times) and all(near(actual,expected) for actual,expected in zip(sampled_times,expected_snapshot_times))
    raw_forces={}
    force_sources=sorted((run/"postProcessing/forceCoefficients").glob("*/forceCoeffs.dat"))
    for path in force_sources:
        for row in force_data(sources.text(path)):
            if 0<=row[0]<=end+1e-8:raw_forces[row[0]]=row
    require(raw_forces,"Missing actual native force coefficient data")
    csvtext=sources.local("force-history.csv")
    force_rows=[]
    for row in csv.DictReader(csvtext.splitlines()):
        values=tuple(float(row[n]) for n in ("time_s","Cd","Cl"))
        require(all(map(math.isfinite,values)),"Non-finite force CSV sample")
        force_rows.append(values)
    checks["force_csv_matches_actual_native_coefficients"]=len(force_rows)==len(raw_forces) and all(row[0] in raw_forces and all(near(a,b,absolute=1e-12,relative=1e-12) for a,b in zip(row,raw_forces[row[0]])) for row in force_rows)
    positive_forces=[r for r in force_rows if r[0]>0]
    checks["native_force_sample_times_complete"]=len(positive_forces)==expected_steps and all(near(row[0],(i+1)*dt) for i,row in enumerate(positive_forces))
    record=dict(solver="OpenFOAM Foundation 14 incompressibleFluid",analysis_completed=True,actual_solver_evidence_present=True,mesh_cells=cells,actual_settings=actual_settings,force_references=force_refs,mesh_quality=quality,initial_native_profile=dict(n_faces=profile["n_faces"],volume_flow_m3_s=profile["discrete_volume_flow_m3_s"],continuum_flow_m3_s=.41,rescaled=False,midpoint_quadrature_error_percent=profile["discrete_flux_error_percent"]),residuals=residual_record,boundary_mass_balance=flux,final_native_field_statistics=fields,final_flux_receipt=receipt,saved_snapshot_times_s=sampled_times,expected_saved_snapshot_times_s=expected_snapshot_times,force_history=dict(samples=len(force_rows),first_time_s=force_rows[0][0] if force_rows else None,last_time_s=force_rows[-1][0] if force_rows else None,raw_sources=[p.relative_to(run).as_posix() for p in force_sources],force_coefficients_required_to_be_constant=False,periodicity_required_for_pilot=False),acceptance_thresholds=dict(final_step_terminal_initial_residual=RESIDUAL_TARGET,final_step_terminal_final_residual=RESIDUAL_TARGET,relative_boundary_mass_imbalance=MASS_TARGET),check_results=checks,failed_checks=[name for name,value in checks.items() if not value],all_declared_checks_passed=all(checks.values()),physical_validation=False,benchmark_validation=False,mesh_independence_established=False,periodic_shedding_established=False,interpretation="Bounded numerical checks on actual startup transient evidence; force oscillations are expected, and this short pilot is not a periodic or benchmark validation")
    return record,step_rows,force_rows


def write_json(path,data):
    with Path(path).open("x",encoding="utf-8",newline="\n") as stream:
        json.dump(data,stream,indent=2,allow_nan=False);stream.write("\n")


def export(out,record,manifest,steps,forces):
    out.mkdir(parents=True)
    write_json(out/"checks.json",record)
    write_json(out/"source-manifest.json",manifest)
    if record.get("actual_solver_evidence_present"):
        with (out/"residual-history.csv").open("x",encoding="ascii",newline="") as stream:
            columns=["time_s","source_log","solve_count","all_printed_residuals_finite","variables_complete","terminal_initial_target_passed","terminal_final_target_passed"]
            columns += [v+"_"+stat for v in RESIDUAL_VARIABLES for stat in ("first_initial","maximum_initial","terminal_initial","terminal_final","maximum_final","linear_solve_count")]
            writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
            for row in steps:
                flat={key:row[key] for key in columns if key in row}
                for variable,data in row["variables"].items():
                    for key,value in data.items():
                        if variable+"_"+key in columns:flat[variable+"_"+key]=value
                writer.writerow(flat)
        with (out/"force-history.csv").open("x",encoding="ascii",newline="") as stream:
            writer=csv.writer(stream);writer.writerow(("time_s","Cd","Cl"));writer.writerows(forces)
        # Small archive contains derived checks and provenance, not native fields/logs.
        with zipfile.ZipFile(out/"evidence.zip","x",compression=zipfile.ZIP_STORED) as archive:
            for name in ("checks.json","source-manifest.json","residual-history.csv","force-history.csv"):
                archive.write(out/name,arcname=name)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    run,out=args.run.resolve(),args.out.resolve()
    require(run.is_dir(),"Native OpenFOAM run folder does not exist")
    require(not out.exists(),"Choose a fresh output directory; nothing is overwritten")
    require(not out.is_relative_to(run),"Keep analysis output outside the native case")
    sources=Sources(run);steps=[];forces=[]
    try:
        record,steps,forces=analyse(run,sources)
    except (ValueError,KeyError,TypeError,IndexError,OSError,AttributeError,OverflowError) as error:
        record=dict(analysis_completed=False,actual_solver_evidence_present=False,all_declared_checks_passed=False,check_results=dict(actual_evidence_complete_and_parseable=False),failed_checks=["actual_evidence_complete_and_parseable"],error=str(error),physical_validation=False,benchmark_validation=False,solver_started=False)
    manifest=sources.finish()
    record["inputs_unchanged_during_analysis"]=manifest["inputs_unchanged_during_analysis"]
    record["check_results"]["inputs_unchanged_during_analysis"]=manifest["inputs_unchanged_during_analysis"]
    record["all_declared_checks_passed"]=all(record["check_results"].values())
    record["failed_checks"]=[name for name,value in record["check_results"].items() if not value]
    export(out,record,manifest,steps,forces)
    print(json.dumps(dict(output_directory=str(out),analysis_completed=record["analysis_completed"],all_declared_checks_passed=record["all_declared_checks_passed"],failed_checks=record["failed_checks"]),indent=2))
    return 0 if record["all_declared_checks_passed"] else 2


if __name__=="__main__":
    raise SystemExit(main())
