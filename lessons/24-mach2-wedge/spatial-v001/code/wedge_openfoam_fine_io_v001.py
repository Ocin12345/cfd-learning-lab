from pathlib import Path
import re
from math import tan
import numpy as np
from wedge_common_fine_v001 import NX,NY,CORNER,DEPTH,H,L,P,R,T,THETA,U
from prepare_wedge_openfoam_fine_v001 import END_TIME,FRAME_DT,MAX_CO
from plot_wedge_fluent_fine_v001 import fields,measure
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

def read_openfoam(run, geo):
    from analyse_nozzle_openfoam_v002 import FIELDS, field, patch_values

    assert (run / "run-status.txt").read_text().strip() == "complete"
    mesh_log = (run / "mesh-check.log").read_text()
    assert "Mesh OK." in mesh_log and re.search(r"\bcells:\s+89600\b", mesh_log)
    log = (run / "solver.log").read_text()
    assert re.search(r"^End\s*$", log, re.M) and "Finalising parallel run" in log
    assert "FOAM FATAL" not in log
    step_times = np.array([float(t) for t in re.findall(r"^Time = (" + NUMBER + r")s?\s*$", log, re.M)])
    assert len(step_times) > 0 and np.all(np.diff(step_times) > 0)
    assert abs(step_times[-1] - END_TIME) < 1e-10
    co = np.array([float(c) for c in re.findall(
        r"^Courant Number mean: " + NUMBER + r" max: (" + NUMBER + r")", log, re.M)])
    assert len(co) == len(step_times) + 1 and co.max() <= MAX_CO * 1.02  # Includes startup report.

    count = NX * NY
    native_centres = field(run / "0.00025/C", count, *FIELDS["C"])
    x, y, z = native_centres.T
    lower = np.maximum(x - CORNER, 0) * tan(THETA)
    i = np.floor(x * NX / L).astype(int)
    j = np.floor((y - lower) * NY / (H - lower)).astype(int)
    order = j * NX + i
    assert ((i >= 0) & (i < NX) & (j >= 0) & (j < NY)).all()
    assert len(np.unique(order)) == count
    coord_error = np.linalg.norm(native_centres[:, :2] - geo["centres"][order], axis=1).max()
    assert coord_error < 2e-7 and np.max(abs(z - DEPTH / 2)) < 1e-9
    to_common = np.argsort(order)

    frame_times = np.arange(round(END_TIME / FRAME_DT) + 1) * FRAME_DT
    assert "internalField uniform 101325;" in (run / "0/p").read_text()
    assert "internalField uniform 300;" in (run / "0/T").read_text()
    assert f"internalField uniform ({U:.16g} 0 0);" in (run / "0/U").read_text()
    zeros = np.zeros((NY, NX))
    frames = [fields(dict(p=zeros + P, T=zeros + T, rho=zeros + P / (R * T),
                          u=zeros + U, v=zeros.copy()))]
    mass_rows = []
    for t in frame_times[1:]:
        directory = run / f"{t:.9g}"
        scalars = {name: field(directory / name, count, *FIELDS[name])[to_common].reshape(NY, NX)
                   for name in ("p", "rho", "T")}
        vector = field(directory / "U", count, *FIELDS["U"])[to_common].reshape(NY, NX, 3)
        assert abs(vector[..., 2]).max() < 1e-9
        frames.append(fields(dict(**scalars, u=vector[..., 0], v=vector[..., 1])))
        phi = directory / "phi"
        assert re.search(r"\bdimensions\s*\[\s*1\s+0\s+-1\s+0\s+0\s+0\s+0\s*\];", phi.read_text())
        fluxes = {name: float(patch_values(phi, name, faces).sum()) / DEPTH for name, faces
                  in (("inlet", NY), ("outlet", NY), ("top", NX), ("ramp", NX))}
        assert fluxes["inlet"] < 0 < fluxes["outlet"]
        mass_rows.append(dict(time_s=float(t), **fluxes,
                              relative_imbalance=abs(sum(fluxes.values()) / fluxes["inlet"])))
    native = dict(final_time_s=float(frame_times[-1]), comparison_time_s=float(frame_times[-5]),
                  time_steps=len(step_times), maximum_reported_courant=float(co.max()),
                  mesh_check="Mesh OK.", native_centroid_alignment_error_m=float(coord_error),
                  final_native_outward_mass_flux_kg_s_per_m_depth=mass_rows[-1],
                  mass_imbalance=mass_rows[-1]["relative_imbalance"],
                  last_five_frames_max_mass_imbalance=max(row["relative_imbalance"] for row in mass_rows[-5:]),
                  initial_frame="Prescribed uniform field")
    return frame_times, frames, mass_rows, native
