"""Audit the resolved-plenum Helmholtz case and plot native probe histories."""
import argparse
import hashlib
import json
import re
import runpy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CELLS = 25344
DT = 1e-4
STATES = ("0.01", "0.02", "0.03", "0.04", "0.05")
PROBES = [(-.045, 0, 0), (-.045, .02, 0), (-.01, 0, 0),
          (.0125, 0, 0), (.0125, .02, 0)]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, reader, out):
    read = runpy.run_path(str(reader))["foam_values"]
    sources = ["input-sha256.txt", "check-mesh.log", "solver.log",
               "postProcessing/probes/0/p"]
    for line in (case/sources[0]).read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        assert sha(case/name) == digest, name
    mesh = (case/sources[1]).read_text(encoding="ascii")
    log = (case/sources[2]).read_text(encoding="ascii")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*25344\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd")
    expected = DT*np.arange(501)
    step_times = np.array(re.findall(r"^Time = ([0-9.eE+-]+)s$", log, re.M), float)
    assert len(step_times) == 500 and np.allclose(step_times, expected[1:], rtol=0, atol=1e-12)

    def field(name, kind="scalar"):
        sources.append(name)
        values = read((case/name).read_text(encoding="ascii"), "internalField", CELLS, kind)
        assert np.all(np.isfinite(values)), name
        return values

    saved = []
    for index, time in enumerate(STATES, 1):
        name = f"{time}/uniform/time"
        sources.append(name)
        clock = (case/name).read_text(encoding="ascii")
        assert int(re.search(r"^index\s+(\d+);", clock, re.M)[1]) == index*100
        assert abs(float(re.search(r"^value\s+([^;]+);", clock, re.M)[1])-float(time)) < 1e-12
        assert float(re.search(r"^deltaT\s+([^;]+);", clock, re.M)[1]) == DT
        u, p, temp = field(f"{time}/U", "vector"), field(f"{time}/p"), field(f"{time}/T")
        assert u.shape == (CELLS, 3) and p.shape == temp.shape == (CELLS,)
        assert np.all(p > 0) and np.all(temp > 0)
        saved.append({"time_s": float(time), "min_pressure_pa": float(p.min()),
                      "max_pressure_pa": float(p.max()), "min_temperature_k": float(temp.min()),
                      "max_temperature_k": float(temp.max()),
                      "max_speed_m_s": float(np.linalg.norm(u, axis=1).max())})
    probe_file = case/sources[3]
    header = probe_file.read_text(encoding="ascii")
    points = re.findall(r"^# Probe (\d+) \(([^)]+)\)\s*$", header, re.M)
    assert [int(i) for i, _ in points] == list(range(5))
    assert np.array_equal(np.array([list(map(float, p.split())) for _, p in points]), PROBES)
    history = np.loadtxt(probe_file)
    assert history.shape == (501, 6) and np.all(np.isfinite(history))
    assert np.all(history[:, 1:] > 0) and np.all(np.diff(history[:, 0]) > 0)
    assert np.allclose(history[:, 0], expected, rtol=0, atol=1e-12)
    # Native sample clocks already coincide with the uniform analysis grid.
    time = history[100:, 0]
    neck = history[100:, 3]
    centered = neck-neck.mean()
    peaks = np.flatnonzero((centered[1:-1] > centered[:-2]) &
                           (centered[1:-1] > centered[2:]) & (centered[1:-1] > 0))+1
    assert np.ptp(neck) > 1 and len(peaks) >= 3 and np.all(np.diff(peaks) >= 3)
    window = np.hanning(len(time))
    frequencies = np.fft.rfftfreq(len(time), DT)
    amplitude = 2*np.abs(np.fft.rfft(centered*window))/window.sum()
    dominant = 1+int(np.argmax(amplitude[1:]))
    spacings = np.diff(time[peaks])

    centres, volumes = field("0/C", "vector"), field("0/Vc")
    assert np.all(volumes > 0)
    x, y, z = centres.T
    transverse = np.maximum(abs(y), abs(z))
    # Plenum box and its outer annulus, excluding the rounded entrance cells.
    plenum = (x < -.035) | ((x < -.020) & (transverse > .0125))
    throat = (x > -.020) & (x < 0) & (transverse < .005)
    volume = float(volumes[plenum].sum())
    length = .020
    area = float(volumes[throat].sum()/length)
    assert np.isclose(area, 1e-4, rtol=1e-6) and volume > 0
    gas_r = 8314.462618/28.9
    gamma = 1+gas_r/712
    initial_temperature = field("0/T")
    assert np.all(initial_temperature == initial_temperature[0]) and initial_temperature[0] > 0
    sound_speed = float(np.sqrt(gamma*gas_r*initial_temperature[0]))
    corrected_length = float(length+1.7*np.sqrt(area/np.pi))
    lumped = {name: sound_speed/(2*np.pi)*np.sqrt(area/(volume*value))
              for name, value in (("straight_length_hz", length),
                                  ("illustrative_end_correction_hz", corrected_length))}

    result = {"cells": CELLS, "steps": 500, "time_step_s": DT, "saved_fields": saved,
              "probe_locations_m": PROBES, "probe_samples": len(history),
              "neck_window_s": [float(time[0]), float(time[-1])],
              "neck_swing_pa": float(np.ptp(neck)),
              "positive_peak_times_s": time[peaks].tolist(), "peak_spacings_s": spacings.tolist(),
              "median_peak_spacing_s": float(np.median(spacings)),
              "fft_dominant_bin_hz": float(frequencies[dominant]),
              "fft_bin_spacing_hz": float(frequencies[1]), "interpolated": False,
              "plenum_volume_m3": volume, "straight_throat_area_m2": area,
              "straight_throat_length_m": length, "illustrative_corrected_length_m": corrected_length,
              "initial_sound_speed_m_s": sound_speed, "lumped_frequency_estimates": lumped,
              "native_sha256": {name: sha(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    np.savetxt(out/"probe-pressure.csv", history, delimiter=",", comments="",
               header="time_s,plenum_centre_pa,plenum_off_axis_pa,neck_pa,outlet_centre_pa,outlet_off_axis_pa",
               fmt="%.10g")
    np.savetxt(out/"neck-spectrum.csv", np.column_stack([frequencies, amplitude]),
               delimiter=",", comments="", header="frequency_hz,window_corrected_amplitude_pa")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.1), layout="constrained")
    for col, label in ((1, "Plenum centre"), (3, "Neck"), (4, "Outlet centre")):
        axes[0].plot(history[:, 0]*1000, history[:, col]-history[0, col], label=label)
    axes[0].set(xlabel="Time (ms)", ylabel="Pressure change (Pa)", title="Pressure relative to start")
    axes[0].legend(frameon=False)
    axes[1].plot(time*1000, centered, color="#bf5831")
    axes[1].scatter(time[peaks]*1000, centered[peaks], s=15, color="#173f5f", zorder=3)
    axes[1].set(xlabel="Time (ms)", ylabel="Pressure minus window mean (Pa)", title="Neck ring-down")
    show = (frequencies > 0) & (frequencies <= 1000)
    axes[2].plot(frequencies[show], amplitude[show], "o-", ms=3, color="#173f5f")
    axes[2].axvline(frequencies[dominant], color="#bf5831", alpha=.75, lw=1)
    axes[2].set(xlabel="Frequency (Hz)", ylabel="Amplitude (Pa)", title="Neck spectrum, 10–50 ms")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle("Pressure oscillations in a Helmholtz resonator", fontweight="bold")
    for extension in ("png", "pdf"):
        fig.savefig(out/f"helmholtz-resonance.{extension}", dpi=180,
                    metadata={"CreationDate": None, "ModDate": None} if extension == "pdf" else None)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("case", "reader", "output"):
        parser.add_argument(argument, type=Path)
    args = parser.parse_args()
    audit(args.case, args.reader, args.output)