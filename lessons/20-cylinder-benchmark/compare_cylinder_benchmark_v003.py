"""Compare saved cylinder forces with published periodic extrema; no CFD launch.

Run from the lesson directory with Python, NumPy and Matplotlib. Pass a fresh
output directory for derived files. The two source histories are read only.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

START, END, TOL = 2.05, 4.0, 1e-8
D, U = 0.1, 1.0
REF = {"Cd": [3.168, 3.230], "Cl": [-1.023, 0.987], "St": 0.300}
REF_URL = "https://riemann.mathematik.uni-dortmund.de/papers/HuebnerMahmoodTurek2012.pdf"
RUNS = {"fluent": ("ANSYS Fluent", "#1765a0", 0.025),
        "openfoam": ("OpenFOAM", "#ae461e", 0.005)}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sample(path, spacing):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        raw = np.array([[float(r[key]) for key in ("time_s", "Cd", "Cl")]
                        for r in csv.DictReader(stream)])
    assert np.isfinite(raw).all() and np.all(np.diff(raw[:, 0]) > 0)
    a = raw[(raw[:, 0] >= START - TOL) & (raw[:, 0] <= END + TOL)]
    expected = START + spacing * np.arange(round((END - START) / spacing) + 1)
    assert a.shape == (len(expected), 3)
    np.testing.assert_allclose(a[:, 0], expected, rtol=0, atol=TOL)
    t, cd, cl = a.T
    j = np.where((cl[:-1] < 0) & (cl[1:] >= 0))[0]
    crossings = t[j] - cl[j] * (t[j + 1] - t[j]) / (cl[j + 1] - cl[j])
    periods = np.diff(crossings)
    assert len(periods) >= 4 and (periods > 0).all()
    return a, {
        "actual_rows": len(a), "native_sample_spacing_s": spacing,
        "sampled_extrema": {"Cd": [float(cd.min()), float(cd.max())],
                            "Cl": [float(cl.min()), float(cl.max())]},
        "rising_zero_crossings_s": crossings.tolist(),
        "cycle_periods_s": periods.tolist(),
        "St_per_cycle": (D / (U * periods)).tolist(),
        "St_from_mean_period": float(D / (U * periods.mean())),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    paths = {k: here / (k + "-force-history.csv") for k in RUNS}
    digests = {k: sha256(p) for k, p in paths.items()}
    data = {k: sample(paths[k], meta[2]) for k, meta in RUNS.items()}
    args.out.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none"})
    fig, ax = plt.subplots(3, 1, figsize=(9.4, 7.7), sharex=True,
                           gridspec_kw={"height_ratios": [1, 1.25, .78]},
                           layout="constrained")
    fig.suptitle("Cylinder wake at Re = 100: startup versus a periodic reference",
                 fontsize=15, fontweight="bold")
    for panel, force, symbol in zip(ax[:2], ("Cd", "Cl"), (r"C_D", r"C_L")):
        panel.axhspan(*REF[force], color="#536579", alpha=.14, zorder=0,
                      label="Published periodic extrema")
        panel.set_ylabel(fr"${symbol}$")
    ax[2].axhline(REF["St"], color="#536579", linestyle="--", lw=1.2,
                  label="Published periodic value")
    for k, (name, color, _) in RUNS.items():
        a, stats = data[k]
        for j in (0, 1):
            ax[j].plot(a[:, 0], a[:, j + 1], color=color, lw=1.6, label=name)
        q = np.array(stats["rising_zero_crossings_s"])
        ax[2].plot((q[1:] + q[:-1]) / 2, stats["St_per_cycle"], color=color,
                   marker="o" if k == "fluent" else "s", ms=4, lw=1.0)
    ax[0].legend(loc="lower right", frameon=False, ncols=3, fontsize=8.8)
    ax[2].legend(loc="upper right", frameon=False, fontsize=8.8)
    ax[2].set(ylabel=r"$St = f D/U$", xlabel="Physical time (s); cycle midpoint in lower panel",
              xlim=(2.02, 4.03))
    for panel in ax:
        panel.grid(True, alpha=.20)
    for suffix in ("png", "svg"):
        fig.savefig(args.out / ("benchmark-force-comparison." + suffix), dpi=220,
                    facecolor="white")
    plt.close(fig)
    assert digests == {k: sha256(p) for k, p in paths.items()}
    summary = {
        "native_window_s": [START, END], "D_m": D, "mean_inlet_U_m_s": U,
        "method": "Rising Cl zero crossings estimated linearly between adjacent saved force samples; no phase alignment.",
        "sources": {k: {"file": p.name, "sha256": digests[k]} for k, p in paths.items()},
        "actual": {k: data[k][1] for k in data},
        "reference": {"values": REF, "source": REF_URL, "locator": "Table 8 reference values"},
        "reference_extrema_are_acceptance_bounds": False,
        "matched_periodic_benchmark_established": False,
        "physical_validation": False,
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)
                                          + "\n", encoding="utf-8")
    for name, stats in summary["actual"].items():
        print(name, stats["actual_rows"], stats["St_from_mean_period"],
              stats["sampled_extrema"])


if __name__ == "__main__":
    main()

