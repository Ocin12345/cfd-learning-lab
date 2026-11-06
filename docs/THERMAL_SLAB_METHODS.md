# Solid-slab transient conduction methods

[Actual transient visuals](../lessons/13-solid-conduction/README.md).
The selected Fluent v004 and OpenFOAM v002 runs contain complete matched histories: 100 implicit Euler steps of 2 s on the same 2,000-cell geometry. Physical properties, time limits and numerical targets remained unchanged through the documented compatibility corrections. Earlier attempts are preserved. The selected figure v002 improves label readability.

The rectangle x[0,0.1], y[0,0.02] m uses 100 by 20 uniform 0.001 m quads.
All four physical patches are walls: hot/cold have 20 faces each and insulated
bottom/top have 100 each. Fluent solves a planar solid of depth 1 m (volume
0.002 m³); OpenFOAM extrudes depth 0.01 m (volume 0.00002 m³), with empty end
planes. The imported shared zone is converted to solid before energy advance.
The generic constant-property model uses rho=1000 kg/m³, c=1000 J/(kg K),
k=10 W/(m K), alpha=1e-5 m²/s. OpenFOAM constSolidThermo defines Cp=Cv here.
Initial interior T=300 K; hot left jumps to 400 K at time zero; cold right stays
300 K; long edges impose zero heat flux. No motion, sources or radiation.

Both codes use fixed first-order implicit Euler, dt=2 s, exactly 100 steps to
200 s. There are 101 saved native temperature states, including the initial
unsolved state. Fluent caps 20 inner iterations per step, using its own scaled
energy-residual target 1e-10; OpenFOAM solves energy once per step with PCG and
checks the printed final linear residual against 1e-10. These norms differ.
A new-time OpenFOAM initial residual is recorded separately and is not this gate.

| Actual quantity | Fluent | OpenFOAM |
|---|---|---|
| Saved states / solved steps | 101 / 100 | 101 / 100 |
| Maximum own thermal residual | 9.822242141807367e-11 | 9.9639575948522e-13 |
| Maximum global Euler energy error | 8.565010678625376e-08 | 6.386218518425234e-12 |
| Final max analytical error (K) | 0.17406020008468204 | 0.1740601233519783 |
| Final weighted RMS analytical error (K) | 0.12297640164163436 | 0.1229763702719236 |
| All-time maximum analytical error (K) | 12.527370562631688 | 12.527370478194882 |
| Time of all-time maximum error (s) | 2 | 2 |
| Final raw native energy above 300 K (J) | 88518.90493197009 | 885.1890498508678 |

Temperature bounds [300,400] K and column y-spread <=1e-6 K are checked over all
states. Each actual solved step checks global Euler energy conservation:
E=sum(rho*c*V*(T-300 K)); R=(E_n-E_previous)/2 s - sum(Q_into_solid,n).
The relative error is abs(R)/max(sum(abs(Q_each_wall)),abs(storage_rate),Qref),
with Qref=200 W for Fluent and 2 W for OpenFOAM; target 1e-5 at every step.
Native powers are positive into the solid. Fluent's SV_HEAT_FLUX already contains
integrated face W and is crosschecked against native wall reports without an
extra area multiplication. OpenFOAM wallHeatFlux W/m² is integrated once using
actual native face areas and crosschecked against its native table Q in W.
Insulated-wall powers are reported separately. The initial unsolved state has
no invented heat-rate or residual check. This establishes a global balance,
not unobserved local all-cell thermal-flux conservation. Energy/power figures
normalize by each native depth; raw J/W values remain in the metrics.

The [frozen analytical source](../scripts/thermal_slab_reference_v001.py) gives
exact cell averages of the continuous 1D Dirichlet solution; 80 and 160 terms
are crosschecked to 1e-9 K. At t=0 the interior is explicitly 300 K rather than a
truncated Fourier series. FINAL t=200 s targets are max absolute error <=0.5 K
and volume-weighted RMS <=0.25 K. Errors at all saved times remain reported;
there is no all-time accuracy threshold. The initial boundary jump is incompatible
at the boundary and poorly resolved by the early 2 s Euler steps. No clipping,
automatic time extension or target relaxation is used.

[Fluent audit](../results/thermal-slab-v001/fluent/fluent-checks.json)
and [OpenFOAM audit](../results/thermal-slab-v001/openfoam/checks.json)
retain actual booleans unchanged. Structural/parser/source errors and incomplete
native histories block selection. Genuine complete-history numerical misses
remain false. [Selected OpenFOAM dictionaries](../cases/openfoam-thermal-slab-selected-v001/README.md)
and small native settings/metrics are byte-identical copies. All 202,000 native
cell temperatures per solver are reread for figures. The 51-frame GIF selects
actual native times 0,4,...,200 s with 50 ms per frame and fixed 300–400 K colour
limits; no time interpolation or fixed-field tracer animation is used.

Fresh versioned E: runs use existing D: Fluent, E: Python, hidden two-process
Fluent/two-rank OpenFOAM, sequential sessions and live resource/licence checks.
Native HDF, meshes, NPZ, full temperature tables, full traces, dependencies and
private launch state remain external. Prior projects and rejected source/API
attempts remain preserved. No experimental validation, mesh/time independence,
research novelty, invented reflection, Git commit, publication or submission is claimed.
