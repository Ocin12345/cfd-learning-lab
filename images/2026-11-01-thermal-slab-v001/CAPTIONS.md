# Actual transient solid-slab temperatures

Generic constant-property model solid: rho1000kg/m3, c1000J/(kgK), k10W/(mK), alpha1e-5m2/s. Rectangle.1m by.02m, common100x20quads. Initial interior300K; hot left jumps to400K, cold right300K, insulated top/bottom. Both codes use matched first-order implicit Euler,dt2s,100steps to200s. No fluid motion, tracers or synthetic physical times.

All101native temperature states are independently reread and matched to auditedNPZ,202000cell values per solver. Native cell connectivity maps by complete vertex sets. Fluent physical time, timestep index and dt are independently parsed from saved native Data Variables; OpenFOAM uses actual native time folders and solved energy-log times. GIF51frames select exactly0,4,...200s,50ms/frame,2.55s playback. Temperature colour limits stay300..400K and the rectangle retains its physical5:1aspect.

Dotted reference profiles are exact cell averages of the continuous1D Dirichlet PDE, not another numerical solver result. The initial interior is explicitly300K;80vs160series terms agree within1e-9K for allsavedtimes. Reference/data values are not clipped. Errors at every actual saved time, including early coarse-implicit-Euler boundary-jump errors, remain in the CSV and separate error figure. Only final200s has accuracy targets:maxabs0.5K, volume-weightedRMS0.25K. No all-time accuracy threshold is imposed.

E=sum(rho*c*V*(T-300K)); nativeQ is signed positive into the solid. Global backward-Euler balance uses(E_n-E_prev)/2s minus summed nativeQ_n. Paired heat-rate/residual figures use only solved times2..200s. Fluent initialQ is explicitly unsolved; actual OpenFOAM time0boundary powers remain in its native audit and are omitted from paired plots. No time0Euler balance is invented. Energy/power plots divide by native depth1m(Fluent)or.01m(OpenFOAM); rawJ/W remain in metrics. Insulation is reported separately and is not sole conservation proof. No all-cell thermal flux conservation is claimed.

This is numerical verification of an ideal material/PDE model, not experimental physical validation or timestep/mesh independence. The user directs the project; Codex selected/prepared the model and workflow; solvers computed fields. No user reflection or independent user execution is invented.

ANSYS Fluent: declared numerical checks passed
OpenFOAM: declared numerical checks passed
