# Actual steady square-obstacle wake

Finite 2D confined domain x=[0,2] m, y=[-0.4,0.4] m minus the solid square x=[0.45,0.55] m, y=[-0.05,0.05] m. Body-centre distances are 5 D from inlet and 15 D to outlet. Height 8 D gives 12.5% blockage. Top and bottom are symmetry/slip; the body is stationary no-slip. Inlet U=0.1 m/s, rho=1 kg/m3, mu=0.0005 Pa s, D=0.1 m, Re=20. Native fluid fields map through complete vertex sets to 15,900 source quads; the actual body hole remains unfilled.

Pressure retains the common outlet gauge reference 0 Pa. No mean subtraction or invented void values. Both codes are steady; RK4 particles only follow the fixed actual saved velocity. All four interpolation centres must be native fluid cells and every stage/segment is rejected if it touches/crosses the body. Inlet colours identify seed origin, not concentration. No forced path closure or transient vortex-shedding claim.

Force histories contain actual native iterations 1..2000. Stability/symmetry windows are exactly 1801..2000. The separately declared paired Cd target is relative final-window mean difference <=1%, normalized by abs(Fluent mean Cd); a zero mean is undefined and false. Native numerical flags and the paired flag remain separate and every miss is shown. No published Cd, experimental validation, mesh independence or unbounded-flow accuracy claim.

Tracer travel clock 40.0 s, injection period 6.0 s, replay 48 frames at 50 ms. These parameters must be selected by an actual native-field timing probe.

The user directs the project; Codex prepares and operates code; Fluent/OpenFOAM compute the fields. Independent user execution/reflection is not invented.

ANSYS Fluent: declared numerical checks passed
OpenFOAM: declared numerical checks passed
Paired Cd: PASS | actual final-window mean difference 0.15405% | target <= 1%
