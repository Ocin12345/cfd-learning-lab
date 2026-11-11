# Annular Couette spin-up: methods

Motion spreads outward from the spinning inner wall; both solvers settle to very similar circular shear profiles.

Both solvers computed 400 momentum steps of 0.05 s to 20 s, retaining 401 native states each. Each continued its own accepted five-step pilot for 395 further steps without reinitialization or pilot replay. Own primary checks: fluent PASS; openfoam PASS. All numerical misses remain in the compact result record.

## Geometry and assumptions

A planar annulus has inner radius 0.05 m and outer radius 0.10 m. The shared mesh has 64 angular by 24 radial quadrilaterals (1,536 cells), a true hole and wrapped angular connectivity. Stationary polygonal wall segments approximate circles: polygon area 0.02352411367909455 m² versus circular area 0.023561944901923454 m² (about -0.16056%). Native quality, positive volumes and source/native vertex-cell-face mappings are audited. Fluent uses unit depth; OpenFOAM uses a 0.01 m empty extrusion, with per-depth quantities compared.

The generic constant-density laminar fluid has rho = 1 kg/m³, mu = 0.000125 Pa s and nu = 0.000125 m²/s. Cells start at rest. The inner wall rotates counterclockwise at +1 rad/s about +z, surface speed 0.05 m/s; the outer wall is stationary. Both are no-slip on a fixed mesh in an inertial frame. The gap Reynolds number is 20. Live native units and projected tangential Omega cross r face velocities are checked. Energy, turbulence, gravity, MRF and moving meshes are absent.

## Numerical methods and checks

Both use first-order backward Euler and first-order upwind momentum. Fluent uses transient SIMPLE with at most 100 inner iterations per step; its own terminal continuity/x/y scaled residual target is 1e-6. OpenFOAM uses PIMPLE with the separately amended budget of 30 outer correctors, two pressure correctors and no nonorthogonal correctors; final-outer initial p/U targets remain 1e-6 and final linear targets remain 1e-8. Actual p/Ux/Uy and Uz only if observed are retained. These solver residual norms differ and are plotted separately.

The original five-step, 10-outer OpenFOAM pilot passed 18/19 checks but missed final-outer pressure initial residuals. It is preserved as a numerical attempt. The explicit amendment changes only the outer budget to 30, retaining all physical settings and thresholds. The fresh amended pilot uses only copied, verified native time-zero fields, never the insufficiently converged step-five result. Source hashes and full original audits remain external; compact failed names and evidence hashes are saved here. Setup/API unit spelling and native boundary parser corrections are separately preserved; zero-step preparation failures are not computed physical results.

The OpenFOAM CFD continuation completed all 395 further steps to 20 s, and the independent native audit retains all 401 authoritative states and 800 solved rank clocks. Its raw wrapper nevertheless remains failed (exit 1) at the post-run copied-input preservation check. Exactly two derived restart vorticity copies changed, by a measured maximum 3.9968028886505635e-13 s⁻¹. Both the original strict all-inputs-unchanged diagnostic and wrapper postcheck remain **false** (245 of 247 copies unchanged). The original accepted pilot and primary U, p, face flux, time, mesh and physics inputs remained unchanged. Original pilot states 0 through 5 remain authoritative. The native-completion receipt, exact old/current hashes, derived-field differences and unchanged raw failed status are source-bound externally; compact diagnostics retain both false outcomes. No CFD was replayed and no native output or status was edited to repair the record. A passed primary numerical audit does not mean the wrapper postcheck passed.

Pressure references lie at source cell 0's actual fluid centroid (0.05092584598549641, 0.00250182638511251, 0) m, with OpenFOAM's unique owning rank proved by native addressing. Native Fluent pressure is Pa; OpenFOAM kinematic pressure is converted with rho. Comparisons subtract each actual field's area/volume-weighted mean, retaining raw gauges separately.

This closed domain has no inlet-relative mass metric. Actual native all-face fluxes give every-cell divergence, scaled by gap/surface-speed = 1 s: maximum target 1e-5, volume-weighted RMS target 1e-6. Each wall face flux and total signed boundary flux have absolute targets 1e-12 in each code's native kg/s or m³/s units. Zero net boundary flux alone does not establish local continuity. Positive angular momentum per depth is tracked for solved states. Torque and an angular-momentum balance are not claimed because native torque integration/sign/units have not been independently proved.

Numerically steady wording requires both the actual step 360-to-400 area-weighted velocity L2 change with final-step denominator and relative angular momentum change to be at most 0.001. Otherwise the last field remains a saved transient. The circular steady limit v_theta = -r/3 + 1/(300r) is descriptive and does not describe startup. There is no analytic or cross-code accuracy gate, mesh/time independence, experimental validation or 3D Taylor-vortex stability claim.

## Figures and reproduction

The GIF uses 101 actual native steps 0, 4,..., 400 (0 to 20 s), at 60 ms playback per frame, with a fixed common speed range 0 to 0.05 m/s and equal xy aspect. The arrow annotates prescribed wall rotation on fixed geometry. No temporal interpolation or particles are added; settled or quantized repeats remain. OpenFOAM maps both native ranks through actual cell/face/point addressing; no serial-native global file is invented. The full all-state native audits and source hashes establish the record, with independent selected-frame native readback for visualization. Root separately inspected the actual ANSYS native velocity PNG and decoded paired assets.

See [source index](../results/annular-couette-v001/source-index.json) for exact public technical source versions and [selected case dictionaries](../cases/openfoam-annular-couette-selected-v001/README.md). Public code records the setup and executed workflow. The OpenFOAM preparation modules also require the external source-evidence manifest at import and the archived native preflight evidence, together with saved pilot and launch authorization records. This checkout is not a self-contained rerun package. Native HDF/mesh/NPZ, full fields/traces, raw audits/manifests, environments, licences and private launch state stay on E: outside GitHub content. The project's installed Fluent remains on D:.

## Contributions

ANSYS Fluent and OpenFOAM computed their own transient fields. No publication or university submission is performed by this record.
