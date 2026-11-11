# Actual passive-tracer figures

The concentration GIF shows 101 actual native states at steps 0:8:800, separated by 0.2 s physical time, with 60 ms playback delay. Flow remains frozen while concentration evolves. Both panels use fixed c=0–1 limits, equal physical axes and actual T-domain cells. Under/over colors expose saturation. No animated particles, interpolated physical states, forced map changes or clipped data are added. Initial c=0 has no solved residual or Euler budget.

The final map is 20 s; snapshots select 2, 6, 12 and 20 s. Transport histories show per-depth mass-weighted dimensionless inventory, flux-weighted outlet concentration and cumulative boundary transfers. Advection/diffusion, storage and Euler defects are RECONSTRUCTED from native fields and fluxes, not saved native matrix flux. Inventory in kg/m is not calibrated dye mass. Distinct native residual norms occupy separate panels; logarithmic drawing omits exact zeros only, retaining them in CSVs. The contact sheet was decoded from the actual GIF.

CSV coordinates/areas describe canonical plotting polygons. NPZ geometry retains authenticated native coordinates/areas; scalar ordering, mappings and selected native values were checked independently. QA verified all 7,048,800 full cell rows, 8,800 final rows, 1,602 history rows, 101 GIF frames and source/output hashes.

Each solver passed 11 primary numerical checks. Interrupted archive/reconstruction status, missing producer records and raw pending flags remain preserved separately from independent audit/media acceptance; no CFD was replayed. First-order smearing and the finite 20 s window limit interpretation. No complete washout, physical validation or mesh/time independence is established.
