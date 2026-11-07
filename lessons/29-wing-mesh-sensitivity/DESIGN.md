# Finite wing mesh sensitivity
Same CAD/profile and 8-degree incidence, Re400, rho1, mu0.0025, U=(1,0,0), area4.
Baseline169800 cells. New target mesh sizes are0.9 times original; hard cap250000.
Candidate native volume-quality proxy margin0.105; final Fluent minimum orthogonal quality0.10.
One failing candidate tetrahedron may be repaired by bounded interior-node relocation using the existing preservation-checked routine. No boundary or connectivity changes.
Run sequential hidden two-process Fluent coupled(CFL20,pseudooff,explicit0.5) and two-rank OpenFOAM SIMPLEconsistent; maximum1500 iterations each.
Residuals<=1e-5, relative mass imbalance<=1e-4, final100 force range/max(abs(mean),1e-3)<=1e-3.
No numerical miss can be relabelled as a pass. No mesh independence or physical validation claimed from a two-grid comparison.
