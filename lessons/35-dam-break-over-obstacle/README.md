# 35. Dam break over an obstacle

![Actual paired development](full-production-v001/dam-break-development.gif)

[Play the MP4](full-production-v001/dam-break-development.mp4) · [Numerical evidence](full-production-v001/numerical-checks.png) · [Full record and methods](full-production-v001/README.md)

The selected ANSYS Fluent and OpenFOAM calculations cover **0.5 seconds, 10,000 steps and 101 native states per code** on9,072 cells. All selected steps/states pass their declared numerical checks. OpenFOAM is a saved-state continuation from a verified0.25s prefix; it is not a separate full rerun from zero.

Water crosses a low obstacle, rises into a jet and meets the far wall. The actual solver fields retain different jet/drop details. This is a coarse2D laminar incompressible model; mesh/time independence, impact-load accuracy and physical validation remain open.

The full record contains actual GIF/MP4, fixed-scale fields, pressure, convergence/mass/Courant plots, numerical CSVs, source copies, small inputs and preserved failures. Earlier qualification-stage files remain below unchanged. No owner reflection, independent execution or novelty is invented. This completed numerical record is prepared locally; publication is a separate step.
