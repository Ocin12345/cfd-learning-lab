# Native energy evidence inventory

Both original codes contain 600 actual thermal step records, but only 32 distinct temperature states: 0:0.1:3 s plus 2.995 s. The original README's 31 frames denotes display frames; the extra penultimate state supports the final-step balance. Fluent has endpoint native heat and mass flux reports. OpenFOAM supplies frozen native phi and sparse temperature fields; its published heat rates are orthogonal/upwind reconstructions, not native boundary heat reports.

The available 2.995/3 s pair supports a final backward-Euler storage-rate check. It does not supply all 600 independent boundary heat samples. A curve interpolated through 0.1 s snapshots would be a quadrature approximation with unresolved intermediate variation. Neither its integral nor a boundary flux inferred from storage can be labelled an independently measured full discrete energy closure.

Ruling: add one bounded thermal-only replay with every-step temperature and Fluent native boundary heat/mass reports, retaining the exact accepted thermal initial fields and each code's own frozen flow. Original runs remain accepted and untouched. This repeats temperature evolution to collect missing evidence; it does not recompute or improve the flow. Cost: two finite thermal-only qualifications and, only if accepted, two finite exact-checkpoint continuations. The later wing convergence repair remains open.

Full idealization: constant properties, frozen velocity/pressure/flux, no buoyancy or viscous/pressure/kinetic work; 300 K jet, prescribed 320 K plate. The plate does not cool and temperature does not affect the flow.
