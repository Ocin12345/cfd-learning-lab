# OpenFOAM pressure-coupling budget amendment v001

Declared after the preserved five-step10-outer pilot, before any new pilot run. Same physical annulus, mesh, rho/nu/Omega, rest initialization, Euler dt0.05s, upwind momentum and all numerical thresholds. Canonical TARGETS-v002.md remains preserved and unchanged. This amendment supersedes only its OpenFOAM outer-corrector iteration budget for the next fresh case.

The actual original pilot native log finished five steps to0.25s with10PIMPLEoutercorrectors each. At step5, the lastouter pressure-correction initial residuals were7.4664243934238e-5 and5.6044434716295e-5, exceeding the declared1e-6 threshold. Terminal matrix residuals were small; they do not establish outer coupling convergence. Full original audit/false checks and all actual native states remain preserved.

Next fresh OpenFOAM pilot: fixed30outercorrectors,2pressurecorrectors,0nonorthogonalcorrectors; otherwise byte-identical physical/numerical dictionaries except generated provenance. No tolerance is relaxed. Start from copied/proven original native initial-rest fields at time0, never from the insufficiently converged step5. Exactlyfivephysicalsteps to0.25s. If all unchanged native/source/physics/residual/continuity gates pass, separately admit395additionalsteps to20s/400total, noautomaticextension. Geometry conversion and time0runtime-registration commands need not be replayed when their actual mesh/initial fields/mappings/installed-source evidence are copied and byte verified.

Original10outer pilot is one numerical attempt, not a separate physical model or validation. Fluent's independently accepted5step pilot and fixed100inner budget are unaffected. No paired pass or experimental validation claim is introduced.
