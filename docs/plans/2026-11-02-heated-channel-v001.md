# Heated channel implementation plan — 2 November 2026

Previous goal turn made concrete progress: square-obstacle Fluent/OpenFOAM completed, native numerical checks passed, actual figures/GIF reviewed, complete local record verified. The overnight objective stays active.

1. Freeze matching geometry, properties, actual time integration, bounds and conservation targets in local/heated-channel-v001/TARGETS.md before CFD results. Preserve prior1017-file portfolio manifest.
2. Build and independently review the shared2,304-cell mesh. Prepare hidden, hash-bound two-process/rank Fluent and OpenFOAM pilots with actual steady flow followed by one frozen-flow thermal step only.
3. Recheck live resources, D: installation and licence/no competing solver. Run sequentially; save native flow, initial thermal and first actual0.005 s fields/reports/settings, then close. Independently check source mapping, flow convergence/mass/divergence/stability, actual frozen U/p/flux and all boundary energy semantics/storage. Preserve every failed setup or numerical miss.
4. Admit fixed600-step production only when first-step physics/units/frozen-flow checks are established. Save every actual thermal state with residual/bounds/energy evidence. Keep both native outputs separate from repository.
5. Render paired temperature evolution at61 matched actual times and a few academic static figures with units/shared scales/true channel aspect. Show actual Fluent image. Review all rendered PNGs and decoded GIF, verify sources and movement in both thermal maps.
6. Prepare complete local methods/code/compact results/limitations/contribution record, inspect before writing, preserve prior files/back up shared docs, verify links/figures/privacy/source hashes and refresh global manifest last. No Git commit/publication.
7. Continue the next evidence-supported model toward the overnight objective. No whole-curriculum launch and no claim that numerical agreement validates physics.
