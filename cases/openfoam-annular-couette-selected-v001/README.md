# Selected OpenFOAM annular Couette dictionaries

These compact dictionaries record the actual selected fixed-mesh, laminar annulus and 30-outer pressure-coupling budget. The shared 64 by 24 quad mesh comes from the public mesher; native mesh/rank files remain external. See the [source index](../../results/annular-couette-v001/source-index.json) for exact public mesher, preparation, Fluent and OpenFOAM driver/analyser versions, wrapper sources and plotting dependencies.

Public code records the setup and executed workflow. The OpenFOAM preparation modules also require the external source-evidence manifest at import and the archived native preflight evidence, together with saved pilot and launch authorization records. This checkout is not a self-contained rerun package.

A fresh reproduction must generate the mesh, verify native conversion/quality and moving-wall units/reference, run and audit a five-step pilot, then separately authorize 395 additional steps. This directory is a configuration record, not a standalone archived restart. The guarded production drivers require the external actual accepted pilot, its native state files and a root authorization receipt; those private launch records and large native results are not copied here. Do not present a missing external restart as runnable from these dictionaries alone.

Native rank U, p, face flux, addressing and clocks underlie all 401 actual saved states. No serial global field is invented. Public histories retain 400 solved rows per code, with component labels and their distinct native residual norms. Initial rest has no solved residual row. [Methods](../../docs/ANNULAR_COUETTE_METHODS.md) state units, assumptions, checks, limitations and contributions.
