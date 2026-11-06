# Rotating-cylinder startup: methods

Both codes show accelerated upper-side flow and a slower asymmetric lower wake at 10 s.

Both codes computed 500 actual steps to 10 s and retained 501 native states. The [compact record](../results/magnus-cylinder-v001/summary.json) retains primary numerical outcomes and separate artifact diagnostics.

## Actual shared geometry

The prepared canonical mesh has 128 angular × 64 graded radial cells: 8,192 positive-area quads, 8,320 distinct nodes and 16,512 edges. It has one connected fluid region and one solid hole, with 128 cylinder facets, 64 upstream inlet facets and 64 downstream outlet facets. The inner and outer boundaries are inscribed 128-segment polygons with nominal radii 0.05 m and 2 m. The actual polygonal fluid area is 12.553473799926623 m²; cell areas, centroids and oriented facet vectors come from that geometry, rather than a circular-area substitute. Actual native mesh quality is included in the selected primary audits.

## Physics, time and separate native norms

Use constant density 1 kg/m³, dynamic viscosity mu = 0.0005 Pa s and kinematic viscosity nu = 0.0005 m²/s and inlet velocity (0.1, 0) m/s, giving Re_D = 20. The outlet gauge is 0 Pa. Prescribe absolute clockwise wall rotation Omega_z = −2 rad/s about the origin, with no mesh motion or MRF. Nominal wall speed is 0.1 m/s and spin ratio is 1. The native wall check uses the actual facet centre and tangent-projected Omega × r; the polygon's facet-centre speed differs slightly from the nominal radius value.

Each code continued five accepted pilot steps at dt = 0.02 s with 495 separately authorized additional steps: 500 solved steps, 501 saved states, physical time 0–10 s. Preserve original accepted pilot states 0–5 as authoritative. Do not reinitialize, replay or extend automatically.

Both codes use first-order Euler time integration and upwind momentum. Fluent uses SIMPLE with a 100-inner-iteration cap and own terminal continuity/x/y scaled residuals ≤ 1e−6. OpenFOAM uses 30 PIMPLE outer correctors, two pressure correctors and no nonorthogonal correctors; each solved step's actual final-outer initial p/U residuals must be ≤ 1e−6 and final linear residuals ≤ 1e−8. Retain the components actually printed. These native norms are different quantities.

## Initial state and predeclared solved-state checks

State 0 is prescribed cell rest with active inlet and rotating-wall boundary conditions. Retain its actual native fields and fluxes without projecting a manufactured conservative state. Initial mass/divergence/residual gates are excluded; an undefined initial ratio stays unavailable. Solved force histories have no invented initial row, while any raw native time-zero output is kept separately.

For every solved state, check native clocks, finite fields, material/BC controls, inlet velocity and wall motion. The outward inlet flux expectation is −0.4 kg/s at Fluent depth 1 m and −0.004 m³/s at OpenFOAM depth 0.01 m, with relative tolerance 1e−10. Cylinder native face flux must be ≤ 1e−12 SI. Relative all-boundary imbalance uses the actual inlet denominator and must be ≤ 1e−5. The **predeclared outlet reverse-flow ratio target is 1e−3**. Local all-face volume-flux divergence, scaled by D/U = 1 s, must have maximum ≤ 1e−5 and volume-weighted RMS ≤ 1e−6. Cell flux Courant number is descriptive, not a time-accuracy or implicit-stability claim.

## Native forces and surface pressure

Retain actual pressure, viscous and total x/y force coefficients at every solved step. Each code uses its actual depth: the coefficient denominator is 0.5 rho U² D depth, or 0.0005 N for Fluent and 0.000005 N for OpenFOAM. Compare coefficients, not the raw forces from different depths. Fluent saved wall shear is already integrated force and receives no additional area multiplier. Require independently checked native pressure/viscous/total endpoint reports at steps **5 and 500**, with coefficient agreement tolerance 1e−6 × max(1, absolute native coefficient). No lift-sign or ideal Cl target is imposed.

Final surface pressure uses 128 actual native facet centres, sorted by canonical source face IDs 16257–16384. Fluent native face IDs are one-based; OpenFOAM native global face indices are zero-based. Compute Cp = (p_surface − p_outlet) / q with p_outlet = 0 Pa and q = 0.005 Pa. Do not append fictitious angular endpoints or interpolate cell values onto invented wall points. If OpenFOAM omits a serialized zeroGradient pressure value, use and label its exact native owner-pressure semantics. Actual final facet pressure, centres and face mappings must be independently read back.

## Visuals, evidence and limits

The main GIF selects exactly steps 0:5:500 at 60 ms per frame. Its near-cylinder view is x = [−0.35, 0.65], y = [−0.30, 0.30] m, with equal xy aspect, the exact polygon hole and fixed speed scale 0–0.25 m/s. Label/count saturation. Preserve actual repeated maps; no temporal interpolation, invented particles or assumed signed vorticity. Physical duration 10 s and playback duration 6.06 s are distinct. A full-domain geometry overview supplies context. Final pressure maps share one explicitly stated range from both complete final native fields, with no mean subtraction.

Select actual accepted audits/manifests first, then verify selected raw U/p frames, native clocks/topology and final native surface pressures. Stream one native state and a proved mesh cache; keep full fields external. The compact records contain all 500 solved history rows, 16,384 final cell rows across both codes and 256 surface rows. Sources and actual visual QA were explicitly selected before this record.

Preserve any false derived-copy diagnostic and original/current hashes separately from primary numerical checks; native CFD completion, raw wrapper status and read-only audit outcome are distinct. Do not claim ideal-flow accuracy, periodicity, physical validation or mesh/time-step independence. This record authorizes no further CFD runtime or external publication.

## Reproduction and evidence

See [frozen targets](../results/magnus-cylinder-v001/TARGETS.md), [selected code](../results/magnus-cylinder-v001/source-index.json), [prerequisites](../results/magnus-cylinder-v001/OPENFOAM-PREREQUISITES.md) and [selected dictionaries](../cases/openfoam-magnus-cylinder-selected-v001/README.md). Public code records the setup and executed workflow. Native reruns and audit replay require the external project archive, native preflight/source evidence, accepted pilot fields and processor addressing, installed solver environments, and explicit launch authorization. This checkout is not a self-contained native rerun package. Native cases, meshes, fields, environments, licences, full traces and raw audit provenance remain external.

The user directed progression and visual priorities. Codex prepared and operated the scripts, checked native evidence and assembled the record. Fluent and OpenFOAM computed their own fields, forces and residuals. Independent user execution and owner reflection are unprovided.
