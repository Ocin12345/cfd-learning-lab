# Fluent vortex continuation v002: numerical correction

The preserved v001 startup pilot covered 0 to 2 s with dt = 0.005 s and a cap of 20 inner iterations per time step. All 400 native terminal records exceeded the declared 1e-5 residual target. At 2 s the continuity residual was 4.471299967178942e-5; x/y momentum residuals were 2.55532384625479e-7 and 1.0560008708284054e-7. The reported relative mass imbalance was 1.1021649292362917e-7. All 81 saved field snapshots and force samples passed finite-field, material, reference and physical-time checks. Those facts do not establish fully converged time stepping for the startup.

The v002 driver is a byte-for-byte copy of v001 apart from two substitutions:

- Increase MAX_INNER from 20 to 60. The dt = 0.005 s, PISO, second-order implicit temporal scheme, spatial schemes and 1e-5 residual target remain unchanged. This provides additional corrections per time step without weakening the declared target.
- Change the optional native vorticity PNG field from unsupported z-vorticity to supported vorticity-mag. This native picture represents vorticity magnitude. It is distinct from any signed curl derived from saved velocity fields for the movie.

The proposed continuation reads the saved v001 case/data at 2 s and advances to an absolute end time of 4 s in a fresh native/fluent-corrected-v002 directory. The existing 8 s absolute limit remains available but does not automatically authorize or require an extension beyond 4 s. No solver is launched by this preparation step.

The saved v001 results and sources remain unchanged. The continuation inherits the v001 state: additional corrections from 2 s onward do not retroactively repair the under-converged 0 to 2 s trajectory. The animation should focus on the corrected continuation, with the inherited startup limitation stated. Visual evidence of wake evolution is separate from a validated periodic shedding frequency or benchmark-accuracy claim.

The versioned v002 analyser checks the native maximum-inner-iteration setting against the prelaunch declared cap and summary, accepting the recorded v001 20 or v002 60 cap. The final current-step residual threshold remains 1e-5. Per-step terminal failures remain visible diagnostics, and no constant Cd/Cl criterion is imposed on this unsteady calculation.

## Frozen preparation sources

- vortex_fluent_v002.py: SHA-256 `76a2e46009915459da9cd9b267527bcb5d0dab787435022b3bbce45587ea1ce1`
- analyse_vortex_fluent_v002.py: SHA-256 `524e0f5517dfb54e024c345e8c219b8ac7115740e9df35b6f87f029c69c61478`

Preserved v001 driver SHA-256: `79e0fe095907119338d393c7f7f98519562e14b5e78273a6a66bab5d59e844d2`.
Preserved v001 analyser SHA-256: `bd61deaa7196ce6bebec0a8bb39f3c896a47c4b5a7186c15eb4bc265183a3a23`.
