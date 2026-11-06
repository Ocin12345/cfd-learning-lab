# OpenFOAM vortex corrected continuation: recorded inputs v003

NUMERICAL POSTCHECKS PASSED AT 4 s. These are the actual dictionaries loaded for the authorized continuation from the retained 2 s checkpoint. All declared v002 numerical checks passed; the analyser returned exit code 0.

Fixed dt 0.005 s; backward time scheme; laminar PIMPLE with 8 outer and 2 pressure corrections. Continuation fields are saved every 0.05 s; force coefficients remain sampled every 0.005 s.

Shared model: 2.2 x 0.41 m channel, cylinder centre (0.2, 0.2), diameter 0.1 m; rho=1 kg/m3 and nu=0.001 m2/s. The parabolic inlet has continuum mean 1 m/s and maximum 1.5 m/s. Force references use U=1 m/s, A=0.1 m2 and L=0.1 m. The eight dictionaries below are exact byte copies from `native/runs/2026-11-01-vortex-shedding-v001/openfoam-extended`.

The corrected continuation uses a separate verified clone. The original pilot dictionaries and results remain preserved. At 4 s, terminal initial residuals are Ux=5.73485455919e-09, Uy=1.57455620707e-08, p=1.22731151267e-08, each below 1e-5. Independent native-phi relative inlet/outlet mass imbalance is 3.82538121839e-11, below 1e-4. Final U, p, C and phi fields are finite.

For restart, use the matching external mesh and retained native 2 s checkpoint. The copied 0/U and 0/p describe inherited startup inputs; controlDict selects latestTime for continuation.

Shared mesh SHA256: `e2f209fa90d32f5c9967f1733c2b2b9b10054e3bb3b408a8ed6a5f6addcbb837`.

The clone was verified against 1,521 source files (546,087,911 bytes): local/vortex-shedding-v001/openfoam-clone-verification-v001.json. The native continuation provenance archives pilot-fvSolution before changing the outer correction count. Dictionary-copy verification is retained in `local/vortex-shedding-v001/openfoam-recorded-dictionaries-v001.json`.

The saved data contain 120 nonzero-time snapshots: 80 pilot snapshots through 2 s at 0.025 s intervals, then 40 through 4 s at 0.05 s intervals. Force data contain 801 samples including t=0. Actual logs confirm 400 steps with 3 outer corrections followed by 400 with 8; all 400 corrected steps meet the terminal-initial residual target. The 400 pilot misses remain recorded. These startup checks do not establish mesh independence, periodic shedding or benchmark validation.

Postcheck evidence: `local/vortex-shedding-v001/openfoam-corrected-checks-v002/checks.json`, its `source-manifest.json`, and `openfoam-corrected-transition-checks-v001.json`. The manifests bind the actual v003 helper and executed wrapper hashes and confirm unchanged source inputs.

| Dictionary | SHA256 |
| --- | --- |
| 0/U | `602bfad0c9ee7708445413db7e4b2e25e996e046e235bea7a9f59fc5876eb1a0` |
| 0/p | `bb9817522929bdaa6ebb2239774819ef994c6f426ef34a85545a096e1dd16a56` |
| constant/physicalProperties | `0c2cccd25f6ff7637351833f35ea078390e21633c99246122dc74f495fa94b0b` |
| constant/momentumTransport | `7498c38a3be52e135c3f5a82c90cc81cbe88a784f678ba1fb5fa24483c164844` |
| system/controlDict | `e295c32c212d2120acdb1eabaf5338f34b14c453dbfee716ae415ac7cbbf75db` |
| system/fvSchemes | `10527a4b72853f3e8fa9d7b15643bac177e2f4df02214473f8a8f821fafc827c` |
| system/fvSolution | `74d8143df5498fcf84cc4aa9bf7f4a58f8d763b1843ec41b8e0fa2344b1eb009` |
| system/decomposeParDict | `5c4d636cb3070a405d3437172ffbf73a3104582ae6347170c91cc11666011797` |
