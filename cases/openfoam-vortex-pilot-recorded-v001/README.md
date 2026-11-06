# OpenFOAM vortex pilot: recorded inputs v001

PILOT UNDERCONVERGED. The completed 2 s startup pilot missed the terminal-initial p/U residual target of 1e-5. The recorded check keeps that criterion false.

Fixed dt 0.005 s; backward time scheme; laminar PIMPLE with 3 outer and 2 pressure corrections. Saved fields every 0.025 s and force coefficients every 0.005 s.

Shared model: 2.2 x 0.41 m channel, cylinder centre (0.2, 0.2), diameter 0.1 m; rho=1 kg/m3 and nu=0.001 m2/s. The parabolic inlet has continuum mean 1 m/s and maximum 1.5 m/s. Force references use U=1 m/s, A=0.1 m2 and L=0.1 m. The eight dictionaries below are exact byte copies from `native/runs/2026-11-01-vortex-shedding-v001/openfoam`.

The read-only v002 pilot checks confirm 400 physical time steps, finite final fields, and relative inlet/outlet mass imbalance 3.115e-11. These checks establish neither periodic shedding nor benchmark validation.

Use the matching external shared mesh and these recorded initial fields to reproduce the startup input state.

Shared mesh SHA256: `e2f209fa90d32f5c9967f1733c2b2b9b10054e3bb3b408a8ed6a5f6addcbb837`.

Pilot evidence: local/vortex-shedding-v001/openfoam-pilot-checks-v002/checks.json. Its source manifest confirms unchanged native inputs during analysis. Dictionary-copy verification is retained in `local/vortex-shedding-v001/openfoam-recorded-dictionaries-v001.json`.

| Dictionary | SHA256 |
| --- | --- |
| 0/U | `602bfad0c9ee7708445413db7e4b2e25e996e046e235bea7a9f59fc5876eb1a0` |
| 0/p | `bb9817522929bdaa6ebb2239774819ef994c6f426ef34a85545a096e1dd16a56` |
| constant/physicalProperties | `0c2cccd25f6ff7637351833f35ea078390e21633c99246122dc74f495fa94b0b` |
| constant/momentumTransport | `7498c38a3be52e135c3f5a82c90cc81cbe88a784f678ba1fb5fa24483c164844` |
| system/controlDict | `aea20f578218449b3b129bf1c1618b5637cdf23f1ce0652da78f341e42ceaab0` |
| system/fvSchemes | `10527a4b72853f3e8fa9d7b15643bac177e2f4df02214473f8a8f821fafc827c` |
| system/fvSolution | `ce259a6bd1ec9e8507155394be8d367d850146c0c4f0c73d4efa10578cae9349` |
| system/decomposeParDict | `5c4d636cb3070a405d3437172ffbf73a3104582ae6347170c91cc11666011797` |
