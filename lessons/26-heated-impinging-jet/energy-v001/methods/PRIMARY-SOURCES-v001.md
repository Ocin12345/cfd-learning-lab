# Primary-source checks (8 November 2026)

ANSYS Fluent v261 User Guide, Fluxes Through Boundaries: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/flu_ug_FluxReport.html . Its total heat report includes conductive heat and convective enthalpy plus enabled work contributions. The accepted case disables the work terms. Treat its heat rates as reported integrated quantities, verify native face sums against actual reports, and shift enthalpy datum using the signed native mass flux.

OpenFOAM Foundation 14 scalarTransport source: https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/src/functionObjects/solvers/scalarTransport/scalarTransport.C . For volumetric phi without MULES, it solves transient transport plus divergence and diffusion, with the declared source/constraint machinery. Installed source was read directly; SHA256 f2997418fb3a3aa0e9869080f0108de0425421a1b4e537a8469200a5aaa96eb8.

Foundation 14 EulerDdtScheme source: https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-14/master/src/finiteVolume/finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.C . On this stationary mesh the matrix uses V/dt and the previous field. Installed source SHA256 5450404723d4d8089a7aa561b4ec4c01744d4dd941d0f998e4de762c293d356e. The online branch is a contextual reference, not an asserted byte match to the installed source.

The proposed sum dt*Q_n follows the configured backward-Euler discrete equation. It is not an exact integral of the continuous physical heat transfer, experimental validation, or time independence.
