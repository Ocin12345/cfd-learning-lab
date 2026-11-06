# Project errata

## 1 November 2026: R40 native Fluent pressure caption

The preserved [v001 native-export caption](../images/2026-11-01-airfoil-domain40-fluent-native-v001/caption.md) incorrectly described the Fluent **operating pressure as 0 Pa**. It confused operating pressure with the **0 Pa gauge pressure reference**. The saved R40 continuation and its solver settings use **operating pressure 101325 Pa**, **outlet gauge pressure 0 Pa** and **force/plot reference pressure 0 Pa gauge**. The [corrected v002 caption and unchanged PNG copies](../images/2026-11-01-airfoil-domain40-fluent-native-v002/caption.md) preserve the original labels, interpolation and displayed ranges. The v001 images, caption and manifest remain unchanged as historical evidence.

This is a documentation correction. The actual native case/data, material, boundary conditions, meshes, iteration histories and measured forces were not changed. The versioned [v002 source manifest](../images/2026-11-01-airfoil-domain40-fluent-native-v002/source-manifest.json) binds the native saved case/data, solver readback `loaded-settings.json`, v001 files and byte-identical output images.

The physical material density is 1 kg/m3 and dynamic viscosity is **0.01 Pa s**. With inlet velocity 1 m/s and chord 1 m, the calculation therefore remains **Re = rho U c / mu = 100**. The separate legacy reference-viscosity entry is 1.7894e-5 Pa s in the original medium, fine and R40 Fluent records; it is not the viscosity used in the saved flow fields or the stated Reynolds number. Numeric force references remain density 1 kg/m3, speed 1 m/s, area 1 m2 and planar depth 1 m.

The reported force coefficients and computed pressure-coefficient figures are unaffected by correcting this caption. Force coefficients use `F / (0.5 rho U^2 A)`; the field plotting script explicitly computes `Cp = (p_gauge - p_ref,gauge) / (0.5 rho U^2)` with `p_ref,gauge = 0 Pa`. Neither expression uses the caption text or the legacy reference-viscosity entry. The original solver rendering is static gauge pressure, not Cp. No new flow iterations or image edits were performed for this correction.

These records describe numerical checks and measured sensitivity. They do not establish physical validation, mesh/domain independence, aircraft performance or independently authored owner work.
