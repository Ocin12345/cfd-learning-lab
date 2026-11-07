# Pressure coupling notes

The two initial attempts retain PISO pressure relaxation 0.3 and momentum relaxation 0.7 inherited from the steady control. The 200-inner attempt missed step 1. Increasing only the cap to 800 still missed y-momentum limits on steps 4, 5 and 9. The third attempt starts at the same t=0 state with pressure/momentum relaxation 1.0, PISO neighbor correction retained, the same cap and strict criteria. This is a numerical-method change, not a new physical condition.

Official guidance: https://ansyshelp.ansys.com/public/views/secured/corp/v252/en/flu_ug/flu_ug_sec_uns_solve_pvel_1.html

A possible later efficiency diagnostic is Fluent floating operating pressure for gradual pressure rise when heating gas in a closed domain. It has NOT been enabled in the current attempts. If tested, independently record the actual evolving operating pressure and use it for absolute-pressure/EOS comparisons; do not keep adding the original 101325 Pa to gauge pressure and do not derive an EOS check from the same EOS used to reconstruct pressure. Compare against saved qualified fixed-operating-pressure fields before accepting it.

Official scope and restrictions: https://ansyshelp.ansys.com/public/views/secured/corp/v252/en/flu_ug/flu_ug_sec_compressible.html

All references were checked on 2026-11-07. These are guidance for choosing numerical methods, not benchmark validation of this heater geometry.
