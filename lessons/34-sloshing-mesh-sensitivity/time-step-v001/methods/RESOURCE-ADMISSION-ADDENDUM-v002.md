# Resource-wait admission, version002

The source/AST/test freeze for fullintervalv001 passed (1277inputs,15helpertests), but a read-only prelaunch probe found RAM below4GiB and correctly withheld before ANYfullcontroller/solver launch. RAM subsequently measured4.25GiB and was fluctuating. No nativefullstep/case was advanced. Keepv001 unlauncheddrafts/manifest; never call this a numericalmiss or a failedfullrun.

Use version002 fullcontrollers/wrapper/gate with tested wait_admission: every memory/process wait, including an immediate-gate resource race or normal Fluent cleanup, repeats source/parenthash checks and executable/process/RAM/disk admission. Unknown/source mismatch errors remainfatal. No lowerthreshold/no processbypass/no silenttimeextension. The Windows MPI gate retains explicit E:scriptpath and waits before MPI if resourceschangeaftercheckMesh. Two regressions passed afterRED: sourcechange duringwaitisrejected beforeprelaunch, andtransientresourcerace retriesallgates.
