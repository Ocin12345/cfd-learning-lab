def configure_residuals(residual,check):
 for name in residual.equations.get_object_names():
  equation=residual.equations[name]
  equation.check_convergence=True
  equation.absolute_criteria=1e-9 if name=='energy' else 1e-6
  equation.check_convergence=bool(check)