import unittest,numpy as np
from wedge_width_columns_v001 import measure_columns
class Columns(unittest.TestCase):
 def sample(self):
  xx=np.linspace(.35,.85,11);yy=np.linspace(0,.82,165);x,y=np.meshgrid((xx[:-1]+xx[1:])/2,(yy[:-1]+yy[1:])/2)
  slope=.8;b=-.16;d=(slope*x+b-y)/np.sqrt(1+slope*slope);q=np.clip(.5+d/.04,0,1)
  geo=dict(y_nodes=np.zeros((165,11)),centres=np.stack([x,y],axis=-1).reshape(-1,2))
  metrics=dict(shock_fit_slope=slope,shock_fit_intercept_m=b,pressure_upstream_pa=1.,pressure_downstream_pa=2.)
  return 1+q,geo,metrics
 def test_linear_width_and_coverage(self):
  p,g,m=self.sample();a=measure_columns(p,g,m);self.assertEqual(a['valid_columns'],10)
  self.assertAlmostEqual(a['median_width_m'],.032,places=12)
  for row in a['columns']:
   self.assertEqual(len(row['bracket_samples']),6);self.assertGreater(row['width_in_projected_column_intervals'],0)
 def test_no_silent_column_drop(self):
  p,g,m=self.sample();p[:,4]=1.
  with self.assertRaisesRegex(ValueError,'column 4'):measure_columns(p,g,m)
if __name__=='__main__':unittest.main()