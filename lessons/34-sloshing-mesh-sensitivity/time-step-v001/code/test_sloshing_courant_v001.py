import unittest
from pathlib import Path
from sloshing_courant_v001 import courant_records
R=Path('E:/CFD-Learning-Lab/native/runs/sloshing-refinement-v001/production')
class Coverage(unittest.TestCase):
 def test_actual_complete_baseline_fluent(self):
  r=courant_records((R/'fluent/solver-transcript.txt').read_text(errors='replace'),'fluent',.0005,0,2400)
  self.assertEqual(len(r),2400);self.assertAlmostEqual(r[-1]['actual_pre_advance_time_s'],1.1995)
 def test_actual_complete_baseline_openfoam(self):
  r=courant_records((R/'openfoam/solver.log').read_text(),'openfoam',.0005,0,2400)
  self.assertEqual(len(r),2401);self.assertEqual([x['actual_pre_advance_time_s'] for x in r[:2]],[0,0])
 def test_missing_or_extra_report_rejected(self):
  t='Global Courant Number [Explicit VOF Criteria] : .01\nFlow time = .00025s, time step = 1\n'
  self.assertEqual(len(courant_records(t,'fluent',.00025,0,1)),1)
  for v in [t.replace('Global Courant Number','Absent'),t+'Global Courant Number [Explicit VOF Criteria] : .02\n']:
   with self.assertRaises(AssertionError):courant_records(v,'fluent',.00025,0,1)
 def test_continuation_reads_real_pre_advance_time(self):
  t='Interface Courant Number mean: 0 max: .01\nInterface Courant Number mean: 0 max: .01\nTime = 0.02025s\n'
  r=courant_records(t,'openfoam',.00025,80,81);self.assertEqual([x['actual_pre_advance_time_s'] for x in r],[.02,.02])
if __name__=='__main__':unittest.main()
