import unittest
from wing_amg_transcript_v001 import parse_blocks
S='''pressure coupled equations:
 tol. 1e-17 2e-17 3e-17 4e-17
 0 4e-5 2e-5 3e-5 4e-5
 1 2e-8 1e-8 2e-8 3e-8
 1501 1.8e-5 2e-9 3e-9 4e-9 0:00:01 1

pressure coupled equations:
 tol. 1e-17 2e-17 3e-17 4e-17
 0 5e-5 3e-5 4e-5 5e-5
 1 3e-8 2e-8 3e-8 4e-8
 1502 2e-5 3e-9 4e-9 5e-9 0:00:00 0
'''
class Transcript(unittest.TestCase):
 def test_exact_outer_and_inner_identity(self):
  r=parse_blocks(S,1501,1502);self.assertEqual([x['iteration'] for x in r],[1501,1502]);self.assertEqual(r[0]['initial'],[4e-5,2e-5,3e-5,4e-5]);self.assertEqual(r[0]['terminal'],[2e-8,1e-8,2e-8,3e-8]);self.assertEqual(r[0]['cycles'],1)
 def test_missing_iteration_rejected(self):
  with self.assertRaises(ValueError):parse_blocks(S,1501,1503)
 def test_missing_initial_rejected(self):
  with self.assertRaises(ValueError):parse_blocks(S.replace(' 0 4e-5 2e-5 3e-5 4e-5\n',''),1501,1502)
 def test_bad_inner_index_rejected(self):
  with self.assertRaises(ValueError):parse_blocks(S.replace(' 1 2e-8',' 2 2e-8'),1501,1502)
 def test_duplicate_iteration_rejected(self):
  with self.assertRaises(ValueError):parse_blocks(S.replace(' 1502 ',' 1501 '),1501,1502)
 def test_nonfinite_rejected(self):
  with self.assertRaises(ValueError):parse_blocks(S.replace(' 1 2e-8',' 1 nan'),1501,1502)
 def test_unclosed_block_rejected(self):
  with self.assertRaises(ValueError):parse_blocks(S+'pressure coupled equations:\n 0 4e-5 2e-5 3e-5 4e-5\n',1501,1502)
if __name__=='__main__':unittest.main()
