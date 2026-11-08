import unittest,json,re
from pathlib import Path
P=Path('E:/CFD-Learning-Lab/repo/cfd-learning-lab/lessons/29-wing-mesh-sensitivity/formulation-investigation-v001')
class PressureEvidence(unittest.TestCase):
 def test_pressure_caption_and_methods_match_actual_binding(self):
  methods=json.loads((P/'checks/simplec-qualification-visual-binding-v002.json').read_text())['pressure_methods']
  self.assertEqual(methods['Fluent'],'saved native wing-face pressure');self.assertIn('owner-cell pressure proxy',methods['OpenFOAM'])
  text=(P/'README.md').read_text();i=text.index('![Actual saved wing pressure]');caption=text[i:text.index('\n\nThe accepted',i)] if '\n\nThe accepted' in text[i:] else text[i:i+650]
  self.assertRegex(caption,r'(?s)Fluent.*saved.*face.*OpenFOAM.*owner-cell pressure proxy')
  for name in ['METHODS.md','RESULTS.md']:
   t=(P/name).read_text();self.assertIn('owner-cell pressure proxy',t);self.assertIn('exact wall pressure',t)
if __name__=='__main__':unittest.main()
