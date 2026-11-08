import unittest,json,numpy as np
from cavity_json_scalars_v001 import scalar_default
class Scalars(unittest.TestCase):
 def test_numpy_scalars_preserve_values(self):
  q=json.loads(json.dumps({'ok':np.bool_(False),'step':np.int64(2000),'value':np.float32(.5)},default=scalar_default));self.assertIs(q['ok'],False);self.assertEqual(q['step'],2000);self.assertEqual(q['value'],.5)
 def test_unknown_not_silently_stringified(self):
  with self.assertRaises(TypeError):json.dumps({'x':object()},default=scalar_default)
if __name__=='__main__':unittest.main()