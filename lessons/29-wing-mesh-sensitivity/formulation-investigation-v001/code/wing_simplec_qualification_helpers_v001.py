import json,numpy as np
from pathlib import Path
from sloshing_fixed_settings_v001 import compare_fixed

def checkpoints():return [1520,1540,1560,1580,1600]

def scale_guard(iterations,raw_continuity,source_scale,end):
 i=np.asarray(iterations);v=np.asarray(raw_continuity)
 if not np.array_equal(i,np.arange(1,end+1)) or v.shape[0]!=len(i) or v.ndim!=2 or v.shape[1]<2 or not np.isfinite(v[:,:2]).all() or not np.all(v[:,1]>0):raise ValueError('Invalid actual residual history coverage')
 if end<=1500 or not np.all(v[i>1500,1]==source_scale):raise ValueError('Original continuity scale changed')
 return True

def query_iteration(st,out,label,expected):
 path=Path(out)/(label+'-live-iteration.json')
 if path.exists():raise FileExistsError(path)
 expressions=st.setup.named_expressions;before=expressions.get_state();name='codex-simplec-current-iteration'
 if name in expressions.get_object_names():raise ValueError('Temporary expression name already exists')
 try:
  expressions.create(name=name);expressions[name].definition='Iteration';value=expressions[name].get_value()
  with path.open('x',encoding='utf-8') as f:json.dump(dict(expression='Iteration',query_value=value,expected=expected),f,indent=2,allow_nan=False)
  number=float(value)
  if not np.isfinite(number) or number!=expected:raise ValueError('Actual live iteration differs from expected global endpoint')
 finally:
  if name in expressions.get_object_names():del expressions[name]
  compare_fixed(expressions.get_state(),before)
 return int(number)
