"""Exact restart equality and complete discrete-step merge; no interpolation."""
import numpy as np

def merge_steps(prefix,continued):
 if [r['step'] for r in prefix]!=list(range(1,21)) or [r['step'] for r in continued]!=list(range(21,601)):raise ValueError('Incomplete, duplicated or reset step sequence')
 rows=prefix+continued
 for r in rows:
  if not np.isfinite(r['time_s']) or abs(r['time_s']-r['step']*.005)>1e-10:raise ValueError('Native global time mismatch')
 return rows

def same_restart(before,after):
 if set(before)!=set(after):raise ValueError('Restart field names differ')
 for n in before:
  a=np.asarray(before[n]);b=np.asarray(after[n])
  if not np.isfinite(a).all() or not np.isfinite(b).all() or not np.array_equal(a,b):raise ValueError('Restart changed: '+n)
 return True
