import re,numpy as np

def phi_wall_flux(text,nwall,ninternal):
 d=re.search(r'dimensions\s*\[([^]]+)\]',text);assert d and np.array_equal(np.fromstring(d[1],sep=' '),[0,3,-1,0,0,0,0])
 m=re.search(r'internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
 if m:
  v=np.fromstring(m[2],sep=' ');assert int(m[1])==ninternal and v.shape==(ninternal,) and np.isfinite(v).all()
 else:
  m=re.search(r'internalField\s+uniform\s+([^;]+);',text);assert m and np.isfinite(float(m[1]))
 b=re.search(r'\btankWalls\s*\{([^{}]*)\}',text,re.S);assert b
 m=re.search(r'\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;',b[1],re.S)
 if m:
  v=np.fromstring(m[2],sep=' ');assert int(m[1])==nwall and v.shape==(nwall,)
 else:
  m=re.search(r'\bvalue\s+uniform\s+([^;]+);',b[1]);assert m;v=np.full(nwall,float(m[1]))
 assert np.isfinite(v).all();return float(abs(v).sum())
