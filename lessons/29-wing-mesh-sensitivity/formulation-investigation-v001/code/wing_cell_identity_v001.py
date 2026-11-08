import numpy as np,itertools
def centroid_bijection(native,exported,tolerance=2e-10):
 a=np.asarray(native,dtype=float);b=np.asarray(exported,dtype=float)
 if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3 or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('Invalid centroid arrays')
 grid=1e-8;keys=np.rint(a/grid).astype(np.int64);lookup={tuple(k):i for i,k in enumerate(keys)}
 if len(lookup)!=len(a):raise ValueError('Duplicate native centroid bins')
 index=[]
 for p in b:
  key=tuple(np.rint(p/grid).astype(np.int64));j=lookup.get(key)
  if j is None or np.linalg.norm(p-a[j])>tolerance:
   matches=[]
   for offset in itertools.product([-1,0,1],repeat=3):
    q=lookup.get(tuple(key[d]+offset[d] for d in range(3)))
    if q is not None and np.linalg.norm(p-a[q])<=tolerance:matches.append(q)
   if len(matches)!=1:raise ValueError('Unmatched/ambiguous exported centroid')
   j=matches[0]
  index.append(j)
 index=np.asarray(index,dtype=np.int64);distance=np.linalg.norm(b-a[index],axis=1)
 if len(np.unique(index))!=len(a) or np.max(distance)>tolerance:raise ValueError('Not a complete cell bijection')
 return index,distance
