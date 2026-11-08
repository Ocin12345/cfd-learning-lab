import math

def compare_fixed(a,b):
 roundoff=[]
 def same(x,y,path):
  if isinstance(x,dict) and isinstance(y,dict):
   assert x.keys()==y.keys(),path
   for k in x:same(x[k],y[k],path+'/'+k)
  elif isinstance(x,list) and isinstance(y,list):
   assert len(x)==len(y),path
   for i,(v,w) in enumerate(zip(x,y)):same(v,w,path+'/'+str(i))
  elif isinstance(x,(int,float)) and not isinstance(x,bool) and isinstance(y,(int,float)) and not isinstance(y,bool):
   assert math.isclose(x,y,rel_tol=1e-14,abs_tol=1e-30),(path,x,y)
   if x!=y:roundoff.append(dict(path=path,actual=x,baseline=y))
  else:assert x==y and type(x)==type(y),(path,x,y)
 same(a,b,'');return dict(passed=True,roundoff=roundoff,relative_serialization_tolerance=1e-14)
