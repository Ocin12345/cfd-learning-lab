import numpy as np
from local_refine_triangles_v002 import refine_mesh,edge
from triangle_quality_v002 import cell_quality

def repair_refinement(xy,tri,boundary,marked,centre,max_extent,screening_floor=.22,max_passes=20):
 marked=np.asarray(marked,dtype=bool).copy();history=[]
 adjacency={}
 for k,(a,b,c) in enumerate(tri):
  for e in [edge(a,b),edge(b,c),edge(c,a)]:adjacency.setdefault(e,[]).append(k)
 neighbours=[set() for _ in tri]
 for values in adjacency.values():
  if len(values)==2:a,b=values;neighbours[a].add(b);neighbours[b].add(a)
 for attempt in range(max_passes):
  q,t,b,p=refine_mesh(xy,tri,boundary,marked);affected=p['affected_parent_ids'];extent=float(np.linalg.norm(xy[tri[affected]]-centre,axis=2).max()) if affected else 0.
  if extent>max_extent:raise ValueError('Quality closure exceeded declared local extent')
  quality=cell_quality(q,t);bad=quality<screening_floor;history.append(dict(pass_index=attempt,marked_cells=int(marked.sum()),minimum_screening_quality=float(quality.min()),bad_children=int(bad.sum()),affected_extent=extent))
  if not bad.any():p['quality_closure_history']=history;p['screening_quality_minimum']=float(quality.min());p['screening_only_native_quality_still_required']=True;return q,t,b,p
  parents=np.unique(np.array(p['parent_ids'])[bad]);new=set(map(int,parents))
  # Red cells cannot improve their own shape again; include adjacent original
  # parents if the poor child is already red or if no new parent was found.
  for parent in parents:
   if marked[parent]:new.update(neighbours[parent])
  proposed=marked.copy();proposed[list(new)]=True
  if np.array_equal(proposed,marked):raise ValueError('No progress in bounded quality closure')
  marked=proposed
 raise ValueError('Quality closure pass budget exhausted')