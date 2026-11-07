"""Same wing mesh, stream aligned with its clockwise-eight-degree chord."""
import math
U=[math.cos(math.radians(8.)),-math.sin(math.radians(8.)),0.]
LIFT=[-U[1],U[0],0.]
assert abs(sum(x*x for x in U)-1)<1e-14
assert abs(sum(a*b for a,b in zip(U,LIFT)))<1e-14
MAX_ITERATIONS=1200
