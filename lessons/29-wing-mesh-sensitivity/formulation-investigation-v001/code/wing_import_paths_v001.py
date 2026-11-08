import os,sys
from pathlib import Path
R=Path('E:/CFD-Learning-Lab');W=Path(__file__).parent
paths=[str(W),str(R/'repo/cfd-learning-lab/scripts'),str(R/'local/finite-wing-continuation-v001'),str(R/'local/finite-wing-mesh-sensitivity-v001'),str(R/'local/showcase-improvements-2026-11-08-v001/sloshing-time-step-v001')]
assert all(Path(p).is_dir() for p in paths)
sys.path[:0]=paths
os.environ['PYTHONPATH']=os.pathsep.join(paths)
