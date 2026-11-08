import sys,json
sys.path.insert(0,'E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts')
from sloshing_prelaunch_v001 import prelaunch
print(json.dumps(prelaunch('full-interval-inputs-v001.json')))
