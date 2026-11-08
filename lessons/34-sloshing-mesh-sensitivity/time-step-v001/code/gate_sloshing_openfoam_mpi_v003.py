import sys
sys.path.insert(0,'E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts')
from sloshing_prelaunch_v001 import prelaunch
import json
print(json.dumps(prelaunch('qualification-inputs-v003.json')))
