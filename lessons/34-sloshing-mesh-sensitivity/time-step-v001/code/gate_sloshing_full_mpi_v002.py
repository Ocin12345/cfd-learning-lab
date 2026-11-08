import sys,json
sys.path.insert(0,'E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts')
from sloshing_wait_admission_v001 import wait_admission
print(json.dumps(wait_admission('full-interval-inputs-v002.json')))
