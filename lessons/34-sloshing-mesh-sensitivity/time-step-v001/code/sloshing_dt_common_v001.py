from pathlib import Path
from sloshing_refinement_common_v001 import NX,NY,DX,L,H,MATERIALS,WALLS,EXPECTED_AREA,RESIDUAL_TARGET,INITIAL_TILT,RAW_TARGETS,write_json,alpha_group,cell_values
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent;N=R/'native/runs/showcase-sloshing-time-step-v001';B=R/'native/runs/sloshing-refinement-v001/production'
DT=.00025
