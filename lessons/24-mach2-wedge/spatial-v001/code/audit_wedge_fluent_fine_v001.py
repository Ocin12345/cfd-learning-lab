from pathlib import Path
import hashlib,json,math
import numpy as np
from wedge_endpoint_v001 import evaluate
from analyse_nozzle_fluent_v005 import residual_rows
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab')
def audit(run):
 assert json.loads((run/'status.json').read_text())['stage']=='closed'
 checkpoints=json.loads((run/'checkpoints.json').read_text());last,old=checkpoints[-2:]
 old,last=checkpoints[-2],checkpoints[-1]
 assert last['iteration']-old['iteration']==200 and last['iteration']<=1520
 row=residual_rows(run/last['data'])[-1];assert row['iteration']==last['iteration']
 for key in ('continuity','x-velocity','y-velocity','energy'):assert math.isclose(row[key],last['residuals'][key],rel_tol=1e-10,abs_tol=1e-20)
 with np.load(run/'geometry.npz') as geo:result=evaluate(run,last,old,geo)
 actual=json.loads((run/'loaded-settings.json').read_text());base=json.loads((R/'native/runs/mach2-wedge-fluent-v001/loaded-settings.json').read_text())
 assert all(actual[key]==base[key] for key in ('general','viscous','gas','inlet','outlet','methods')),'Fine actual physical/numerical setup differs from baseline'
 assert actual['controls']['courant_number']==3
 assert actual['viscous']['model']=='inviscid'
 wall=result['native_mass'];result['check_results']['wall_impermeability']=max(abs(wall[k]) for k in ('top','ramp'))/abs(wall['inlet'])<=1e-12
 paths=[run/n for n in ('geometry.npz','wedge.msh','loaded-settings.json','native-mesh-quality.json','status.json','checkpoints.json','executed-script.py')]+[run/c[k] for c in checkpoints for k in ('case','data','mass_report')]+[Path(__file__),S/'wedge_endpoint_v001.py',S/'plot_wedge_fluent_fine_v001.py']
 result.update(numerical_endpoint_pass=all(result['check_results'].values()),failed_checks=[k for k,v in result['check_results'].items() if not v],source_case=str(run),cells=89600,native_settings_equal_to_baseline=True,native_terminal_residuals=row,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},physical_validation=False)
 out=run/'endpoint-audit.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n')
 return result
