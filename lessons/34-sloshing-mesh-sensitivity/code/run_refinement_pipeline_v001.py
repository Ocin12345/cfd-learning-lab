"""Sequential, prerequisite/resource-guarded fine-grid qualification and production."""
from pathlib import Path
import datetime,json,os,subprocess,sys,traceback
from sloshing_refinement_admission_v001 import verify_prerequisites,wait_resources
R=Path('E:/CFD-Learning-Lab');P=R/'local/sloshing-refinement-v001';N=R/'native/runs/sloshing-refinement-v001'
STATUS=P/'pipeline-v001-status.json';OWNED=False
STAGES=[('fluent','flat','qualification/fluent-flat'),('openfoam','flat','qualification/openfoam-flat'),('fluent','tilted','qualification/fluent-tilted'),('openfoam','tilted','qualification/openfoam-tilted'),('fluent','production','production/fluent'),('openfoam','production','production/openfoam')]
def mark(stage,**details):
    temp=STATUS.with_suffix('.tmp');temp.write_text(json.dumps(dict(stage=stage,observed_at=datetime.datetime.now().astimezone().isoformat(),**details),indent=2)+'\n');temp.replace(STATUS)
def call(command,log,env=None):
    with log.open('x',encoding='utf-8') as f:
        result=subprocess.run(command,env=env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
    if result.returncode:raise RuntimeError(f'Child exited {result.returncode}; see {log}')
def run():
    global OWNED
    assert not STATUS.exists();OWNED=True
    assert Path(os.environ['FLUENT_EXECUTABLE']).is_file()
    for solver,mode,relative in STAGES:
        stage=solver+'-'+mode;out=N/relative;audit=P/(stage+'-audit.json')
        assert not audit.exists()
        if solver=='fluent':assert not out.exists()
        else:assert not (out/'solver.log').exists()
    for solver,mode,relative in STAGES:
        stage=solver+'-'+mode;out=N/relative;audit=P/(stage+'-audit.json')
        mark('waiting_for_resources',next_stage=stage)
        resources=wait_resources(stage);prerequisites=verify_prerequisites(stage)
        mark('solver_running_or_driver_resource_wait',current_stage=stage,output=str(out),ranks=2,resources=resources,prerequisite_sha256=prerequisites)
        if solver=='fluent':call([sys.executable,'-B',str(P/'run_refined_fluent_v001.py'),str(out)],P/(stage+'-solver.log'),dict(os.environ,CFD_FINE_MODE=mode))
        else:call(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/sloshing-refinement-v001/run_refined_openfoam_v001.sh','/mnt/e/CFD-Learning-Lab/native/runs/sloshing-refinement-v001/'+relative],P/(stage+'-wrapper.log'))
        mark('auditing_native_fields',current_stage=stage)
        call([sys.executable,'-B',str(P/'audit_refined_run_v001.py'),solver,str(out),str(audit),'--mode',mode],P/(stage+'-audit.log'))
        result=json.loads(audit.read_text());assert result['numerical_checks_pass']
        mark('stage_accepted',completed_stage=stage,audit=str(audit),next_stage_pending=True)
    mark('fine_pair_numerically_accepted',physical_time_s=1.2,comparison_visuals_complete=False,physical_validation=False)
if __name__=='__main__':
    try:run()
    except Exception:
        if OWNED:mark('failed',error=traceback.format_exc())
        raise
