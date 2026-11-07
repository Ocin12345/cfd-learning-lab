"""Run one guarded background OpenFOAM case, leaving numerical audit separate."""
from pathlib import Path
import datetime,json,subprocess,sys,traceback
from dam_break_admission_v002 import verify,wait_resources
P=Path('E:/CFD-Learning-Lab/local/dam-break-over-obstacle-v001')
def main():
    stage=sys.argv[1];mapping={'openfoam-flat':'qualification/openfoam-flat-v002','openfoam-collapse':'qualification/openfoam-collapse-v002','openfoam-production':'production/openfoam-v002'};assert stage in mapping
    status=P/(stage+'-run-v002-status.json');assert not status.exists()
    def mark(state,**kw):
        temp=status.with_suffix('.tmp');temp.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),stage=stage,**kw),indent=2)+'\n');temp.replace(status)
    try:
        mark('waiting_for_resources');resources=wait_resources(stage);proof=verify(stage)
        mark('native_job_running',resources=resources,prerequisite_sha256=proof,ranks=2)
        out=P/(stage+'-wrapper-v002.log')
        with out.open('x',encoding='utf-8') as f:result=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/dam-break-over-obstacle-v001/run_dam_break_openfoam_v002.sh','/mnt/e/CFD-Learning-Lab/native/runs/dam-break-over-obstacle-v001/'+mapping[stage]],stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,creationflags=subprocess.CREATE_NO_WINDOW)
        assert result.returncode==0,(result.returncode,str(out))
        mark('computed_saved_audit_pending',numerical_acceptance=False,physical_validation=False)
    except Exception:
        mark('failed',error=traceback.format_exc());raise
if __name__=='__main__':main()