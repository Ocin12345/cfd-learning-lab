"""Resource and immutable-source admission for sequential dam-break jobs."""
from pathlib import Path
import csv,hashlib,io,json,os,shutil,subprocess,sys,time
R=Path('E:/CFD-Learning-Lab' if os.name=='nt' else '/mnt/e/CFD-Learning-Lab');P=R/'local/dam-break-over-obstacle-v001'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def project_path(name):
    name=name.replace('\\','/');prefix='E:/CFD-Learning-Lab/';assert name.startswith(prefix);p=(R/name[len(prefix):]).resolve();assert p.is_relative_to(R.resolve());return p
def verify(stage):
    assert stage in [s+'-'+m for s in ['fluent','openfoam'] for m in ['flat','collapse','production']]
    receipt=P/'solver-inputs-v002.json';a=json.loads(receipt.read_text())
    for name,h in a['source_sha256'].items():assert sha(project_path(name))==h,'Input/source changed: '+name
    needed=[] if stage.endswith('-flat') else ['fluent-flat','openfoam-flat']
    if stage.endswith('-production'):needed+=['fluent-collapse','openfoam-collapse']
    proof={str(receipt):sha(receipt)}
    for key in needed:
        audit=P/(key+'-audit.json');a=json.loads(audit.read_text());assert a['numerical_checks_pass'] and a['cells']==9072 and a['steps']==200 and a['mode']==key.split('-',1)[1]
        for name,h in a['native_source_sha256'].items():assert sha(project_path(name))==h,'Prerequisite changed: '+name
        proof[str(audit)]=sha(audit)
    return proof
def busy():
    windows=subprocess.run(['tasklist.exe','/FO','CSV','/NH'],capture_output=True,text=True,encoding='utf-8',errors='replace',check=True,creationflags=subprocess.CREATE_NO_WINDOW)
    fluent=any(row and row[0].lower().startswith('fluent') for row in csv.reader(io.StringIO(windows.stdout)))
    foam=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','pgrep','-x','foamRun'],capture_output=True,timeout=30,creationflags=subprocess.CREATE_NO_WINDOW);assert foam.returncode in [0,1]
    return fluent,foam.returncode==0
def wait_resources(stage):
    from wing_fluent_v001 import available_memory
    while True:
        ram=available_memory()/2**30;disk=shutil.disk_usage(R).free/2**30;fluent,foam=busy()
        info=dict(stage=stage,available_ram_GiB=ram,E_free_GiB=disk,fluent_running=fluent,openfoam_running=foam)
        if ram>=4 and disk>=10 and not(fluent or foam):return info
        print(json.dumps(dict(waiting=True,**info)),flush=True);time.sleep(30)
if __name__=='__main__':print(json.dumps(verify(sys.argv[1])))