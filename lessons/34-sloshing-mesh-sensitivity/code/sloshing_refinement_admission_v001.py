"""Guarded resource and native-prerequisite admission for the fine-grid comparison."""
from pathlib import Path
import csv,hashlib,io,json,os,shutil,subprocess,sys,time
ROOT=Path('E:/CFD-Learning-Lab' if os.name=='nt' else '/mnt/e/CFD-Learning-Lab');LOCAL=ROOT/'local/sloshing-refinement-v001'
NAMES={'fluent-flat':'fluent-flat-audit.json','openfoam-flat':'openfoam-flat-audit.json','fluent-tilted':'fluent-tilted-audit.json','openfoam-tilted':'openfoam-tilted-audit.json'}
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
def fluent_running():
    r=subprocess.run(['tasklist.exe','/FO','CSV','/NH'],capture_output=True,text=True,encoding='utf-8',errors='replace',check=True,creationflags=subprocess.CREATE_NO_WINDOW)
    return any(row and row[0].lower().startswith('fluent') for row in csv.reader(io.StringIO(r.stdout)))
def wait_resources(stage):
    from wing_fluent_v001 import available_memory
    while True:
        ram=available_memory()/2**30;disk=shutil.disk_usage(ROOT).free/2**30;busy=fluent_running()
        if ram>=4 and disk>=10 and not busy:return dict(available_ram_GiB=ram,E_free_GiB=disk)
        print(json.dumps(dict(stage='waiting_for_resources',next_stage=stage,available_ram_GiB=ram,E_free_GiB=disk,fluent_session_detected=busy)),flush=True);time.sleep(30)
def source_path(name):
    text=name.replace('\\','/');prefix='E:/CFD-Learning-Lab/'
    assert text.startswith(prefix);p=(ROOT/text[len(prefix):]).resolve();assert p.is_relative_to(ROOT.resolve());return p
def verify_prerequisites(stage):
    receipt=json.loads((LOCAL/'prepared-inputs-check-v001.json').read_text())
    for source,expected in receipt['source_sha256'].items():assert sha(source_path(source))==expected,'Prepared input changed: '+source
    needed=[]
    if stage=='openfoam-flat':needed=['fluent-flat']
    elif stage in ['fluent-tilted','openfoam-tilted']:needed=['fluent-flat','openfoam-flat']+(['fluent-tilted'] if stage=='openfoam-tilted' else [])
    elif stage.endswith('production'):needed=list(NAMES)
    elif stage!='fluent-flat':raise ValueError(stage)
    hashes={}
    for key in needed:
        p=LOCAL/NAMES[key];a=json.loads(p.read_text());assert a['numerical_checks_pass'] is True,p
        assert a['cells']==12800 and a['steps']==40 and a['dt_s']==.0005
        assert a['case_mode']==('flat' if key.endswith('flat') else 'tilted')
        for name,expected in a['native_source_sha256'].items():assert sha(source_path(name))==expected,'Prerequisite source changed: '+name
        hashes[str(p)]=sha(p)
    baseline=ROOT/'local/sloshing-production-v001/accepted-pair-selection-v001.json'
    plan=json.loads((LOCAL/'PLAN.json').read_text());assert sha(baseline)==plan['coarse_selection_sha256'],'Coarse selection changed'
    if stage.endswith('production'):
        selection=json.loads(baseline.read_text())
        for name,key in [('fluent-production-audit-v001.json','fluent_audit_sha256'),('openfoam-production-audit-v003.json','openfoam_audit_sha256')]:
            p=ROOT/'local/sloshing-production-v001'/name;assert sha(p)==selection[key];a=json.loads(p.read_text());assert a['production_checks_pass']
            for source,expected in a['native_source_sha256'].items():assert sha(source_path(source))==expected,'Coarse native source changed: '+source
    hashes[str(baseline)]=sha(baseline)
    return hashes

if __name__=='__main__':print(json.dumps(verify_prerequisites(sys.argv[1])))
