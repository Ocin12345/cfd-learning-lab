"""Final numerical and artifact checks after the paired plume renderer finishes."""
from pathlib import Path
import os,json,time,subprocess,traceback,hashlib,ast
from PIL import Image
R=Path('E:/CFD-Learning-Lab');L=R/'local/hot-patch-plume-production-v001';OUT=R/'exports/figures/2026-11-07-hot-patch-plume-20s-paired-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mark(stage,**kw):(L/'finish-checks-status.json').write_text(json.dumps(dict(stage=stage,**kw),indent=2)+'\n')
def main():
    receipt=L/'production-final-checks-v001.json';assert not receipt.exists()
    script=L/'audit_native_wall_flux_v001.py';ast.parse(script.read_text(encoding='utf-8-sig'))
    reference=json.loads((L/'native-wall-flux-reference-v001.json').read_text());assert reference['checks_pass']
    assert reference['native_source_sha256'][str(script)]==sha(script)
    script_hash=sha(script)
    mark('waiting_for_rendered_native_fields')
    while True:
        try:s=json.loads((L/'render-watch-status.json').read_text())
        except json.JSONDecodeError:time.sleep(2);continue
        assert s['stage']!='failed',s
        if s['stage']=='rendered_awaiting_visual_review':break
        time.sleep(10)
    assert sha(script)==script_hash,'Auditor changed after reference verification'
    mark('auditing_wall_heat_and_mass_flux')
    audit_path=OUT/'native-wall-audit.json'
    with (L/'final-wall-audit.log').open('w',encoding='utf-8') as log:
        p=subprocess.run([os.sys.executable,'-B',str(script),'--out',str(audit_path)],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,creationflags=subprocess.CREATE_NO_WINDOW,env=os.environ)
    assert p.returncode==0,p.returncode
    field=json.loads((OUT/'native-field-audit.json').read_text());wall=json.loads(audit_path.read_text())
    assert field['accepted_steps']=={'ANSYS Fluent':2000,'OpenFOAM':2000} and field['all_saved_field_checks_pass'] and wall['checks_pass']
    for audit in [field,wall]:
        assert all(sha(Path(p))==h for p,h in audit['native_source_sha256'].items()),'Native source hash changed'
    gif=OUT/'plume-development.gif'
    with Image.open(gif) as im:
        assert im.n_frames==41
        for i in range(im.n_frames):im.seek(i);im.load();assert im.info['duration']==250
    for name in ['paired-final.png','numerical-checks.png']:
        with Image.open(OUT/name) as im:im.verify()
    record=dict(physical_duration_s=20,steps_each=2000,numerical_checks_pass=True,native_hashes_verified=True,gif_frames_decoded=41,wall_checks_pass=True,final_comparison=field['comparisons'][-1],heat_transfer={n:{k:v for k,v in x.items() if k!='walls'} for n,x in wall['results'].items()},artifacts={p.name:sha(p) for p in OUT.iterdir() if p.is_file()},visual_review_pending=True,portfolio_record_pending=True,github_uploaded=False,physical_validation=False)
    receipt.write_text(json.dumps(record,indent=2)+'\n')
    mark('numerical_checks_complete_visual_review_pending',receipt=str(receipt),output=str(OUT))
try:main()
except Exception:
    mark('failed',error=traceback.format_exc());raise
