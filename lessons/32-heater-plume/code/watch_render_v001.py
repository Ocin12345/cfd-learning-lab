"""Render only after the sequential production pipeline completes successfully."""
from pathlib import Path
import os,json,time,subprocess,traceback
R=Path('E:/CFD-Learning-Lab'); L=R/'local/hot-patch-plume-production-v001'
OUT=R/'exports/figures/2026-11-07-hot-patch-plume-20s-paired-v001'
def mark(stage,**kw):
    (L/'render-watch-status.json').write_text(json.dumps(dict(stage=stage,**kw),indent=2)+'\n')
try:
    assert not OUT.exists()
    mark('waiting_for_both_saved_runs')
    while True:
        try: state=json.loads((L/'pipeline-status.json').read_text())
        except json.JSONDecodeError:time.sleep(2);continue
        assert state['stage']!='failed',state
        if state['stage']=='both_runs_saved_ready_for_native_field_audit':break
        time.sleep(10)
    mark('auditing_and_rendering_native_fields')
    with (L/'production-render.log').open('w',encoding='utf-8') as log:
        p=subprocess.run([os.sys.executable,'-B',str(L/'render_native_plume_v002.py'),'--out',str(OUT)],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,creationflags=subprocess.CREATE_NO_WINDOW,env=os.environ)
    assert p.returncode==0,p.returncode
    audit=json.loads((OUT/'native-field-audit.json').read_text())
    assert audit['all_saved_field_checks_pass'] and audit['rendering']['frames']==41
    mark('rendered_awaiting_visual_review',output=str(OUT),github_uploaded=False)
except Exception:
    mark('failed',error=traceback.format_exc());raise
