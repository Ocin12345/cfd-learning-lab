"""Wait for numerically accepted results, then export actual-data figures and video."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time,traceback
R=Path('E:/CFD-Learning-Lab');P=R/'local/sloshing-refinement-v001';OLD=R/'local/sloshing-production-v001';STATUS=P/'visual-pipeline-v001-status.json';OWNED=False
FIG=R/'exports/figures/2026-11-07-sloshing-fine-paired-v001';CMP=R/'exports/figures/2026-11-07-sloshing-mesh-comparison-v001';VID=R/'exports/videos/2026-11-07-sloshing-fine-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mark(stage,**details):
    temp=STATUS.with_suffix('.tmp');temp.write_text(json.dumps(dict(stage=stage,observed_at=datetime.datetime.now().astimezone().isoformat(),**details),indent=2)+'\n');temp.replace(STATUS)
def resources():
    from wing_fluent_v001 import available_memory
    import shutil
    return available_memory()/2**30>=4 and shutil.disk_usage(R).free/2**30>=10
def call(args,name):
    receipt=json.loads((P/'visual-inputs-v002.json').read_text())
    for f,h in receipt['source_sha256'].items():assert sha(Path(f))==h,'Visualization input code changed: '+f
    with (P/(name+'.log')).open('x',encoding='utf-8') as log:
        r=subprocess.run([sys.executable,'-B',*map(str,args)],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,creationflags=subprocess.CREATE_NO_WINDOW)
    assert r.returncode==0,(name,r.returncode)
def main():
    global OWNED
    assert not STATUS.exists();OWNED=True
    assert not any(p.exists() for p in [FIG,CMP,VID])
    previews=set()
    while True:
        state=json.loads((P/'pipeline-v001-status.json').read_text());assert state['stage']!='failed',state
        for mode in ['flat','tilted']:
            audit=P/('fluent-'+mode+'-audit.json')
            if mode not in previews and audit.exists() and resources():
                a=json.loads(audit.read_text());assert a['numerical_checks_pass']
                out=R/('exports/figures/2026-11-07-sloshing-fine-'+mode+'-v002')
                mark('rendering_qualification',mode=mode);call([P/'show_fluent_qualification_v002.py',mode,out],'qualification-'+mode+'-render-v002');previews.add(mode)
                mark('waiting_for_full_pair',qualification_previews=sorted(previews))
        if state['stage']=='fine_pair_numerically_accepted' and resources():break
        time.sleep(20)
    mark('rendering_fine_pair');call([P/'render_fine_pair_v001.py','--openfoam-case',R/'native/runs/sloshing-refinement-v001/production/openfoam','--openfoam-audit',P/'openfoam-production-audit.json','--out',FIG],'fine-pair-render-v001')
    mark('comparing_matching_strips');call([P/'compare_coarse_fine_v001.py',FIG,CMP],'mesh-comparison-v001')
    mark('encoding_and_decoding_video');call([OLD/'export_sloshing_video_v001.py','--frames-dir',FIG,'--out',VID],'video-export-v001')
    mark('exporting_gif_preview');call([P/'export_fine_sloshing_gif_v001.py'],'gif-export-v001')
    mark('media_exported_review_pending',fine_figures=str(FIG),comparison_figures=str(CMP),video=str(VID/'sloshing-development.mp4'),numerical_checks_pass=True,physical_validation=False,repository_record_complete=False)
if __name__=='__main__':
    try:main()
    except Exception:
        if OWNED:mark('failed',error=traceback.format_exc())
        raise
