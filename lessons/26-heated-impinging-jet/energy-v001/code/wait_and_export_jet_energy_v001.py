import json,time,subprocess,traceback,os
from jet_energy_common_v001 import *
FIG=R/'exports/figures/2026-11-08-heated-jet-energy-full-v001';VIDEO=R/'exports/videos/2026-11-08-heated-jet-energy-v001';STATUS=S/'full-visual-pipeline-v001-status.json'
def verify_media():
 for p,h in json.loads((S/'media-inputs-v001.json').read_text()).items():assert sha(p)==h,p
def run(script,args,label):
 with (S/(label+'-v001.log')).open('xb') as f:p=subprocess.run([str(PY),'-B',str(S/script),*map(str,args)],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,(label,p.returncode)
def main():
 assert not STATUS.exists();verify_media()
 while True:
  a=json.loads((S/'production-pair-v001-status.json').read_text());stage=a['stage']
  if stage=='failed':raise RuntimeError('Native production controller failed; media withheld')
  if stage=='paired_full_600step_audits_pass':break
  mark(STATUS,'waiting_both_complete_native_audits',solver_stage=stage);time.sleep(30)
 verify_media()
 from jet_energy_full_prelaunch_v001 import wait_admission
 wait_admission(lambda **v:mark(STATUS,'waiting_media_resources',**v));verify_media()
 mark(STATUS,'rendering_actual_native_fields_and_allstep_energy');run('render_jet_energy_full_v001.py',[],'full-native-render');verify_media()
 mark(STATUS,'encoding_and_verifying_every_MP4_frame');run('export_jet_energy_video_v001.py',['--frames-dir',FIG/'frames','--out',VIDEO],'full-video-export');verify_media()
 mark(STATUS,'encoding_and_verifying_every_GIF_frame');run('export_jet_energy_gif_v001.py',['--frames-dir',FIG/'frames','--out',VIDEO/'gif-v001'],'full-gif-export');verify_media()
 mark(STATUS,'visuals_ready_for_actual_inspection',figures=str(FIG),video=str(VIDEO/'heated-jet-energy-comparison.mp4'),gif=str(VIDEO/'gif-v001/heated-jet-energy-comparison.gif'),whole_improvement_complete=False,fresh_review_not_yet_started=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
