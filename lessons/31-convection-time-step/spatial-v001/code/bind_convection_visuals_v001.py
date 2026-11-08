from pathlib import Path
import json,hashlib,datetime
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent;f=R/'exports/figures/2026-11-08-convection-spatial-full-v003';v=R/'exports/videos/2026-11-08-convection-spatial-v003'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
old=f.with_name('2026-11-08-convection-spatial-full-v001')
for name in ['field-and-heat-history.csv','mesh-response.csv']:
 assert (f/name).read_bytes()==(old/name).read_bytes()
a=json.loads((f/'selected-native-analysis.json').read_text());assert a==json.loads((old/'selected-native-analysis.json').read_text())
for path in [f/'native-field-audit.json',v/'video-export-record.json',v/'gif-v001/gif-export-record.json']:assert path.exists()
video=json.loads((v/'video-export-record.json').read_text());gif=json.loads((v/'gif-v001/gif-export-record.json').read_text());assert video['decoded_frames']==gif['decoded_frames']==41;assert video['maximum_rgb_rmse_0_255']<5 and gif['maximum_rgb_rmse_0_255']<8
files=[f/'frames/frame-000.png',f/'frames/frame-040.png',f/'coarse-fine-at20s.png',f/'mesh-response-histories.png',f/'full-numerical-checks.png',v/'decoded-020.png',v/'gif-v001/decoded-040.png']
receipt={'passed':True,'observed_at':datetime.datetime.now().astimezone().isoformat(),'actual_inspection':'Codex viewed actual initial/end PNGs, decoded MP4 middle and GIF end, four-case map and numerical plots. Uniform initial T, clockwise rising-warm-left circulation, fixed scales, unclipped readable axes, no central arrowhead blob in final display. Histories are unchanged from prior inspected version.','preferred_figures':str(f),'preferred_video':str(v/'cavity-refined-development.mp4'),'inspected_sha256':{str(p):sha(p) for p in files},'numeric_analysis_equal_to_preserved_v001':True,'csv_bytes_equal_to_v001':True,'all_frame_decode_pass':True,'source_bound_audit_sha256':sha(f/'native-field-audit.json'),'video_record_sha256':sha(v/'video-export-record.json'),'gif_record_sha256':sha(v/'gif-v001/gif-export-record.json'),'streamline_local_display_cutoff_mm_s':.5,'native_fields_modified':False}
(S/'convection-visual-review-v001.json').write_text(json.dumps(receipt,indent=2)+'\n');print('Visual receipt and numeric identity verified')
