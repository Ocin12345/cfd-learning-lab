from pathlib import Path
import sys,json,datetime
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent;sys.path.insert(0,str(S));from jet_energy_common_v001 import sha,dump
F=R/'exports/figures/2026-11-08-heated-jet-energy-full-v001';B=R/'exports/figures/2026-11-08-heated-jet-energy-full-v002';V=R/'exports/videos/2026-11-08-heated-jet-energy-v001';selection=json.loads((S/'accepted-full-pair-selection-v001.json').read_text());bound={}
for code in ['fluent','openfoam']:
 p=Path(selection[code]['audit']);assert sha(p)==selection[code]['audit_sha256'];a=json.loads(p.read_text());assert a['numerical_pass'] and a['actual_steps']==600 and a['unique_states']==601;bound.update(a['source_sha256'])
for n in ['accepted-baseline-persistent-inputs-v001.json','full-interval-inputs-v001.json','media-inputs-v001.json']:bound.update(json.loads((S/n).read_text()))
for p,h in bound.items():assert sha(p)==h,p
a=json.loads((F/'native-field-audit.json').read_text());assert a['all_saved_field_checks_pass'] and a['plot']['times_s']==[s*.005 for s in range(0,601,10)]
for code,frames in a['native_frame_binding'].items():
 assert [f['step'] for f in frames]==list(range(0,601,10))
 assert [f['time_s'] for f in frames]==a['plot']['times_s']
for p,h in a['rendered_png_sha256'].items():assert sha(p)==h,p
for q in [F,B]:
 for p,h in json.loads((q/'native-field-audit.json').read_text())['static_png_sha256'].items():assert sha(p)==h,p
video=json.loads((V/'video-export-record.json').read_text());gif=json.loads((V/'gif-v001/gif-export-record.json').read_text());assert video['frames']==video['decoded_frames']==gif['frames']==gif['decoded_frames']==61
assert video['input_audit_sha256']==sha(F/'native-field-audit.json') and gif['source_native_field_audit_sha256']==sha(F/'native-field-audit.json');assert sha(V/'heated-jet-energy-comparison.mp4')==video['video_sha256'] and sha(V/'gif-v001/heated-jet-energy-comparison.gif')==gif['gif_sha256']
assert video['maximum_rgb_rmse_0_255']<5 and gif['maximum_rgb_rmse_0_255']<8
record=dict(passed=True,observed_at=datetime.datetime.now().astimezone().isoformat(),source_bindings_verified=len(bound),selection_sha256=sha(S/'accepted-full-pair-selection-v001.json'),figures_preferred=str(B),media=str(V),actual_inspected=[str(F/'frames/frame-000.png'),str(F/'frames/frame-030.png'),str(F/'temperature-comparison-3s.png'),str(V/'decoded-060.png'),str(V/'gif-v001/decoded-030.png')],static_figures_inspected={str(p):sha(p) for p in B.glob('*.png')},all_MP4_frames_decoded=61,all_GIF_frames_decoded=61,video_max_rgb_RMSE=video['maximum_rgb_rmse_0_255'],gif_max_rgb_RMSE=gif['maximum_rgb_rmse_0_255'],display_change='v002 fixes only residual title spacing; data/frame/video remain v001',observations=['Actual initial/mid/endpoint temperature evolution, equal color/domain scales and native static streamlines inspected','All wallheat, independentenergy and perstep residual plots inspected, axes/limits/labels readable','Video/GIF decoded representatives agree visually; everyframe source/RMSE/count receipts verified','No platecooling, thermalfeedback, physicalvalidation or continuous-integral inference'])
dump(S/'heated-jet-visual-review-v001.json',record);dump(S/'final-native-source-verification-v001.json',dict(passed=True,unique_bindings=len(bound),native_source_sha256=bound));print(json.dumps({k:record[k] for k in ['passed','source_bindings_verified','video_max_rgb_RMSE','gif_max_rgb_RMSE']}))
with (S/'progress.md').open('a',encoding='utf-8') as f:f.write('\nBOTHfullauditsaccepted600steps/601states, Ocompleted18:01:47. Existingmediawatchercompleted18:03:28,61actualframesMP4/GIFdecoded/sourcechecked. Actualinitial/mid/end/decodedrepresentativesandallstaticplotsinspected. Displayv002correctsmissingtitle spaceOnly; data/frames/videoUNCHANGEDv1. Postprocessingverification initiallycompared i*.05 against(10*i)*.005 withbitwisefloat equality; maxdifference4.44e-16, no nativeclockmiss. Correct canonicalintegerstep*dt metadata sequence checkedalongnativebindings; original1e-10nativeclockcriteriaunchanged. SourceverificationPASS'+str(len(bound))+'uniqueinputs; visualacceptancebound. No solver/media rerun, fullrecord/freshreviewpending.\n')
