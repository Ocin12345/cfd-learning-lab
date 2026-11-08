from pathlib import Path
import hashlib,json
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-cylinder-convergence-v001';F=R/'exports/figures/2026-11-08-cylinder-half-second-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
p=N/'fluent-coupled-window-v001';a=json.loads((p/'window-audit.json').read_text());assert a['numerical_window_pass'] and a['all_runtime_guards_pass'] and a['steps']==100 and a['saved_states']==101
for name,h in a['native_source_sha256'].items():assert sha(Path(name))==h
(S/'actual-half-second-image-review-v001.json').write_text(json.dumps(dict(reviewed_by='Codex',actual_native_figure=str(F/'cylinder-window.png'),sha256=sha(F/'cylinder-window.png'),visually_inspected=True,shown_in_chat=True,readable_axes_scales_no_clipped_values=True,force_history_not_assumed_developed=True),indent=2)+'\n')
(S/'DEVELOPED-INTERVAL-v001.md').write_text('''# Finite cylinder development interval

The accepted half-second continuation covers4.05..4.55s,100steps/101states, zeroresidualmisses, maxrelativeboundaryimbalance3.2761735869050655e-6(<1e-4), all100nativeHDFguards pass. Its actualfield/force/residualfigurewasinspected/shown.

Continue its exact4.55s checkpoint to absolute8s,690steps atdt0.005s,120innercap, unchangedCoupledbundle/physics/mesh/boundaries/spatialandtimeschemes/residualdefinitions. Eachnative state/force/massreport saved andcheckedbeforeadvancing. No extension beyond8s is admitted by this version. On anymisssave/stop/preserve; do notrelaxcriteria.

Combine accepted4..4.05,4.05..4.55 and4.55..8s only after all audits pass; exclude duplicate initialrows at segment boundaries. Expected800actualsteps/801uniquestates. Earlier0..4startup and numericalmisses remain inherited, not retrospectivelyrepaired. Quantifycyclewiseperiod, extrema, liftamplitude and drift from nativeforces. A settledwindow will be an observed metric, not assumed fromendtime or agreementinfrequency.

MatchingOpenFOAMcontinuationbegins only afterFluentacceptance, sequential2ranks. Practicalpressureoutlets/tractionequivalence, polygonalgeometry, inheritedstartup, mesh/timeindependence andexperimentalvalidation remainlimits. The publishedreference's longer development protocol andpressurepointdefinition are not silently replaced by this8srun. Anyadjacentcellpressureestimate must stateitscoordinates andsamplingmethod andmust notbe called exactbenchmarkwallpressure.
''')
s=(S/'cylinder_fluent_coupled_guarded_v001.py').read_text().replace('cylinder-coupled-window-inputs-v001.json','cylinder-developed-inputs-v001.json');(S/'cylinder_fluent_coupled_guarded_v002.py').write_text(s,encoding='utf-8')
s=(S/'run_cylinder_guarded_window_v001.py').read_text().replace('cylinder_fluent_coupled_guarded_v001','cylinder_fluent_coupled_guarded_v002').replace('guarded-window-controller-inputs-v001.json','developed-controller-inputs-v001.json').replace('cylinder-guarded-window-v001-status','cylinder-developed-fluent-v001-status').replace("parent=N/'fluent-coupled-pilot-v001';out=N/'fluent-coupled-window-v001'","parent=N/'fluent-coupled-window-v001';out=N/'fluent-developed-v001'").replace("a['steps']==10","a['steps']==100 and a['all_runtime_guards_pass']").replace("'--end-time','4.55'","'--end-time','8.0'").replace('start_time_s=4.05,end_time_s=4.55,steps=100','start_time_s=4.55,end_time_s=8.0,steps=690').replace("'cylinder-window-stdout.log'","'cylinder-developed-stdout.log'").replace("'cylinder-window-stderr.log'","'cylinder-developed-stderr.log'").replace('audit_window(out,4.05,100)','audit_window(out,4.55,690)').replace("'guarded_half_second_window_pass'","'fluent_development_interval_pass'").replace('steps=100,saved_states=101','steps=690,saved_states=691')
compile(s,'developedcontroller','exec');(S/'run_cylinder_developed_fluent_v001.py').write_text(s)
files=json.loads((S/'cylinder-coupled-window-inputs-v001.json').read_text());files.pop(str(S/'cylinder_fluent_coupled_guarded_v001.py'))
for n in ['cylinder_fluent_coupled_guarded_v002.py','actual-half-second-image-review-v001.json','DEVELOPED-INTERVAL-v001.md']:
 f=S/n;files[str(f)]=sha(f)
for name,h in a['native_source_sha256'].items():files[name]=h
for n in ['window-audit.json','code-provenance.json','geometry-setup.json','loaded-settings.json','final-settings.json','snapshots.json','force-history.csv']:
 f=p/n;files[str(f)]=sha(f)
(S/'cylinder-developed-inputs-v001.json').write_text(json.dumps(files,indent=2)+'\n')
for n in ['run_cylinder_developed_fluent_v001.py','audit_cylinder_guarded_window_v001.py','run_cylinder_pressure_pair_v003.py','audit_cylinder_window_v002.py','cylinder_fluent_pressure_v002.py','cylinder-developed-inputs-v001.json']:
 f=S/n;files[str(f)]=sha(f)
(S/'developed-controller-inputs-v001.json').write_text(json.dumps(files,indent=2)+'\n')
with (S/'progress.md').open('a') as f:f.write('\nTask2 halfsecond: ACCEPTED100steps/101states/100runtimeguards, maxrelativeboundaryimbalance3.27617e-6; all415nativesourcehashesrechecked. Actualfield/residual/forcesfigureinspected/shown. Ruling: finiteextensionfromexact4.55scheckpointto8s(690steps) withinexistingbound, thenmeasurecycledevelopment; noassumedbenchmarkpass orretroactivestartuprepair. Newsource002/receiptcylinder-developed-inputs001/controllerdeveloped001frozen.\n')
print('Finite4.55..8s continuation frozen')