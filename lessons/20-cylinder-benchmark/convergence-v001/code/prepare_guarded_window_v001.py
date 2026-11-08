from pathlib import Path
import hashlib,json
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-cylinder-convergence-v001';F=R/'exports/figures/2026-11-08-cylinder-convergence-diagnostic-v002'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a=json.loads((N/'fluent-coupled-pilot-v001/window-audit.json').read_text());assert a['numerical_window_pass'] and a['steps']==10 and a['saved_states']==11
for name,h in a['native_source_sha256'].items():assert sha(Path(name))==h
(S/'actual-coupled-diagnostic-image-review-v001.json').write_text(json.dumps(dict(reviewed_by='Codex',visually_inspected=True,shown_in_chat=True,figure=str(F/'cylinder-diagnostic.png'),sha256=sha(F/'cylinder-diagnostic.png'),notes=['Actual native topology/field, physical time and speed units clear.','All residual points visible with unchanged threshold; legend no longer covers data.','Native force traces nearly coincide; this does not establish accuracy or developed cycles.'],prior_plot_clipping_preserved=True),indent=2)+'\n')
(S/'COUPLED-ADDENDUM-v001.md').write_text('''# Bounded Coupled continuation

The pressure-termination-only comparison failed10/10newsteps in both variants; exactnumbers in pressure-hypothesis-result-v001.json. Coupledalgorithm/Courant20/explicitpressure-momentumrelaxation0.5 then passesall10steps and11states fromsame4s checkpoint atdt0.005s/120inner/unchanged1e-5/1e-4limits. Actualimagewasinspected/shown. Nativephysicalinputs, spatial/time schemes andresidualsettings matchthebaseline; no normalizedresidualreset.

Nextwindowcontinues its actual4.05scheckpoint to4.55s for100steps. Eachnewstep savesnativefield andforce/massreport, thenchecks its exactHDFnativecurrentclock, terminalraw/scale residuals, finiteU/V/P/density/viscosity andboundarymass beforeadvancingagain. Afailedstepissavedandstops; nofullclaim. Runtimeguardtests onrealaccepted/missed/wrongtimestates REDthenGREEN3/3. Originalstartup/parentmisses remain. SubsequentdevelopmentwindowandOpenFOAMrequirethiswindowpass andfreshsource/resource/processgates. No benchmarkpass orphysicalvalidation.
''')
s=(S/'cylinder_fluent_coupled_v001.py').read_text().replace('cylinder-coupled-inputs-v001.json','cylinder-coupled-window-inputs-v001.json')
s=s.replace('        completed = 0','        completed = 0\n        runtime_guard_records=[]')
needle="            export('forces_'+str(step), lambda: sample_forces(step, actual_time))"
assert s.count(needle)==1
s=s.replace(needle,needle+"\n            from cylinder_runtime_guard_v001 import check_native_step\n            guard=check_native_step(out/'snapshots'/f'time-{step:06d}.dat.h5',out/'force-reports'/f'mass-{step:06d}.txt',actual_time)\n            runtime_guard_records.append(guard)\n            write_json(out/'runtime-guard-records.json',runtime_guard_records)")
compile(s,'guardeddriver','exec');(S/'cylinder_fluent_coupled_guarded_v001.py').write_text(s,encoding='utf-8')
files=json.loads((S/'cylinder-coupled-inputs-v001.json').read_text());files.pop(str(S/'cylinder_fluent_coupled_v001.py'))
for name in ['cylinder_fluent_coupled_guarded_v001.py','cylinder_runtime_guard_v001.py','test_cylinder_runtime_guard_v001.py','runtime-green-v001.log','COUPLED-ADDENDUM-v001.md','actual-coupled-diagnostic-image-review-v001.json']:
 p=S/name;files[str(p)]=sha(p)
for name,h in a['native_source_sha256'].items():files[name]=h
for name in ['window-audit.json','code-provenance.json','geometry-setup.json','snapshots.json','force-history.csv','loaded-settings.json','final-settings.json']:
 p=N/'fluent-coupled-pilot-v001'/name;files[str(p)]=sha(p)
(S/'cylinder-coupled-window-inputs-v001.json').write_text(json.dumps(files,indent=2)+'\n')
with (S/'progress.md').open('a') as f:f.write('\nTask1 coupledpilot: PASS10/10newsteps,11states; physicalinputs/spatialandtimeschemes/residualsettings exactlymatchbaseline. Actualdiagnosticfigurev002inspected/shown; v001clippedrange/legendretained. Noteearlierpracticallyidenticalpressurehistorieswording: bothfail,buttightfinalcontinuity1.2519e-5versusbaseline1.2003e-5, exactdifferencesretained.\nTask2: nativeHDFtime/residual/material/massruntimeguardREDmissingimplementationthenGREEN3/3realnativepass/MISS/wrongtime. Guardeddriver001preparedfor4.05..4.55s100steps; acceptancepending. Alloldsources/casespreserved.\n')
print('Guarded0.5swindow prepared, not yet launched')