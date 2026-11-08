from pathlib import Path
import hashlib,json,subprocess
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001')
cmd=['wsl.exe','-d','CFD-Foundation14','-e','pgrep','-af','[f]oamRun|[s]impleFoam|[p]impleFoam|[m]pirun'];x=subprocess.run(cmd,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW);assert x.returncode in [0,1];assert x.returncode==1,'A solver is running; do not launch'
s=(S/'run_cylinder_pressure_pair_v002.py').read_text();s=s.replace('cylinder-diagnostic-pair-v001-status','cylinder-diagnostic-pair-v002-status').replace('diagnostic-controller-inputs-v001.json','diagnostic-controller-inputs-v002.json')
s=s.replace("['wsl','-d','CFD-Foundation14','--','pgrep'","['wsl.exe','-d','CFD-Foundation14','-e','pgrep'")
s=s.replace("capture_output=True,text=True,creationflags=subprocess.CREATE_NO_WINDOW)\n if x.returncode not in [0,1]","capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)\n if x.returncode not in [0,1]")
s=s.replace("'Cannot establish OpenFOAM process state: '+x.stderr","'Cannot establish OpenFOAM process state: '+repr(x.stderr)")
compile(s,'controller003','exec');(S/'run_cylinder_pressure_pair_v003.py').write_text(s)
a=json.loads((S/'diagnostic-controller-inputs-v001.json').read_text());a.pop(str(S/'run_cylinder_pressure_pair_v002.py'));a[str(S/'run_cylinder_pressure_pair_v003.py')]=hashlib.sha256((S/'run_cylinder_pressure_pair_v003.py').read_bytes()).hexdigest();(S/'diagnostic-controller-inputs-v002.json').write_text(json.dumps(a,indent=2)+'\n')
(S/'wsl-process-gate-regression-v001.json').write_text(json.dumps(dict(old_failure='WSL -- interpreted unquoted pipe-separated pattern through shell; localized UTF16 warning broke UTF8 text decoding. Controller failed before native output/session.',new_direct_exec=cmd,actual_no_solver_returncode=x.returncode,raw_stderr_bytes=len(x.stderr),passed=True),indent=2)+'\n')
with (S/'progress.md').open('a') as f:f.write('\nLaunchgatecorrection: controller002stoppedbeforecreatingnativecase/launchingsolver. WSL -- sentregexpipeviashell; localizedUTF16warningalsofailedUTF8decode. Actualold127commandfailurepreserved, direct -e pgrep bytes-modeGREENreturn1(noactiveprocess). Newcontroller003/manifest002/statusv002, no scientificlimit/sourcechange.\n')
print('Direct-exec WSL process gate: no active solver; controller003 ready')