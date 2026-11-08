from pathlib import Path
import hashlib,json
import numpy as np
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');N=Path('E:/CFD-Learning-Lab/native/runs/showcase-cylinder-convergence-v001')
b=json.loads((N/'fluent-pressure-baseline-v001/window-audit.json').read_text());t=json.loads((N/'fluent-pressure-tight-v001/window-audit.json').read_text())
for a in [b,t]:
 for name,h in a['native_source_sha256'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==h
v=lambda a:np.array([[r[n] for n in ['continuity','x-velocity','y-velocity']] for r in a['step_records']])
dictout=dict(baseline_misses=b['residual_miss_steps'],tight_misses=t['residual_miss_steps'],baseline_terminal_residuals=v(b)[-1].tolist(),tight_terminal_residuals=v(t)[-1].tolist(),maximum_absolute_difference_terminal_history=float(np.max(np.abs(v(b)-v(t)))),maximum_mass_imbalance_each=[b['maximum_relative_mass_imbalance'],t['maximum_relative_mass_imbalance']],native_hashes_rechecked=True,hypothesis_rejected='Pressure AMG termination .1 to .001 alone does not remove continuity plateau.',next_hypothesis='Change pressure-velocity solution algorithm from segregated PISO to coupled block, Courant20/explicit pressure-momentum relaxation.5, same4s checkpoint/dt.005/120cap/physics/1e-5 target. This is an algorithm bundle, not isolated causality or retroactive startup repair.',source='https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/flu_th/flu_th_sec_uns_solve_pv.html')
(S/'pressure-hypothesis-result-v001.json').write_text(json.dumps(dictout,indent=2)+'\n')
s=(S/'cylinder_fluent_pressure_v002.py').read_text().replace('cylinder-inputs-v002.json','cylinder-coupled-inputs-v001.json').replace("methods.p_v_coupling.flow_scheme = 'PISO'","methods.p_v_coupling.flow_scheme = 'Coupled'")
start=s.index("        pressure_control=st.solution.controls.advanced.multi_grid.mg_controls['pressure']");end=s.index('        run = st.solution.run_calculation',start)
s=s[:start]+'''        pv=st.solution.controls.p_v_controls
        pv.flow_courant_number=20.
        pv.explicit_pressure_under_relaxation=.5
        pv.explicit_momentum_under_relaxation=.5
        if methods.p_v_coupling.flow_scheme()!='Coupled' or float(pv.flow_courant_number())!=20 or float(pv.explicit_pressure_under_relaxation())!=.5 or float(pv.explicit_momentum_under_relaxation())!=.5:raise ValueError('Coupled controls not retained')
        write_json(out/'pressure-controls.json',dict(original=None,actual=None,algorithm='Coupled',flow_courant_number=20.,explicit_pressure_under_relaxation=.5,explicit_momentum_under_relaxation=.5,controls=st.solution.controls.get_state(),same_parent=str(resume)))
''' +s[end:]
s=s.replace("coupling='PISO'","coupling='Coupled'");compile(s,'coupleddriver','exec');(S/'cylinder_fluent_coupled_v001.py').write_text(s,encoding='utf-8')
a=json.loads((S/'cylinder-inputs-v002.json').read_text());a.pop(str(S/'cylinder_fluent_pressure_v002.py'));a[str(S/'cylinder_fluent_coupled_v001.py')]=hashlib.sha256((S/'cylinder_fluent_coupled_v001.py').read_bytes()).hexdigest();a[str(S/'pressure-hypothesis-result-v001.json')]=hashlib.sha256((S/'pressure-hypothesis-result-v001.json').read_bytes()).hexdigest();(S/'cylinder-coupled-inputs-v001.json').write_text(json.dumps(a,indent=2)+'\n')
with (S/'progress.md').open('a') as f:f.write('\nTask1 pressurehypothesis: COMPLETE as preservednegativecomparison. Both10newstepsMISS1e-5continuity, baseline.1/tight.001 terminal histories practicallyidentical; fields/materials/nativeclocks/mass pass. Pressuretolerancealonehypothesisrejected; no extensionadmitted. Ruling: testCoupledsolutionbundle20Courant/.5explicitrelaxation on same4s checkpoint for10steps, preservingphysicalmodelandcriteria. Cost: shortdiagnostic, notuniquelyisolatedcause. Newalgorithmqualifiesonlyafterallstepspass.\n')
print(json.dumps(dictout,indent=2))