from pathlib import Path
import hashlib,json
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001')
s=(S/'run_cylinder_pressure_pair_v003.py').read_text().replace('from cylinder_fluent_pressure_v002','from cylinder_fluent_coupled_v001').replace('diagnostic-controller-inputs-v002.json','coupled-controller-inputs-v001.json').replace('cylinder-diagnostic-pair-v002-status','cylinder-coupled-pilot-v001-status').replace("[('baseline',None),('tight',.001)]","[('coupled',None)]").replace("out=N/('fluent-pressure-'+label+'-v001')","out=N/'fluent-coupled-pilot-v001'").replace('cylinder_fluent_pressure_v002.py','cylinder_fluent_coupled_v001.py')
start=s.index(' # Control comparison:');end=s.index("if __name__=='__main__':",start)
s=s[:start]+''' baseline=N/'fluent-pressure-baseline-v001';selected=N/'fluent-coupled-pilot-v001'
 b=json.loads((baseline/'loaded-settings.json').read_text());c=json.loads((selected/'loaded-settings.json').read_text())
 for key in ['general','viscous','material','cell_zones','inlet','outlet','walls','references']:
  assert b[key]==c[key],'Unexpected physical input change: '+key
 for key in ['spatial_discretization','transient_formulation']:
  assert b['methods'][key]==c['methods'][key]
 assert c['methods']['p_v_coupling']['flow_scheme']=='Coupled'
 assert b['residual']==c['residual'],'Residual criteria/normalization changed'
 controls=json.loads((selected/'pressure-controls.json').read_text());assert controls['flow_courant_number']==20 and controls['explicit_pressure_under_relaxation']==controls['explicit_momentum_under_relaxation']==.5
 (S/'cylinder-coupled-comparison-v001.json').write_text(json.dumps(dict(results=results,physical_inputs_spatial_time_schemes_and_residual_settings_unchanged=True,algorithm_bundle=controls,full_extension_launched=False),indent=2)+'\\n')
 mark('coupled_pilot_complete',results=results,full_extension_launched=False)
''' +s[end:]
compile(s,'controller','exec');(S/'run_cylinder_coupled_pilot_v001.py').write_text(s)
a=json.loads((S/'cylinder-coupled-inputs-v001.json').read_text())
for name in ['run_cylinder_coupled_pilot_v001.py','audit_cylinder_window_v002.py']:
 p=S/name;a[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
for folder in ['fluent-pressure-baseline-v001','fluent-pressure-tight-v001']:
 p=Path('E:/CFD-Learning-Lab/native/runs/showcase-cylinder-convergence-v001')/folder
 for name in ['window-audit.json','loaded-settings.json']:
  f=p/name;a[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
(S/'coupled-controller-inputs-v001.json').write_text(json.dumps(a,indent=2)+'\n')
print('Coupled pilot controller frozen')