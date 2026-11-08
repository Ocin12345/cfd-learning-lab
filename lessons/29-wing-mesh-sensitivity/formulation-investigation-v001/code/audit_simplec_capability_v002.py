import wing_import_paths_v001
from wing_repair_common_v001 import *
from wing_simplec_prelaunch_v002 import verify
from wing_repair_native_v003 import raw_residual,same_native,state_keys
from wing_repair_helpers_v001 import same_residual_prefix
from wing_simplec_capability_helpers_v002 import validate_state,diff_inventory,assert_zero_advance_source
from sloshing_fixed_settings_v001 import compare_fixed
import numpy as np
O=N/'fluent-simplec-capability-v002'

def audit_capability():
 verify();status=json.loads((O/'status.json').read_text());assert status['stage']=='closed' and status['iteration']==1500 and status['new_iterations']==0
 for x,y in [('wing-final.cas.h5','source.cas.h5'),('wing-final.dat.h5','source.dat.h5'),('final-settings.json','expected-source-settings.json')]:assert sha(F/x)==sha(O/y)
 loaded=json.loads((O/'loaded-source-settings.json').read_text());source=json.loads((O/'expected-source-settings.json').read_text());compare_fixed(loaded,{k:source[k] for k in state_keys()});candidate=json.loads((O/'candidate-settings.json').read_text());validate_state(loaded,candidate);compare_fixed(candidate,json.loads((O/'final-settings.json').read_text()))
 assert json.loads((O/'native-state-differences.json').read_text())==diff_inventory(loaded,candidate)
 for name in ['loaded-targets.json','candidate-targets.json']:assert json.loads((O/name).read_text())=={k:1e-5 for k in EQ}
 for name in ['loaded','candidate']:
  assert float(json.loads((O/(name+'-live-iteration-raw.json')).read_text())['query_value'])==1500
  same_native(O/(name+'-restart.cas.h5') if name=='loaded' else O/'candidate.cas.h5',O/(name+'-restart.dat.h5') if name=='loaded' else O/'candidate.dat.h5',O/'source.cas.h5',O/'source.dat.h5')
 pi,pr=raw_residual(O/'source.dat.h5');ci,cr=raw_residual(O/'candidate.dat.h5');assert np.array_equal(pi,ci) and ci[-1]==1500;same_residual_prefix(pi,pr,ci,cr,1500)
 proof=json.loads((O/'zero-advance-proof.json').read_text());assert proof['new_iterations']==0 and proof['live_loaded_iteration']==proof['live_candidate_iteration']==proof['native_current_iteration']==1500 and proof['energy_disabled_after']
 q=json.loads((O/'geometry-quality.json').read_text());assert q['minimum_orthogonal_quality']>=.10 and q['minimum_cell_volume_m3']>0
 assert_zero_advance_source((W/'run_wing_simplec_capability_v002.py').read_text(encoding='utf-8-sig'));version=json.loads((O/'solver-version.json').read_text());assert version['processor_count']==2 and version['connected_before_any_advance']
 names=['source.cas.h5','source.dat.h5','expected-source-settings.json','loaded-restart.cas.h5','loaded-restart.dat.h5','candidate.cas.h5','candidate.dat.h5','loaded-source-settings.json','after-switch-defaults-settings.json','candidate-settings.json','final-settings.json','native-state-differences.json','loaded-targets.json','candidate-targets.json','loaded-live-iteration-raw.json','candidate-live-iteration-raw.json','loaded-field-proof.json','candidate-field-proof.json','zero-advance-proof.json','geometry-quality.json','allowed-flow-schemes.json','active-under-relaxation-names.json','active-AMG-names.json','solver-version.json','solver-transcript.txt','status.json']
 result=dict(zero_advance_capability_pass=True,new_iterations=0,source_and_candidate_iteration=1500,initial_all_seven_fields_and_topology_exact=True,raw_residual_history_exact=True,original_stored_criteria_unchanged=True,energy_disabled=True,declared_SIMPLEC_controls_verified=True,physical_spatial_reference_flux_projection_exact=True,all_native_changes=diff_inventory(loaded,candidate),candidate_native_state=candidate,original_source_still_underconverged=True,no_numerical_qualification_yet=True,no_new_residual_map=True,source_sha256={str(O/n):sha(O/n) for n in names})
 verify();dump(O/'capability-audit-v001.json',result);return result
if __name__=='__main__':print(json.dumps(audit_capability(),indent=2))
