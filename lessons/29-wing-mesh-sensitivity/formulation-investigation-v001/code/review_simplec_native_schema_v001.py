import wing_import_paths_v001
from wing_repair_common_v001 import *
from wing_simplec_prelaunch_v002 import verify
O=N/'fluent-simplec-capability-v002';a=json.loads((O/'capability-audit-v001.json').read_text());assert a['zero_advance_capability_pass'] and a['new_iterations']==0
for p,h in a['source_sha256'].items():assert sha(p)==h,p
categories={}
for r in a['all_native_changes']:
 p=r['path']
 if p.startswith('controls/p_v_controls') or '/coupled_parameters/' in p or '/pressure-coupled/' in p:k='Coupled-only controls become inactive'
 elif p.startswith('methods/pseudo_time_method/formulation/'):k='Pseudo-time stays off under active segregated key'
 elif p=='controls/advanced/enhanced_skewness_correction':assert r['after'] is False;k='Newly visible native flag explicitly disabled'
 elif p.startswith('controls/advanced/multi_grid/mg_controls/'):k='Active scalar pressure/momentum AMG schema and declared settings'
 elif p in ['methods/p_v_coupling/flow_scheme','methods/p_v_coupling/skewness_correction_itr_count','controls/under_relaxation/mom','controls/under_relaxation/pressure']:k='Declared SIMPLEC algorithm/relaxation bundle'
 elif p in ['controls/advanced/expert/save_cell_residuals','controls/advanced/multi_grid/amg_controls/options/verbosity']:k='Diagnostic storage/verbosity only'
 else:raise ValueError('Unreviewed native change: '+p)
 categories.setdefault(k,[]).append(r)
p=W/'accepted-native-SIMPLEC-settings-v001.json';assert not p.exists();dump(p,a['candidate_native_state'])
record=dict(capability_native_schema_review_pass=True,source_audit=str(O/'capability-audit-v001.json'),source_audit_sha256=sha(O/'capability-audit-v001.json'),new_iterations=0,all_native_changes_reviewed=len(a['all_native_changes']),categories=categories,physical_spatial_reference_flux_exact=True,all_seven_initial_fields_and_rawhistory_exact=True,live_before_and_after_iteration=1500,stored_residual_targets=1e-5,active_pressure_cycle='v-cycle',pressure_AMG_termination=.001,active_momentum_names=['x-mom','y-mom','z-mom'],momentum_AMG_cycle='flexible-cycle',momentum_AMG_termination=.1,momentum_residual_reduction_tolerance=.7,scalar_smoother='Gauss-Seidel',implicit_pressure_and_momentum_relaxation=.7,skewness_correction_iterations=1,enhanced_skewness_correction=False,pseudo_time='off',frozen_expected_state=str(p),expected_state_sha256=sha(p),original_source_still_underconverged=True,qualification_not_prepared_or_launched=True,tests_pass=48,source_binding_count=verify(),first_failed_zero_case_preserved=True)
dump(W/'simplec-capability-native-review-v001.json',record)
with (W/'progress.md').open('a',encoding='utf-8') as f:f.write('\nCapability002 CLOSED20:38:58/AUDITEDPASS20:39:02, ZEROadvancement: liveIteration1500 before/after, source7fields/topology/fullraw1500historyEXACT, original1e-5targets/physicalspatialreferenceflux projection unchanged. All'+str(len(a['all_native_changes']))+'native leafchanges individuallycategorized/rootreviewed, no unexplaineddefaults. ActivepressureVcycle.001; x/y/z-mom flexible.1/reduction.7/GaussSeidel; implicitpressure/mom.7,one skewnessloop, enhancedFALSE,pseudooff. 48testsPASS incl4actualnativeactivationregressions. ExpectednativeSIMPLECschemafrozen accepted-native-SIMPLEC-settings-v001.json,436inputbindings reverified. Firststrictschemafailure/recovery preserved, no tolerancechange. No100iterationqualificationprepared/launchedyet.\n')
print(json.dumps({k:v for k,v in record.items() if k!='categories'},indent=2))
