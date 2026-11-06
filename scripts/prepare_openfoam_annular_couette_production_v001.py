"""Prepare exactly 395 additional annular Couette steps only after a selected accepted five-step pilot; no CFD."""
from pathlib import Path
import argparse,json,re
import prepare_openfoam_annular_couette_v003 as pilot
from prepare_openfoam_annular_couette_v003 import *
PILOT_PREPARATION_SHA='576583b56881b4b9b01dc141b9acbd2dd0e86ab66cc9532650fa8b056705c7d8'
PILOT_ANALYSIS_SHA='fb790d9e753b334f32b2673009156b902c6b71f514fdb6c49e24b78d5bc013e7'
DEPENDENCIES={'prepare_openfoam_annular_couette_v003.py':PILOT_PREPARATION_SHA,'prepare_openfoam_annular_couette_v002.py':ANCESTOR_PREP_SHA,'analyse_annular_couette_openfoam_v002.py':PILOT_ANALYSIS_SHA}

def checked_case(path):
 p=Path(path).resolve();require(p==(run_root()/'openfoam-production-v001').resolve(),'Dedicated 395-step annular continuation only');return p

def canonical_project_path(raw):
 text=str(raw).replace('\\','/');suffix=None
 for prefix in ('E:/CFD-Learning-Lab','/mnt/e/CFD-Learning-Lab'):
  if text==prefix:suffix='';break
  if text.startswith(prefix+'/'):suffix=text[len(prefix)+1:];break
 require(suffix is not None and ':' not in suffix and chr(0) not in suffix and (suffix=='' or all(p not in ('','.','..') for p in suffix.split('/'))),'Exact unambiguous project-native path required')
 base=root().resolve();result=(base/suffix).resolve();require(result==base or base in result.parents,'Project-native path containment');return result

def verify_dependencies():
 for name,h in DEPENDENCIES.items():require(sha(Path(__file__).with_name(name))==h,'Frozen selected dependency changed '+name)

def dictionaries():
 d=pilot.dictionaries();s=d['system/controlDict'];require(s.count('startFrom startTime;')==1 and s.count('endTime 0.25;')==1,'Exact five-step source control');d['system/controlDict']=s.replace('startFrom startTime;','startFrom latestTime;').replace('endTime 0.25;','endTime 20;');return d

def verify_inputs(case):
 case=checked_case(case);verify_dependencies();m=json.loads((case/'provenance/prepared-input-sha256.json').read_text());require(all(sha(case/n)==h for n,h in m.items()),'Prepared continuation inputs changed')
 for n,t in dictionaries().items():
  if n not in ('0/U','0/p'):require((case/n).read_text()==t,'Unchanged approved continuation physics/control '+n)
 receipt=json.loads((case/'provenance/inherited-accepted-pilot.json').read_text())
 for r in receipt['files']:require(sha(canonical_project_path(r['source']))==r['sha256']==sha(case/r['destination']),'Accepted pilot native/source byte copy changed')
 return m

def prepare(case,audit,checks_sha,manifest_sha,fields_sha):
 case=checked_case(case);audit=Path(audit).resolve();source=run_root()/'openfoam-pilot-v003';require(audit==(private()/'openfoam-pilot-checks-v002').resolve() and not case.exists(),'Exact selected accepted audit and fresh continuation');verify_dependencies();require(all(re.fullmatch('[a-f0-9]{64}',h) for h in (checks_sha,manifest_sha,fields_sha)),'Actual selected audit SHA arguments required')
 require(sha(audit/'checks.json')==checks_sha and sha(audit/'source-manifest.json')==manifest_sha and sha(audit/'fields.npz')==fields_sha,'Actual selected accepted pilot outputs changed');c=json.loads((audit/'checks.json').read_text());require(c['all_checks_passed'] is True and not c['failed_checks'] and all(v is True for v in c['checks'].values()) and len(c['checks'])==20 and c['actual_steps']==5 and c['native_states']==6 and c['outer_corrector_budget']==30 and c['analysis_source_sha256']==PILOT_ANALYSIS_SHA and c['preparation_source_sha256']==PILOT_PREPARATION_SHA and c['audited_fields_npz_sha256']==fields_sha and canonical_project_path(c['native_case'])==source.resolve(),'Only genuinely accepted 30-outer five-step pilot')
 m=json.loads((audit/'source-manifest.json').read_text());require(m['inputs_unchanged_during_analysis'] is True and all(sha(canonical_project_path(r['path']))==r['sha256_after'] for r in m['sources']),'Every selected original native input still unchanged');pilot.verify_inputs(source);status=json.loads((source/'run-status.json').read_text());require(status['stage']=='finished_pilot_execution' and status['actual_steps']==5 and status['physical_time_s']==.25,'Actual completed pilot status');require(sha(private()/AMENDMENT_NAME)==AMENDMENT_SHA and sha(private()/'TARGETS.md')==TARGET_SHA,'Frozen original targets and 30-outer amendment')
 case.mkdir(parents=True);rows=[]
 def one(a,b):
  h=sha(a);copied(a,b);rows.append(dict(source=str(a),destination=b.relative_to(case).as_posix(),sha256=h))
 def tree(a,b):
  require(a.is_dir(),'Actual inherited native folder missing '+str(a))
  for p in sorted(a.rglob('*')):
   if p.is_file():one(p,b/p.relative_to(a))
 for folder in ('0','constant'):tree(source/folder,case/folder)
 for rank in (0,1):
  for i in range(6):
   label=f'{i*.05:.10g}';tree(source/f'processor{rank}'/label,case/f'processor{rank}'/label)
  raw=(case/f'processor{rank}/0.25/uniform/time').read_text();require(all(abs(float(re.search(r'\b'+n+r'\s+('+FLOAT+r')\s*;',raw)[1])-v)<=1e-12 for n,v in [('value',.25),('deltaT',.05),('index',5)]),'Actual inherited latest native time/index5')
 for p in sorted((source/'system').rglob('*')):
  if p.is_file() and p.relative_to(source).as_posix()!='system/controlDict':one(p,case/p.relative_to(source))
 write_new(case/'system/controlDict',dictionaries()['system/controlDict'])
 tree(source/'provenance',case/'provenance/accepted-pilot-v003')
 for name in ('TARGETS.md',AMENDMENT_NAME,'setup.json','geometry.npz','source-geometry.json','official-research-source-manifest.json','native-mapping.json','preflight-validation.json','preflight-readonly-recovery.json','native-registration-process-status.json','installed-source-sha256.txt','native-initial-copy-receipt.json','validate_annular_openfoam_preflight_v001.py'):one(source/'provenance'/name,case/'provenance'/name)
 tree(source/'provenance/installed-sources',case/'provenance/installed-sources')
 for name in ('mesh-check.log','initial-native-registration.log','parallel-native-registration.log'):one(source/name,case/name)
 for name in ('solver.log','run-status.json','declared-comparison.json'):one(source/name,case/'provenance/accepted-pilot-v003'/name)
 for name in ('checks.json','fields.npz','source-manifest.json','step-history.csv','residual-history.csv'):one(audit/name,case/'provenance/accepted-pilot-audit'/name)
 for p in (Path(__file__),private()/'run-openfoam-production-v001.sh',private()/'Launch-OpenFOAM-Production-v001.ps1',*[Path(__file__).with_name(n) for n in DEPENDENCIES]):one(p,case/'provenance'/p.name)
 for r in rows:require(sha(canonical_project_path(r['source']))==r['sha256']==sha(case/r['destination']),'Every actual accepted pilot/native source copied byte-identically')
 before=pilot.dictionaries();after=dictionaries();require([n for n in after if after[n]!=before[n]]==['system/controlDict'] and after['system/controlDict'].replace('startFrom latestTime;','startFrom startTime;').replace('endTime 20;','endTime 0.25;')==before['system/controlDict'],'Only startFrom/endTime controller changes')
 selection=dict(pilot_native_case=str(source),pilot_audit_path=str(audit),pilot_checks_sha256=checks_sha,pilot_manifest_sha256=manifest_sha,pilot_fields_sha256=fields_sha,pilot_preparation_sha256=PILOT_PREPARATION_SHA,pilot_analysis_sha256=PILOT_ANALYSIS_SHA);json_new(case/'provenance/accepted-pilot-selection.json',selection)
 receipt=dict(accepted_pilot=selection,files=rows,inherited_actual_physical_steps=5,inherited_native_states=6,original_pilot_authoritative_for_steps_0_to_5=True,source_restart_time_s=.25,source_restart_index=5,required_first_new_time_s=.30,required_first_new_index=6,additional_steps=395,total_steps=400,final_time_s=20,deltaT_s=.05,outer_correctors=30,pressure_correctors=2,controller_changes_only=['startFrom startTime -> latestTime','endTime 0.25 -> 20'],reinitialization=False,accepted_pilot_replayed=False,serial_reconstruction_requested=False,production_admitted=False);json_new(case/'provenance/inherited-accepted-pilot.json',receipt)
 json_new(case/'declared-comparison.json',dict(phase='bounded_395_step_annular_continuation',**receipt,canonical_TARGETS_sha256=TARGET_SHA,outer_budget_amendment_sha256=AMENDMENT_SHA,physical_validation=False))
 json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file()});json_new(case/'run-status.json',dict(stage='prepared_accepted_pilot_continuation_only',solver_started=False,inherited_actual_steps=5,additional_steps_authorized=395,total_steps_authorized=400,physical_time_s=.25,production_admitted=False));verify_inputs(case)
 return dict(case=str(case),prepared=True,solver_started=False,copied_files=len(rows),inherited_actual_steps=5,additional_steps=395,total_steps=400,prepared_input_manifest_sha256=sha(case/'provenance/prepared-input-sha256.json'),inherited_accepted_pilot_receipt_sha256=sha(case/'provenance/inherited-accepted-pilot.json'),accepted_pilot_selection_sha256=sha(case/'provenance/accepted-pilot-selection.json'))

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,required=True);ap.add_argument('--verify-inputs',action='store_true');ap.add_argument('--pilot-audit',type=Path);ap.add_argument('--expected-pilot-checks-sha256');ap.add_argument('--expected-pilot-manifest-sha256');ap.add_argument('--expected-pilot-fields-sha256');a=ap.parse_args()
 if a.verify_inputs:result=dict(verified_files=len(verify_inputs(a.case)))
 else:
  require(all((a.pilot_audit,a.expected_pilot_checks_sha256,a.expected_pilot_manifest_sha256,a.expected_pilot_fields_sha256)),'Actual accepted pilot and all3selected SHA arguments required');result=prepare(a.case,a.pilot_audit,a.expected_pilot_checks_sha256,a.expected_pilot_manifest_sha256,a.expected_pilot_fields_sha256)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()