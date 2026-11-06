"""Prepare exactly 795 additional tandem-wake steps after a source-bound accepted five-step pilot; no CFD."""
from pathlib import Path
import argparse,json,re
import prepare_openfoam_tandem_offset_wake_v002 as pilot
import tandem_offset_openfoam_common_v002 as c
PILOT_PREPARATION_SHA='8c92ec571ec735e886fe7e464425654fb92f5d16fd7e649d1e7890fefd8eac48'
PILOT_ANALYSIS_SHA='8875723bd4df1af8f99ebc0573712b66db4fcd400238972acef3806267cca7a8'
COMMON_SHA='d0f0a34cc0bc50533f9072fee1d4329156bd52da5a178318003d857ff8fd98e9'
DEPENDENCIES={'prepare_openfoam_tandem_offset_wake_v002.py':PILOT_PREPARATION_SHA,'analyse_tandem_offset_wake_openfoam_v002.py':PILOT_ANALYSIS_SHA,'tandem_offset_openfoam_common_v002.py':COMMON_SHA}
def checked_case(path):
 p=Path(path).resolve();c.require(p==(pilot.run_root()/'openfoam-production-v001').resolve(),'Dedicated bounded tandem production only');return p
def canonical_project_path(raw):
 text=str(raw).replace('\\','/');suffix=None
 for prefix in ('E:/CFD-Learning-Lab','/mnt/e/CFD-Learning-Lab'):
  if text==prefix:suffix='';break
  if text.startswith(prefix+'/'):suffix=text[len(prefix)+1:];break
 c.require(suffix is not None and ':' not in suffix and chr(0) not in suffix and (suffix=='' or all(v not in ('','.','..') for v in suffix.split('/'))),'Exact contained project path required');base=pilot.root().resolve();p=(base/suffix).resolve();c.require(p==base or base in p.parents,'Resolved project containment');return p
def verify_dependencies():
 for n,h in DEPENDENCIES.items():c.require(c.sha(Path(__file__).with_name(n))==h,'Frozen dependency changed '+n)
def dictionaries():
 d=pilot.dictionaries();s=d['system/controlDict'];c.require(s.count('startFrom startTime;')==s.count('endTime 0.25;')==1,'Exact selected pilot controller');d['system/controlDict']=s.replace('startFrom startTime;','startFrom latestTime;').replace('endTime 0.25;','endTime 40;');return d
def is_derived_native(rel):
 p=Path(rel);return p.name in ('vorticity','wallShearStress') and (p.parts[0]=='0' or p.parts[0] in ('processor0','processor1'))
def verify_inputs(case,allow_derived_changes=False):
 case=checked_case(case);verify_dependencies();m=json.loads((case/'provenance/prepared-primary-sha256.json').read_text())
 for n,h in m.items():c.require(c.sha(case/n)==h,'Primary prepared input changed '+n)
 r=json.loads((case/'provenance/inherited-accepted-pilot.json').read_text())
 for row in r['files']:
  c.require(c.sha(canonical_project_path(row['source']))==row['sha256'],'Original accepted pilot/source changed '+row['source'])
  if row['classification']=='primary' or not allow_derived_changes:c.require(c.sha(case/row['destination'])==row['sha256'],'Inherited input changed '+row['destination'])
 for n,t in dictionaries().items():
  if n not in ('0/U','0/p'):c.require((case/n).read_text()==t,'Unchanged approved physics/controller '+n)
 return m

def prepare(case,audit,checks_sha,manifest_sha,fields_sha):
 case=checked_case(case);audit=Path(audit).resolve();source=pilot.checked_case(pilot.run_root()/'openfoam-pilot-v002');c.require(audit==(pilot.private()/'openfoam-pilot-checks-v001').resolve() and not case.exists(),'Actual selected pilot audit and fresh production');verify_dependencies()
 c.require(all(isinstance(h,str) and re.fullmatch('[a-f0-9]{64}',h) for h in (checks_sha,manifest_sha,fields_sha)),'Three actual selected SHA arguments required')
 for n,h in [('checks.json',checks_sha),('source-manifest.json',manifest_sha),('fields.npz',fields_sha)]:c.require(c.sha(audit/n)==h,'Actual accepted audit output changed '+n)
 a=json.loads((audit/'checks.json').read_text());c.require(a['all_checks_passed'] is True and not a['failed_checks'] and a['checks'] and all(v is True for v in a['checks'].values()) and a['actual_steps']==5 and a['native_states']==6 and a['analysis_source_sha256']==PILOT_ANALYSIS_SHA and a['preparation_source_sha256']==PILOT_PREPARATION_SHA and a['common_source_sha256']==COMMON_SHA and a['audited_fields_npz_sha256']==fields_sha and canonical_project_path(a['native_case'])==source,'Only genuinely accepted selected tandem pilot')
 m=json.loads((audit/'source-manifest.json').read_text());c.require(m['inputs_unchanged_during_analysis'] is True,'Accepted unchanged audit inputs')
 for row in m['sources']:c.require(c.sha(canonical_project_path(row['path']))==row['sha256_after'],'Original accepted native input changed')
 pilot.verify_inputs(source);s=json.loads((source/'run-status.json').read_text());c.require(s['stage']=='finished_pilot_execution' and s['actual_steps']==5 and s['physical_time_s']==.25,'Actual completed five-step pilot');c.require(c.sha(pilot.private()/'TARGETS.md')==pilot.TARGET_SHA,'Frozen unchanged targets')
 case.mkdir(parents=True);rows=[];destinations={}
 def one(a,b):
  h=c.sha(a);rel=b.relative_to(case).as_posix()
  if rel in destinations:c.require(destinations[rel]==h and c.sha(b)==h,'Exact duplicate-copy identity');return
  c.copied(a,b);destinations[rel]=h;derived=is_derived_native(rel);row=dict(source=str(a),destination=rel,sha256=h,classification='derived_native_output' if derived else 'primary')
  if derived:
   baseline='provenance/derived-restart-baselines/'+rel;c.copied(a,case/baseline);row['baseline_path']=baseline
  rows.append(row)
 def tree(a,b):
  c.require(a.is_dir(),'Actual inherited folder absent '+str(a))
  for x in sorted(a.rglob('*')):
   if x.is_file():one(x,b/x.relative_to(a))
 for folder in ('0','constant'):tree(source/folder,case/folder)
 for rank in (0,1):
  for i in range(6):
   label=f'{i*.05:.10g}';tree(source/f'processor{rank}'/label,case/f'processor{rank}'/label)
  c.time_metadata(c.Evidence(),case/f'processor{rank}/0.25',5)
 for x in sorted((source/'system').rglob('*')):
  if x.is_file() and x.name!='controlDict':one(x,case/x.relative_to(source))
 c.write_new(case/'system/controlDict',dictionaries()['system/controlDict']);tree(source/'provenance',case/'provenance/accepted-pilot-v002')
 for n in ('TARGETS.md','EXECUTION-PLAN-v001.md','setup.json','geometry.npz','source-geometry.json','official-primary-source-urls.json','native-mapping.json','preflight-validation.json','native-registration-process-status.json','installed-source-sha256.txt','installed-source-list.json'):one(source/'provenance'/n,case/'provenance'/n)
 tree(source/'provenance/installed-sources',case/'provenance/installed-sources')
 for n in ('mesh-check.log','initial-native-registration.log','parallel-native-registration.log'):one(source/n,case/n)
 for n in ('solver.log','run-status.json','declared-comparison.json'):one(source/n,case/'provenance/accepted-pilot-v002'/n)
 for n in ('checks.json','fields.npz','source-manifest.json','step-history.csv','residual-history.csv','force-history.csv'):one(audit/n,case/'provenance/accepted-pilot-audit'/n)
 for x in (Path(__file__),pilot.private()/'run-openfoam-production-v001.sh',pilot.private()/'Launch-OpenFOAM-Production-v001.ps1',*[Path(__file__).with_name(n) for n in DEPENDENCIES]):one(x,case/'provenance'/x.name)
 for row in rows:c.require(c.sha(canonical_project_path(row['source']))==row['sha256']==c.sha(case/row['destination']),'Actual inherited byte equality')
 before=pilot.dictionaries();after=dictionaries();c.require([n for n in after if after[n]!=before[n]]==['system/controlDict'] and after['system/controlDict'].replace('startFrom latestTime;','startFrom startTime;').replace('endTime 40;','endTime 0.25;')==before['system/controlDict'],'Only two approved controller changes')
 selection=dict(pilot_native_case=str(source),pilot_audit_path=str(audit),pilot_checks_sha256=checks_sha,pilot_manifest_sha256=manifest_sha,pilot_fields_sha256=fields_sha,pilot_preparation_sha256=PILOT_PREPARATION_SHA,pilot_analysis_sha256=PILOT_ANALYSIS_SHA,common_sha256=COMMON_SHA);c.json_new(case/'provenance/accepted-pilot-selection.json',selection)
 receipt=dict(accepted_pilot=selection,files=rows,inherited_actual_physical_steps=5,inherited_native_states=6,original_pilot_authoritative_for_steps_0_to_5=True,source_restart_time_s=.25,source_restart_index=5,required_first_new_time_s=.30,required_first_new_index=6,additional_steps=795,total_steps=800,final_time_s=40,deltaT_s=.05,outer_correctors=30,pressure_correctors=2,controller_changes_only=['startFrom startTime -> latestTime','endTime 0.25 -> 40'],reinitialization=False,accepted_pilot_replayed=False,serial_reconstruction_requested=False,derived_restart_outputs_classified_before_CFD=True,production_admitted=False);c.json_new(case/'provenance/inherited-accepted-pilot.json',receipt)
 c.json_new(case/'declared-comparison.json',dict(phase='bounded_795_step_tandem_continuation',**receipt,target_sha256=pilot.TARGET_SHA,physical_validation=False))
 c.json_new(case/'provenance/prepared-primary-sha256.json',{x.relative_to(case).as_posix():c.sha(x) for x in sorted(case.rglob('*')) if x.is_file() and not is_derived_native(x.relative_to(case).as_posix())});c.json_new(case/'run-status.json',dict(stage='prepared_accepted_pilot_continuation_only',solver_started=False,inherited_actual_steps=5,additional_steps_authorized=795,total_steps_authorized=800,physical_time_s=.25,production_admitted=False));verify_inputs(case)
 return dict(case=str(case),prepared=True,solver_started=False,copied_files=len(rows),inherited_actual_steps=5,additional_steps=795,total_steps=800,prepared_primary_manifest_sha256=c.sha(case/'provenance/prepared-primary-sha256.json'),inherited_accepted_pilot_receipt_sha256=c.sha(case/'provenance/inherited-accepted-pilot.json'),accepted_pilot_selection_sha256=c.sha(case/'provenance/accepted-pilot-selection.json'))

def finalize_native(case):
 case=checked_case(case);verify_inputs(case,True);raw=(case/'solver.log').read_text();times=[float(t) for t in re.findall(r'^Time = ('+c.FLOAT+r')s?\s*$',raw,re.M)];c.require(len(times)==795 and all(abs(t-i*.05)<=1e-10 for i,t in enumerate(times,6)) and re.search(r'^End\s*$',raw,re.M) and 'Finalising parallel run' in raw and 'FOAM FATAL' not in raw,'Actual completed 795 new native steps')
 for rank in (0,1):
  for i in range(1,801):
   folder=case/f'processor{rank}'/f'{i*.05:.10g}';c.require(all((folder/n).is_file() for n in ('U','p','phi','vorticity','wallShearStress')),'Literal all-rank fields');c.time_metadata(c.Evidence(),folder,i)
 rows=[]
 for row in json.loads((case/'provenance/inherited-accepted-pilot.json').read_text())['files']:
  if row['classification']!='derived_native_output':continue
  p=case/row['destination'];b=case/row['baseline_path'];c.require(c.sha(b)==row['sha256'],'Preserved derived baseline');n=9472 if row['destination'].startswith('0/') else c.counted(case/Path(row['destination']).parts[0]/'0/polyMesh/cellProcAddressing')[0];dims=[0,0,-1,0,0,0,0] if p.name=='vorticity' else [0,2,-2,0,0,0,0];x=c.field(b.read_text(),n,3,dims);y=c.field(p.read_text(),n,3,dims);delta=max((abs(v-w) for a,z in zip(x,y) for v,w in zip(a,z)),default=0.);rows.append(dict(path=row['destination'],original_sha256=row['sha256'],current_sha256=c.sha(p),byte_preserved=c.sha(p)==row['sha256'],native_internal_max_abs_difference=delta))
 result=dict(actual_solver_completed=True,log_End=True,additional_actual_steps=795,total_actual_steps=800,native_rank_states=801,rank_time_records=1600,first_new_time_s=.30,first_new_index=6,final_time_s=40,primary_prepared_and_original_pilot_inputs_unchanged=True,original_pilot_authoritative_for_steps_0_to_5=True,strict_all_derived_copies_preserved=all(r['byte_preserved'] for r in rows),derived_outputs=rows,serial_reconstruction_requested=False,CFD_replayed=False);c.json_new(case/'provenance/native-production-completion.json',result);return result

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,required=True);g=ap.add_mutually_exclusive_group();g.add_argument('--verify-inputs',action='store_true');g.add_argument('--finalize-native',action='store_true');ap.add_argument('--pilot-audit',type=Path);ap.add_argument('--expected-pilot-checks-sha256');ap.add_argument('--expected-pilot-manifest-sha256');ap.add_argument('--expected-pilot-fields-sha256');a=ap.parse_args()
 if a.verify_inputs:r=dict(verified_files=len(verify_inputs(a.case)))
 elif a.finalize_native:r=finalize_native(a.case)
 else:
  c.require(all((a.pilot_audit,a.expected_pilot_checks_sha256,a.expected_pilot_manifest_sha256,a.expected_pilot_fields_sha256)),'Actual selected audit and three hashes required');r=prepare(a.case,a.pilot_audit,a.expected_pilot_checks_sha256,a.expected_pilot_manifest_sha256,a.expected_pilot_fields_sha256)
 print(json.dumps(r,indent=2))
if __name__=='__main__':main()