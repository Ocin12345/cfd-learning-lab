"""Prepare the amended 30-outer Couette pilot by byte-copying proven native rest inputs; no native commands or CFD."""
from pathlib import Path
import argparse,json,re
import prepare_openfoam_annular_couette_v002 as ancestor
from prepare_openfoam_annular_couette_v002 import *
ANCESTOR_PREP_SHA='9e07e2982beca3a34242478bb23bf98e7d8ccb14b00836e69831ffc713bdfe34'
ANCESTOR_CHECKS_SHA='5a65e0a88fff60890ec43c477240fc935410143ff374e40bc561dbd075574a66'
ANCESTOR_AUDIT_MANIFEST_SHA='1a6fcec03181f6fb6eb893dbb019ef9f7e7c870df8490df5a3f502fc5e3cced7'
AMENDMENT_NAME='OPENFOAM-PRESSURE-COUPLING-AMENDMENT-v001.md'
AMENDMENT_SHA='157b314383f9f70b780445bc773a07f0348a3b659fcbd1be54629965447f874c'
OUTER_CORRECTORS=30

def checked_case(p):
 p=Path(p).resolve();require(p==(run_root()/'openfoam-pilot-v003').resolve(),'Dedicated amended annular pilot only');return p

def dictionaries():
 d=ancestor.dictionaries();old=d['system/fvSolution'];require(old.count('nOuterCorrectors 10;')==1,'Exact original outer budget');d['system/fvSolution']=old.replace('nOuterCorrectors 10;','nOuterCorrectors 30;');return d

def verify_inputs(c):
 c=checked_case(c);require(sha(Path(ancestor.__file__))==ANCESTOR_PREP_SHA,'Frozen ancestral geometry/source helper');require(sha(c/'provenance/TARGETS.md')==TARGET_SHA and sha(c/'provenance'/AMENDMENT_NAME)==AMENDMENT_SHA,'Unchanged targets and explicit budget amendment');m=json.loads((c/'provenance/prepared-input-sha256.json').read_text());require(all(sha(c/n)==h for n,h in m.items()),'Prepared native-initial copy changed');require((c/'system/fvSolution').read_text()==dictionaries()['system/fvSolution'],'Only approved 30-outer coupling controls');return m

def prepare(source,case):
 case=checked_case(case);source=Path(source).resolve();require(source==(run_root()/'openfoam-pilot-v002').resolve() and not case.exists(),'Exact preserved ancestor and fresh case');require(sha(Path(ancestor.__file__))==ANCESTOR_PREP_SHA,'Frozen ancestral source');ancestor.verify_inputs(source)
 amendment=private()/AMENDMENT_NAME;require(sha(amendment)==AMENDMENT_SHA and sha(private()/'TARGETS.md')==TARGET_SHA,'Predeclared budget-only amendment')
 audit=private()/'openfoam-pilot-checks-v001';require(sha(audit/'checks.json')==ANCESTOR_CHECKS_SHA and sha(audit/'source-manifest.json')==ANCESTOR_AUDIT_MANIFEST_SHA,'Preserved actual 10-outer numeric miss');old=json.loads((audit/'checks.json').read_text());require(old['all_checks_passed'] is False and old['failed_checks']==['every_step_last_outer_initial_p_U_residuals_le1e_6'] and len(old['checks'])==19 and sum(old['checks'].values())==18,'Preserved 18/19 original pilot');am=json.loads((audit/'source-manifest.json').read_text());require(am['inputs_unchanged_during_analysis'] is True and all(sha(Path(r['path']))==r['sha256_after'] for r in am['sources']),'All original audited input hashes still match')
 case.mkdir(parents=True);rows=[]
 def copy_one(a,b):
  h=sha(a);copied(a,b);rows.append(dict(source=str(a),destination=b.relative_to(case).as_posix(),sha256=h));return h
 def tree(a,b):
  require(a.is_dir(),'Missing inherited native directory '+str(a))
  for p in sorted(a.rglob('*')):
   if p.is_file():copy_one(p,b/p.relative_to(a))
 for folder in ('0','constant','processor0/0','processor1/0'):tree(source/folder,case/folder)
 for p in sorted((source/'system').rglob('*')):
  if p.is_file() and p.relative_to(source).as_posix()!='system/fvSolution':copy_one(p,case/p.relative_to(source))
 write_new(case/'system/fvSolution',dictionaries()['system/fvSolution'])
 tree(source/'provenance',case/'provenance/ancestor-v002')
 for name in ('TARGETS.md','setup.json','geometry.npz','source-geometry.json','official-research-source-manifest.json','native-mapping.json','preflight-validation.json','preflight-readonly-recovery.json','native-registration-process-status.json','installed-source-sha256.txt','validate_annular_openfoam_preflight_v001.py'):
  copy_one(source/'provenance'/name,case/'provenance'/name)
 tree(source/'provenance/installed-sources',case/'provenance/installed-sources')
 for name in ('mesh-check.log','initial-native-registration.log','parallel-native-registration.log'):copy_one(source/name,case/name)
 copy_one(source/'input-annular-couette.msh',case/'input-annular-couette.msh');copy_one(source/'run-status.json',case/'provenance/ancestor-v002/run-status.json');copy_one(source/'solver.log',case/'provenance/ancestor-v002/solver.log');copy_one(source/'declared-comparison.json',case/'provenance/ancestor-v002/declared-comparison.json')
 for name in ('checks.json','source-manifest.json','step-history.csv','residual-history.csv'):copy_one(audit/name,case/'provenance/ancestor-v002/audit'/name)
 for p in (amendment,Path(__file__),private()/'run-openfoam-v003.sh',private()/'Launch-OpenFOAM-v003.ps1'):copy_one(p,case/'provenance'/p.name)
 olddecl=json.loads((source/'declared-comparison.json').read_text());olddecl.update(outer_correctors=30,outer_budget_amendment_sha256=AMENDMENT_SHA,ancestral_native_case=str(source),ancestral_10_outer_pilot_passed=False,ancestral_10_outer_audit_checks_sha256=ANCESTOR_CHECKS_SHA,initial_native_state='byte-identical original actual rest time0, no step5 restart',native_conversion_or_registration_replayed=False);json_new(case/'declared-comparison.json',olddecl)
 for r in rows:require(sha(Path(r['source']))==r['sha256']==sha(case/r['destination']),'Every inherited native/source byte copy')
 for prefix in ('','processor0/','processor1/'):
  require(sorted(p.name for p in (case/prefix).iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name))==['0'],'Only original rest0 inherited, no solved ancestor times')
 require(max(abs(v) for r in field(case/'0/U',CELLS,3,[0,1,-1,0,0,0,0]) for v in r)<=1e-12 and max(abs(v) for v in field(case/'0/p',CELLS,1,[0,2,-2,0,0,0,0]))<=1e-12,'Copied actual cell rest0/p0')
 original=ancestor.dictionaries();current=dictionaries();require([n for n in current if current[n]!=original[n]]==['system/fvSolution'] and current['system/fvSolution'].replace('nOuterCorrectors 30;','nOuterCorrectors 10;')==original['system/fvSolution'],'Budget-only physics/control difference')
 receipt=dict(prepared_native_initial_copy=True,ancestor_native_case=str(source),ancestor_audit_checks_sha256=ANCESTOR_CHECKS_SHA,ancestor_audit_manifest_sha256=ANCESTOR_AUDIT_MANIFEST_SHA,ancestor_pilot_checks_passed=18,ancestor_pilot_checks_total=19,ancestor_numeric_miss_preserved=True,ancestor_solved_steps_used_for_restart=0,initial_native_time_s=0,initial_native_fields_byte_identical=True,native_mesh_and_rank_addressing_byte_identical=True,copied_files=rows,native_conversion_replayed=False,native_registration_replayed=False,solver_started=False,outer_correctors=30,pressure_correctors=2,amendment_sha256=AMENDMENT_SHA,canonical_TARGETS_sha256=TARGET_SHA,physical_or_threshold_change=False)
 json_new(case/'provenance/native-initial-copy-receipt.json',receipt);json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file()});json_new(case/'run-status.json',dict(stage='prepared_native_initial_copy_30_outer',solver_started=False,outer_correctors=30,physical_time_s=0,production_admitted=False));verify_inputs(case)
 return dict(case=str(case),prepared=True,solver_started=False,copied_files=len(rows),outer_correctors=30,ancestor_numeric_miss_preserved=True,prepared_input_manifest_sha256=sha(case/'provenance/prepared-input-sha256.json'),native_initial_copy_receipt_sha256=sha(case/'provenance/native-initial-copy-receipt.json'))

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,required=True);g=ap.add_mutually_exclusive_group(required=True);g.add_argument('--from-case',type=Path);g.add_argument('--verify-inputs',action='store_true');a=ap.parse_args();print(json.dumps(prepare(a.from_case,a.case) if a.from_case else dict(verified_files=len(verify_inputs(a.case))),indent=2))
if __name__=='__main__':main()