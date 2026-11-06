"""Read-only balanced-dictionary validation continuation on already completed annular native preflight; never invokes native commands."""
from pathlib import Path
import argparse,json,re,types
import prepare_openfoam_annular_couette_v002 as prep
PREPARATION_SHA='9e07e2982beca3a34242478bb23bf98e7d8ccb14b00836e69831ffc713bdfe34'
def balanced_blocks(path):
 text=prep.clean(Path(path).read_text());i=text.index('{',text.index('boundaryField'))+1;result={}
 while i<len(text):
  m=re.match(r'\s*("[^"]+"|\w+)\s*\{',text[i:])
  if not m:break
  name=m[1].strip('"');a=i+m.end();j=a;depth=1
  while depth:
   prep.require(j<len(text),'Unclosed native patch dictionary');depth+=(text[j]=='{')-(text[j]=='}');j+=1
  body=text[a:j-1];flat=[];nested=0
  for c in body:
   if c=='{':nested+=1
   elif c=='}':nested-=1
   elif not nested:flat.append(c)
  result[name]=''.join(flat);i=j
 return result

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();case=prep.checked_case(a.run);out=a.out.resolve();prep.require(out.drive.upper()=='E:' and not out.exists() and case not in out.parents,'Fresh external E output');prep.require(prep.sha(Path(prep.__file__))==PREPARATION_SHA,'Original frozen preparer changed')
 for n in ('native-mapping.json','preflight-validation.json','preflight-readonly-recovery.json'):prep.require(not (case/'provenance'/n).exists(),'Preserve existing derived preflight output '+n)
 before={p.relative_to(case).as_posix():prep.sha(p) for p in case.rglob('*') if p.is_file()};prep.require(not (case/'solver.log').exists(),'NoCFDmusthavebeenadvanced')
 # Isolate only the field-parser binding. All geometry/physics/targets and original validation bytecode are unchanged.
 scope=dict(prep.validate_native.__globals__);scope['blocks']=balanced_blocks;validate=types.FunctionType(prep.validate_native.__code__,scope,'validate_native_balanced',prep.validate_native.__defaults__,prep.validate_native.__closure__)
 result=validate(case)
 prep.require(all(prep.sha(case/n)==h for n,h in before.items()),'Original native/source/log/status input changed during read-only continuation')
 recovery=dict(readonly_validation_passed=True,native_preflight_commands_completed=True,native_preflight_initial_validator_failed='Flat patch regex omitted inner with nested omega Function1 dictionary',original_preparation_source_sha256=PREPARATION_SHA,readonly_validator_source_sha256=prep.sha(Path(__file__)),preflight_receipt_sha256=prep.sha(case/'provenance/preflight-validation.json'),native_mapping_sha256=prep.sha(case/'provenance/native-mapping.json'),all_original_input_files_unchanged=True,original_input_file_count=len(before),raw_failed_run_status_preserved=json.loads((case/'run-status.json').read_text()),solver_started=False,CFD_replayed=False,native_commands_replayed=False)
 prep.json_new(case/'provenance/preflight-readonly-recovery.json',recovery);prep.copied(Path(__file__),case/'provenance'/Path(__file__).name);out.mkdir(parents=True);prep.json_new(out/'recovery.json',recovery);prep.json_new(out/'preflight-validation.json',result);prep.json_new(out/'source-manifest.json',dict(sources=[dict(path=str(case/n),sha256_before=h,sha256_after=prep.sha(case/n)) for n,h in before.items()],inputs_unchanged_during_analysis=True,solver_started_or_modified=False));print(json.dumps(dict(result=result,recovery=recovery),indent=2))
if __name__=='__main__':main()
