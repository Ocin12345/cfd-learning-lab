from pathlib import Path
import argparse,json,re,shutil
import analyse_openfoam_airfoil_v002 as av

def rewrite_control(raw,prior,end):
 if not isinstance(end,int) or end<=prior or av.scalar(av.dictionary_text(raw) if hasattr(av,'dictionary_text') else raw,'endTime')!=prior:raise ValueError('Wrong prior/finite continuation limit')
 revised=raw
 for key,value in [('startFrom','latestTime'),('startTime',str(prior)),('endTime',str(end)),('writeInterval',str(end))]:
  revised,count=re.subn(r'(?<![\w])'+key+r'\s+[^;{}]+;',key+' '+value+';',revised,count=1)
  if count!=1:raise ValueError('Missing root '+key)
 if av.control_invariants(raw)!=av.control_invariants(revised):raise ValueError('Changed physical/force controls')
 return revised

def validate_miss(a,source,prior):
 if prior!=24000 or a.get('iteration')!=prior or a.get('source_case')!=str(source) or a.get('numerical_checks_passed') is not False or a.get('failed_checks')!=['openfoam.residuals']:raise ValueError('Only the verified 24000 residual-only miss is admitted')

def prepare(source,target,audit,endpoint):
 prior=24000
 root=Path('E:/CFD-Learning-Lab/native/runs/showcase-airfoil-trailing-edge-v001').resolve()
 assert source.resolve()==root/'openfoam-level2-v002' and target.resolve()==root/'openfoam-level2-v003' and endpoint==32000
 assert not target.exists() and target.drive.upper()=='E:'
 a=json.loads(audit.read_text());validate_miss(a,source,prior)
 assert json.loads((source/'run-status.json').read_text())['stage']=='finished_execution'
 for p,h in a['source_sha256'].items():assert av.sha256(Path(p))==h,p
 for part in ('processor0','processor1'):
  for f in ('U','p','phi'):assert (source/part/str(prior)/f).is_file()
 hashes={p.relative_to(source).as_posix():av.sha256(p) for p in source.rglob('*') if p.is_file()}
 revised=rewrite_control((source/'system/controlDict').read_text(),prior,endpoint)
 shutil.copytree(source,target)
 assert all(av.sha256(source/n)==h and av.sha256(target/n)==h for n,h in hashes.items())
 shutil.copy2(target/'run-status.json',target/'parent-v002-run-status.json')
 (target/'system/controlDict').write_text(revised,newline='\n')
 record=dict(source_case=str(source),source_start_iteration=prior,end_iteration=endpoint,no_reinitialization=True,source_manifest=hashes,parent_audit=str(audit),parent_audit_sha256=av.sha256(audit),revised_control_sha256=av.sha256(target/'system/controlDict'),only_top_level_execution_controls_changed=True)
 (target/'continuation-inputs-v003.json').write_text(json.dumps(record,indent=2)+'\n')
 (target/'run-status.json').write_text(json.dumps(dict(stage='prepared_continuation_v003',solver_started=False))+'\n')
 print(json.dumps(dict(copied_files=len(hashes),source=str(source),target=str(target),start=prior,end=endpoint,no_reinitialization=True)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--target',type=Path,required=True);p.add_argument('--audit',type=Path,required=True);p.add_argument('--end',type=int,required=True);a=p.parse_args();prepare(a.source,a.target,a.audit,a.end)
