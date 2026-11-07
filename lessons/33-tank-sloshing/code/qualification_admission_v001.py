"""Read-only qualification admission, shared by Windows and WSL production launchers."""
import hashlib,json,os
from pathlib import Path
ROOT=Path('E:/CFD-Learning-Lab' if os.name=='nt' else '/mnt/e/CFD-Learning-Lab')
QUAL=ROOT/'local/sloshing-tank-v001'
NAMES=('fluent-flat-qualification-v005.json','openfoam-flat-qualification-v001.json','fluent-tilted-qualification-v002.json','openfoam-tilted-qualification-v001.json')
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
def source_path(name):
    normalized=name.replace('\\','/')
    prefix='E:/CFD-Learning-Lab/'
    assert normalized.startswith(prefix), 'Unexpected qualification source root'
    path=(ROOT/normalized[len(prefix):]).resolve()
    assert path.is_relative_to(ROOT.resolve())
    return path
def verify_qualifications():
    paths=[QUAL/name for name in NAMES]
    missing=[p.name for p in paths if not p.is_file()]
    if missing:raise RuntimeError('Production denied; missing qualification audits: '+', '.join(missing))
    verified={}
    for p in paths:
        audit=json.loads(p.read_text())
        assert audit['qualification_pass'] is True and audit['native_field_checks_pass'] is True,p
        assert audit['native_source_sha256'],p
        for name,expected in audit['native_source_sha256'].items():
            assert sha(source_path(name))==expected, 'Qualification source changed: '+name
        verified[p.name]=sha(p)
    return verified
if __name__=='__main__':print(json.dumps(dict(qualification_audit_sha256=verify_qualifications(),production_admission_pass=True)))
