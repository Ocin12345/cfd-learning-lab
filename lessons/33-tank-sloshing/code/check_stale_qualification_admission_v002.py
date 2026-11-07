"""Dry-run regression: refuse qualification sources changed during a simulated RAM wait.
Uses copied actual audits and one copied native input. Never launches Fluent.
"""
from pathlib import Path
import ast,datetime,importlib.util,json,os,shutil,sys,traceback,types
R=Path('E:/CFD-Learning-Lab');P=R/'local/sloshing-production-v001';driver=Path(sys.argv[1])
spec=importlib.util.spec_from_file_location('qualification_admission_v001',P/'qualification_admission_v001.py');admission=importlib.util.module_from_spec(spec);spec.loader.exec_module(admission)
fixture=R/'temp/stale-qualification-regression-v001'/driver.stem/datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f')
fixture.mkdir(parents=True,exist_ok=False)
for name in admission.NAMES:shutil.copyfile(admission.QUAL/name,fixture/name)
audit_path=fixture/'openfoam-flat-qualification-v001.json';audit=json.loads(audit_path.read_text())
key=next(k for k in audit['native_source_sha256'] if k.replace('\\','/').endswith('/0/U'))
copy=fixture/'copied-initial-U';shutil.copyfile(admission.source_path(key),copy)
expected=audit['native_source_sha256'].pop(key);audit['native_source_sha256'][str(copy)]=expected;audit_path.write_text(json.dumps(audit))
admission.QUAL=fixture;sys.modules['qualification_admission_v001']=admission
admission.verify_qualifications() # Real archive verification passes before the simulated wait.
waiting={'released':False};wing=types.ModuleType('wing_fluent_v001');wing.available_memory=lambda:5*2**30 if waiting['released'] else 0;sys.modules['wing_fluent_v001']=wing
launched=[]
class InterceptedLaunch(Exception):pass
def launch(**kwargs):
    launched.append(True);raise InterceptedLaunch('Intercepted before real Fluent launch')
def wait(seconds):
    waiting['released']=True
    copy.write_bytes(copy.read_bytes()+b'\n// test mutation during memory wait\n')
def write_json(path,payload):path.write_text(json.dumps(payload,default=str))
function=next(n for n in ast.parse(driver.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='main')
namespace=dict(__file__=str(driver),os=os,shutil=shutil,time=types.SimpleNamespace(sleep=wait),pf=types.SimpleNamespace(launch_fluent=launch),traceback=traceback,write_json=write_json,make_mesh=lambda out:out/'NOT_LAUNCHED_TEST_MESH')
exec(compile(ast.Module(body=[function],type_ignores=[]),str(driver),'exec'),namespace)
blocked=False
try:namespace['main'](fixture/'dry-run-output')
except AssertionError as exc:
    assert 'Qualification source changed:' in str(exc),str(exc);blocked=True
except InterceptedLaunch:pass
assert blocked and not launched,'Fluent launch attempted before stale qualifications were refused'
print(json.dumps(dict(production_refused_stale_sources=True,fluent_launch_calls=len(launched),test_fixture=str(fixture),original_qualification_files_changed=False)))
