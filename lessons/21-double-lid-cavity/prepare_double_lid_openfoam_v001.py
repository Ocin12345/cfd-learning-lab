from pathlib import Path
import hashlib,json,shutil

lab=Path("E:/CFD-Learning-Lab")
source=lab/"native/runs/2026-11-01-cavity-v001/openfoam"
case=lab/"native/runs/2026-11-06-double-lid-cavity-v001/openfoam"
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert not case.exists()
pre=json.loads((source/"provenance/preflight-validation.json").read_text(encoding="utf-8-sig"))
assert pre["passed"] and pre["mesh_location"]=="0/polyMesh"
relative=["0/U","0/p","system/controlDict","system/fvSchemes","system/fvSolution",
          "system/decomposeParDict","constant/physicalProperties","constant/momentumTransport"]
relative+=["0/polyMesh/"+p.name for p in (source/"0/polyMesh").iterdir() if p.is_file()]
before={name:h(source/name) for name in relative}
for name,digest in pre["native_mesh_sha256"].items():
    assert before["0/polyMesh/"+name]==digest
for name in relative:
    target=case/name
    target.parent.mkdir(exist_ok=True,parents=True)
    shutil.copyfile(source/name,target)
    assert h(target)==before[name]
for name,old,new in [
  ("0/U","bottom { type noSlip; }","bottom { type fixedValue; value uniform (1 0 0); }"),
  ("system/controlDict","writeInterval 100;","writeInterval 200;")]:
    p=case/name
    text=p.read_text(encoding="ascii")
    assert text.count(old)==1
    p.write_text(text.replace(old,new),encoding="ascii",newline="\n")
assert before=={name:h(source/name) for name in relative}
(case/"input-manifest.json").write_text(json.dumps(
    {"source_unchanged":True,"source_hashes":before,"prepared_hashes":{name:h(case/name) for name in relative}},indent=2)+"\n",encoding="utf-8")
print("Prepared fresh OpenFOAM case; copied and matched",len(relative),"files.")
