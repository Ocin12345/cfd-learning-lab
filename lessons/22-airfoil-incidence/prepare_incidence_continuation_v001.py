from pathlib import Path
import hashlib,json,shutil
root=Path("E:/CFD-Learning-Lab/native/runs/2026-11-06-airfoil-incidence-paired-v001")
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for angle in (0,5):
    source=root/f"openfoam-alpha{angle}"
    copy=root/f"openfoam-alpha{angle}-v002"
    assert not copy.exists()
    assert json.loads((source/"run-status.json").read_text())["stage"]=="finished_execution"
    for rank in ("processor0","processor1"):
        for field in ("U","p","phi"):
            assert (source/rank/"3000"/field).is_file()
    watched={p.relative_to(source).as_posix():h(p) for p in source.rglob("*") if p.is_file()}
    shutil.copytree(source,copy)
    assert all(h(copy/n)==d and h(source/n)==d for n,d in watched.items())
    control=copy/"system/controlDict"
    s=control.read_text()
    changes=(("startFrom startTime;","startFrom latestTime;"),
             ("endTime 3000;","endTime 4800;"),("writeInterval 3000;","writeInterval 4800;"))
    for a,b in changes:
        assert s.count(a)==1
        s=s.replace(a,b)
    control.write_text(s,encoding="ascii",newline="\n")
    (copy/"continuation-inputs.json").write_text(json.dumps({"angle_degrees":angle,
        "source_start_iteration":3000,"end_iteration":4800,"no_reinitialization":True,
        "source_manifest":watched,"revised_control_sha256":h(control)},indent=2)+"\n")
    print(angle,"degrees; exact native copies:",len(watched))
