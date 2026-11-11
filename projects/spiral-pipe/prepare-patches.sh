#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
root=${1:?Pass an absolute path to an existing workspace directory}
[[ $root == /* && -d $root ]]
case="$root/spiral-pipe-mesh-v002"
parent="$root/spiral-pipe-mesh-v001"
[[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
[[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|extrudeMesh|snappyHexMesh|potentialFoam|foamPostProcess|foamToVTK)$/ {busy=1} END {exit !busy}'; then
    echo "A conflicting solver or mesh process is active" >&2
    exit 3
fi
(
    cd "$parent"
    sha256sum -c mesh-sha256.txt
    grep -Fq 'Mesh OK.' check-mesh.log
)
mkdir "$case"
cd "$case"
mkdir -p constant/geometry system 0
cp "$parent/constant/geometry/curve.obj.gz" constant/geometry/
cp "$parent/system/"{controlDict,blockMeshDict,extrudeMeshDict} system/
python3 - <<'PY'
from pathlib import Path
p = Path("system/blockMeshDict")
text = p.read_text()
old = "defaultPatch\n{\n    type patch;\n}"
new = "defaultPatch\n{\n    name walls;\n    type wall;\n}"
assert text.count(old) == 1
text = text.replace(old, new)
old = "boundary\n(\n    outlet"
new = """boundary
(
    inlet
    {
        type patch;
        faces
        (
            (0 1 5 4)
            (6 2 0 4)
            (3 7 5 1)
            (6 7 3 2)
            (2 3 1 0)
        );
    }
    outlet"""
assert text.count(old) == 1
p.write_text(text.replace(old, new))
PY
cmp "$parent/constant/geometry/curve.obj.gz" constant/geometry/curve.obj.gz
cmp "$parent/system/extrudeMeshDict" system/extrudeMeshDict
find constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
extrudeMesh > extrude.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
[[ $(awk '$1=="cells:" {print $2; exit}' check-mesh.log) == 112448 ]]
cmp "$parent/constant/polyMesh/points" constant/polyMesh/points
foamToVTK -ascii -noInternal -noPointValues -noFunctionObjects -excludePatches '(walls)' -time 0 > patch-export.log 2>&1
find constant system VTK -type f -print0 | sort -z | xargs -0 sha256sum > mesh-sha256.txt
sha256sum -c mesh-sha256.txt > mesh-verified.log
