#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
root=${1:?Pass an absolute path to an existing workspace directory}
[[ $root == /* && -d $root ]]
case="$root/radiating-room-v001"
admit() {
    [[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
    [[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|extrudeMesh|snappyHexMesh|potentialFoam|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
src=/opt/openfoam14/tutorials/fluid/hotRadiationRoom
mkdir "$case"
cd "$case"
cp -a "$src/0" "$src/constant" "$src/system" .
foamDictionary system/controlDict -entry writeInterval -set 25
foamDictionary system/controlDict -entry writePrecision -set 12
foamDictionary system/controlDict -entry runTimeModifiable -set no
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 1000 ]]
for file in 0/{G,T,U,alphat,epsilon,k,nut,p,p_rgh} constant/{fvModels,g,momentumTransport,pRef,physicalProperties,radiationProperties} system/{blockMeshDict,fvSchemes,fvSolution}; do
    cmp "$file" "$src/$file"
done
python3 - <<'PY'
from pathlib import Path
import re
p=Path("system/blockMeshDict")
counts={5:3,10:6,15:8,40:20,80:40}
text=p.read_text()
def reduce(match):
    value=tuple(int(i) for i in match[2].split())
    assert len(value)==3
    return match[1]+"("+" ".join(str(counts[i]) for i in value)+")"
new,n=re.subn(r"(^\s*hex\s+\([0-9\s]+\)\s+)\(([0-9\s]+)\)",reduce,text,flags=re.M)
assert n==17
p.write_text(new)
PY
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
[[ $(awk '$1=="cells:" {print $2; exit}' check-mesh.log) == 15523 ]]
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
