#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
root=${1:?Pass an absolute path to an existing workspace directory}
[[ $root == /* && -d $root ]]
case="$root/spiral-pipe-mesh-v001"
[[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
[[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|extrudeMesh|snappyHexMesh|potentialFoam|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
    echo "A conflicting solver or mesh process is active" >&2
    exit 3
fi
src=/opt/openfoam14/tutorials/mesh/spiralPipe
mkdir "$case"
cd "$case"
cp -a "$src/constant" "$src/system" .
cmp constant/geometry/curve.obj.gz "$src/constant/geometry/curve.obj.gz"
cmp system/blockMeshDict "$src/system/blockMeshDict"
cmp system/extrudeMeshDict "$src/system/extrudeMeshDict"
[[ $(foamDictionary system/extrudeMeshDict -entry extrudeModel/nLayers -value | xargs) == 500 ]]
foamDictionary system/extrudeMeshDict -entry extrudeModel/nLayers -set 250
[[ $(foamDictionary system/extrudeMeshDict -entry extrudeModel/nLayers -value | xargs) == 250 ]]
find constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
extrudeMesh > extrude.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
cells=$(awk '$1=="cells:" {print $2; exit}' check-mesh.log)
[[ $cells =~ ^[0-9]+$ ]] && [[ $cells -gt 100000 && $cells -le 120000 ]]
find constant system -type f -print0 | sort -z | xargs -0 sha256sum > mesh-sha256.txt
sha256sum -c mesh-sha256.txt > mesh-verified.log
