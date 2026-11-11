#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
case=${1:?Pass an absolute path to a new case directory}
[[ $case == /* ]]
root=$(dirname "$case")
admit() {
    [[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
    [[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
src=/opt/openfoam14/tutorials/fluid/helmholtzResonance
mkdir "$case"
cd "$case"
cp -a "$src/0" "$src/constant" "$src/system" .
cp system/blockMeshDict.resolvedBlocks system/blockMeshDict.caseBlocks
cp system/blockMeshDict.resolvedBoundary system/blockMeshDict.caseBoundary
cmp system/blockMeshDict.resolvedBlocks system/blockMeshDict.caseBlocks
cmp system/blockMeshDict.resolvedBoundary system/blockMeshDict.caseBoundary
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
set -o noclobber
blockMesh > mesh.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
test -s 0.05/p
test -s postProcessing/probes/0/p
