#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
root=${1:?Pass an absolute path to an existing workspace directory}
[[ $root == /* && -d $root ]]
parent="$root/radiating-room-v001"
case="$root/radiating-room-v002"
admit() {
    [[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
    [[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|extrudeMesh|snappyHexMesh|potentialFoam|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
(
    cd "$parent"
    sha256sum -c input-sha256.txt
    grep -qx End solver.log
    grep -qx 'Time = 1000s' solver.log
    [[ -s 1000/T && -s 1000/G && -s 1000/phi ]]
)
mkdir "$case"
cd "$case"
cp -a "$parent/"{0,1000,constant,system} .
cp "$parent/"{check-mesh.log,mesh.log,centres.log,volumes.log,pre-mesh-input-sha256.txt,pre-mesh-inputs-verified.log} .
cp "$parent/solver.log" initial-1000-solver.log
set -o noclobber
(
    cd "$parent"
    find 0 1000 constant system -type f -print0 | sort -z | xargs -0 sha256sum
) > inherited-input-sha256.txt
sha256sum -c inherited-input-sha256.txt > inherited-inputs-verified.log
foamDictionary system/controlDict -entry startFrom -set startTime
foamDictionary system/controlDict -entry startTime -set 1000
foamDictionary system/controlDict -entry endTime -set 1250
[[ $(foamDictionary system/controlDict -entry startTime -value | xargs) == 1000 ]]
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 1250 ]]
find 0 1000 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
