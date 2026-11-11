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
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|snappyHexMesh|surfaceFeatures)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
src=/opt/openfoam14/tutorials/incompressibleFluid/rotorDisk
mkdir "$case"
cd "$case"
cp -a "$src/0" "$src/constant" "$src/system" .
for file in constant/fvModels constant/momentumTransport constant/physicalProperties constant/geometry/fixed.obj constant/geometry/rotatingZone.obj system/blockMeshDict system/snappyHexMeshDict system/fvSchemes system/fvSolution system/controlDict; do
    cmp "$file" "$src/$file"
done
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 1000 ]]
[[ $(foamDictionary system/controlDict -entry writeInterval -value | xargs) == 50 ]]
[[ $(foamDictionary constant/fvModels -entry disk/type -value | xargs) == rotorDisk ]]
[[ $(foamDictionary constant/fvModels -entry disk/cellZone -value | xargs) == rotatingZone ]]
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
surfaceFeatures > features.log 2>&1
snappyHexMesh > snappy.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
cells=$(awk '$1=="cells:" {print $2; exit}' check-mesh.log)
[[ $cells =~ ^[0-9]+$ ]] && [[ $cells -gt 0 && $cells -le 250000 ]]
test -s constant/polyMesh/cellZones
grep -aq 'rotatingZone' constant/polyMesh/cellZones
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
latest=$(foamListTimes -latestTime | tail -n 1)
[[ $latest =~ ^[1-9][0-9]*$ ]]
test -s "$latest/U"
