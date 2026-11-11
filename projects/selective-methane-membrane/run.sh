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
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|snappyHexMesh|createBaffles|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
src=/opt/openfoam14/tutorials/multicomponentFluid/membrane
mkdir "$case"
cd "$case"
cp -a "$src/0" "$src/constant" "$src/system" .
foamDictionary system/controlDict -entry endTime -set 0.5
foamDictionary system/controlDict -entry writeInterval -set 0.05
foamDictionary system/controlDict -entry writePrecision -set 10
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 0.5 ]]
[[ $(foamDictionary system/controlDict -entry writeInterval -value | xargs) == 0.05 ]]
for file in 0/CH4 0/O2 0/N2 0/U 0/T 0/include constant/physicalProperties constant/thermophysicalTransport constant/geometry/boundaries.stl.gz constant/geometry/membrane.stl.gz system/blockMeshDict system/snappyHexMeshDict system/createBafflesDict system/functions system/fvSchemes system/fvSolution; do
    cmp "$file" "$src/$file"
done
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
snappyHexMesh > snappy.log 2>&1
createBaffles > baffles.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
cells=$(awk '$1=="cells:" {print $2; exit}' check-mesh.log)
[[ $cells =~ ^[0-9]+$ ]] && [[ $cells -gt 0 && $cells -le 150000 ]]
grep -aq 'membranePipe' constant/polyMesh/boundary
grep -aq 'membraneSleeve' constant/polyMesh/boundary
test -s constant/polyMesh/cellZones
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
test -s 0.5/CH4
