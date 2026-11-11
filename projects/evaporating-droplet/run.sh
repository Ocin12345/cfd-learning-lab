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
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|createZones|transformPoints|createNonConfo|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
src=/opt/openfoam14/tutorials/multicomponentFluid/parcelInBox
mkdir "$case"
cd "$case"
cp -a "$src/0" "$src/constant" "$src/system" .
foamDictionary system/controlDict -entry endTime -set 0.1
foamDictionary system/controlDict -entry writeInterval -set 0.001
foamDictionary system/controlDict -entry writePrecision -set 12
foamDictionary system/controlDict -entry runTimeModifiable -set no
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 0.1 ]]
[[ $(foamDictionary system/controlDict -entry writeInterval -value | xargs) == 0.001 ]]
for file in 0/{G,H2O,T,U,air,p} constant/{cloudPositions,cloudProperties,fvModels,g,momentumTransport,physicalProperties,speciesThermo} system/{blockMeshDict,fvSchemes,fvSolution}; do
    cmp "$file" "$src/$file"
done
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
[[ $(awk '$1=="cells:" {print $2; exit}' check-mesh.log) == 125 ]]
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
test -s 0.1/T
