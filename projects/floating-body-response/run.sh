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
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|subsetMesh|snappyHexMesh|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
src=/opt/openfoam14/tutorials/incompressibleVoF/floatingObject
mkdir "$case"
cd "$case"
cp -a "$src/0" "$src/constant" "$src/system" .
foamDictionary system/controlDict -entry endTime -set 2
foamDictionary system/controlDict -entry writeFormat -set ascii
foamDictionary system/controlDict -entry runTimeModifiable -set no
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 2 ]]
[[ $(foamDictionary system/controlDict -entry writeInterval -value | xargs) == 0.1 ]]
printf '\n#includeFunc writeCellCentres\n#includeFunc writeCellVolumes\n' >> system/functions
for file in 0/{U,alpha.water,epsilon,k,nut,p_rgh} constant/{dynamicMeshDict,g,momentumTransport,phaseProperties,physicalProperties.air,physicalProperties.water} system/{blockMeshDict,subsetMeshDict,fvSchemes,fvSolution}; do
    cmp "$file" "$src/$file"
done
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
subsetMesh -noFields > subset.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
[[ $(awk '$1=="cells:" {print $2; exit}' check-mesh.log) == 11640 ]]
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
test -s 2/alpha.water
test -s 2/C
test -s 2/Vc
