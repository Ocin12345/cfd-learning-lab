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
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|snappyHexMesh|potentialFoam|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
src=/opt/openfoam14/tutorials/incompressibleFluid/boxTurb16
mkdir "$case"
cd "$case"
cp -a "$src/0" "$src/constant" "$src/system" .
cp 0/U.orig 0/U
foamDictionary system/controlDict -entry writePrecision -set 12
foamDictionary system/controlDict -entry runTimeModifiable -set no
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 10 ]]
[[ $(foamDictionary system/controlDict -entry deltaT -value | xargs) == 0.025 ]]
[[ $(foamDictionary system/controlDict -entry writeInterval -value | xargs) == 0.25 ]]
cmp 0/U "$src/0/U.orig"
for file in 0/p constant/{fvModels,momentumTransport,physicalProperties} system/{blockMeshDict,functions,fvSchemes,fvSolution}; do
    cmp "$file" "$src/$file"
done
[[ $(grep -Ec '^[[:space:]]*tolerance[[:space:]]+1e-05;' system/fvSolution) == 1 ]]
sed -i -E 's/(^[[:space:]]*tolerance[[:space:]]+)1e-05;/\11e-08;/' system/fvSolution
[[ $(grep -Ec '^[[:space:]]*tolerance[[:space:]]+1e-08;' system/fvSolution) == 1 ]]
[[ $(grep -Ec '^[[:space:]]*tolerance[[:space:]]+1e-06;' system/fvSolution) == 1 ]]
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
[[ $(awk '$1=="cells:" {print $2; exit}' check-mesh.log) == 4096 ]]
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
test -s 10/U
test -s 10/enstrophy
