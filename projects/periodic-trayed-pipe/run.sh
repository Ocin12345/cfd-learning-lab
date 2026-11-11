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
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|createZones|createBaffles|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
src=/opt/openfoam14/tutorials/incompressibleVoF/trayedPipe
mkdir "$case"
cd "$case"
cp -a "$src/0" "$src/constant" "$src/system" .
foamDictionary system/controlDict -entry endTime -set 0.25
foamDictionary system/controlDict -entry writeInterval -set 0.025
foamDictionary system/controlDict -entry writeFormat -set ascii
foamDictionary system/controlDict -entry writePrecision -set 10
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 0.25 ]]
[[ $(foamDictionary system/controlDict -entry writeInterval -value | xargs) == 0.025 ]]
[[ $(foamDictionary system/blockMeshDict -entry nBaffles -value | xargs) == 6 ]]
for file in 0/U 0/p_rgh 0/alpha.water constant/phaseProperties constant/g system/blockMeshDict system/createZonesDict system/createBafflesDict system/fvSchemes system/fvSolution; do
    cmp "$file" "$src/$file"
done
[[ $(grep -Ec '^[[:space:]]*MULESCorr[[:space:]]+yes;' system/fvSolution) == 1 ]]
sed -i -E 's/(MULESCorr[[:space:]]+)yes;/\1no;/' system/fvSolution
foamDictionary system/controlDict -entry maxCo -set 0.45
foamDictionary system/controlDict -entry maxAlphaCo -set 0.45
foamDictionary system/fvSolution -entry solvers/p_rghFinal/tolerance -set 1e-11
solver_settings=$(foamDictionary system/fvSolution -entry solvers -value)
grep -Eq 'MULESCorr[[:space:]]+no;' <<< "$solver_settings"
[[ $(foamDictionary system/fvSolution -entry solvers/p_rghFinal/tolerance -value | xargs) == 1e-11 ]]
[[ $(foamDictionary system/fvSolution -entry solvers/p_rghFinal/relTol -value | xargs) == 0 ]]
[[ $(foamDictionary system/controlDict -entry maxCo -value | xargs) == 0.45 ]]
[[ $(foamDictionary system/controlDict -entry maxAlphaCo -value | xargs) == 0.45 ]]
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
createZones > zones.log 2>&1
createBaffles > baffles.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
grep -Eq '^ *cells: +4800 *$' check-mesh.log
test -s constant/polyMesh/boundary
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
test -s 0.25/alpha.water
