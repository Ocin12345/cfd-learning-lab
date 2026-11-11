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
tank=/opt/openfoam14/tutorials/fluid/decompressionTank/decompressionTank
gas=/opt/openfoam14/tutorials/fluid/helmholtzResonance
mkdir "$case"
cd "$case"
sha256sum "$tank/0/p" "$tank/0/T" "$tank/0/U" "$tank/system/blockMeshDict" \
    "$tank/system/fvSchemes" "$tank/system/fvSolution" "$tank/system/controlDict" \
    "$tank/constant/momentumTransport" "$gas/constant/physicalProperties" > tutorial-sha256.txt
cp -a "$tank/0" "$tank/constant" "$tank/system" .
cp "$gas/constant/physicalProperties" constant/physicalProperties
cmp "$gas/constant/physicalProperties" constant/physicalProperties
foamDictionary 0/p -entry internalField -set 'uniform 200000'
foamDictionary 0/p -entry boundaryField/nozzle/value -set 'uniform 100000'
foamDictionary system/controlDict -entry endTime -set 0.0015
foamDictionary system/controlDict -entry writeInterval -set 200
foamDictionary system/controlDict -entry writePrecision -set 10
[[ $(foamDictionary 0/p -entry internalField -value | xargs) == 'uniform 200000' ]]
[[ $(foamDictionary 0/p -entry boundaryField/nozzle/value -value | xargs) == 'uniform 100000' ]]
[[ $(foamDictionary 0/p -entry boundaryField/nozzle/type -value | xargs) == fixedValue ]]
[[ $(foamDictionary constant/physicalProperties -entry thermoType/equationOfState -value | xargs) == perfectGas ]]
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 0.0015 ]]
[[ $(foamDictionary system/controlDict -entry writeInterval -value | xargs) == 200 ]]
[[ $(foamDictionary system/controlDict -entry deltaT -value | xargs) == 5e-07 ]]
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
sha256sum -c tutorial-sha256.txt > tutorial-inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
test -s 0.0015/p
