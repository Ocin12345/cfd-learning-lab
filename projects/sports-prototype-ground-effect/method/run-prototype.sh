#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
case=${1:?Pass absolute new case directory}
angle=${2:?Pass diffuser angle in degrees}
[[ $case == /* ]]
root=$(dirname "$case")
src=/opt/openfoam14/tutorials/incompressibleFluid/motorBikeSteady
generator=$(dirname "$(realpath "$0")")/prepare-prototype.py
admit() {
    [[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
    [[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|snappyHexMesh|surfaceCheck|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "Conflicting solver or mesh process" >&2
        exit 3
    fi
}
admit
python3 "$generator" "$case" "$angle"
cd "$case"
cp "$src/system/"{fvSchemes,fvSolution} system/
set -o noclobber
sha256sum "$generator" > generator-sha256.txt
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-sha256.txt
for part in chassis frontWheels rearWheels; do
    surfaceCheck "constant/geometry/$part.stl" > "surface-$part.log" 2>&1
    grep -Fq 'Surface is closed.' "surface-$part.log"
done
blockMesh > block-mesh.log 2>&1
snappyHexMesh -overwrite > snappy.log 2>&1
grep -Eq '^[[:space:]]*End[[:space:]]*$' snappy.log
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
cells=$(awk '$1=="cells:" {print $2; exit}' check-mesh.log)
[[ $cells -ge 50000 && $cells -le 240000 ]]
for patch in inlet outlet road farfield chassis frontWheels rearWheels; do
    grep -Eq "^[[:space:]]+$patch[[:space:]]*$" constant/polyMesh/boundary
done
sha256sum -c pre-mesh-sha256.txt > original-inputs-verified.log
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
sha256sum -c generator-sha256.txt >> inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -Eq '^[[:space:]]*End[[:space:]]*$' solver.log
test -s 600/U
test -s 600/phi
