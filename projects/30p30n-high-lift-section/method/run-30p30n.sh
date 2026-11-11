#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
case=${1:?Provide a new absolute case directory}
source_grid=${2:?Provide the pinned compressed APC-IV SU2 grid}
[[ $case == /* && $source_grid == /* ]]
root=$(dirname "$case")
generator=$(dirname "$(realpath "$0")")/prepare-30p30n.py
admit() {
    [[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
    [[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|snappyHexMesh|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "Conflicting native process" >&2
        exit 3
    fi
}
admit
python3 "$generator" "$case" "$source_grid"
cd "$case"
cp /opt/openfoam14/tutorials/incompressibleFluid/airFoil2D/system/fvSchemes system/
set -o noclobber
sha256sum "$generator" "$source_grid" > source-input-sha256.txt
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
[[ $(awk '$1=="cells:" {print $2; exit}' check-mesh.log) == 63957 ]]
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
sha256sum -c source-input-sha256.txt >> inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -Eq '^[[:space:]]*End[[:space:]]*$' solver.log
test -s 1600/U
test -s 1600/phi
