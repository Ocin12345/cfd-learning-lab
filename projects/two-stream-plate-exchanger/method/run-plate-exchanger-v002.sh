#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
case=${1:?Provide new absolute case directory}
orientation=${2:?counter or parallel}
level=${3:?base or fine}
[[ $case == /* ]]
root=$(dirname "$case")
generator=$(dirname "$(realpath "$0")")/prepare-plate-exchanger.py
admit() {
    [[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
    [[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|splitMeshRegions|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "Conflicting native process" >&2
        exit 3
    fi
}
admit
python3 "$generator" "$case" "$orientation" "$level"
cd "$case"
set -o noclobber
sha256sum "$generator" > generator-sha256.txt
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-sha256.txt
blockMesh > mesh.log 2>&1
splitMeshRegions -cellZones all -noFields > split.log 2>&1
sha256sum -c pre-mesh-sha256.txt > original-inputs-verified.log
for region in cold plate hot; do
    checkMesh -region "$region" > "check-$region.log" 2>&1
    grep -Fq 'Mesh OK.' "check-$region.log"
    count=$(awk '$1=="cells:" {print $2; exit}' "check-$region.log")
    expected=1920
    if [[ $level == fine ]]; then expected=4320; fi
    if [[ $region == plate ]]; then expected=$((expected/4)); fi
    [[ $count == "$expected" ]]
    foamPostProcess -region "$region" -func writeCellCentres -time 0 > "centres-$region.log" 2>&1
    foamPostProcess -region "$region" -func writeCellVolumes -time 0 > "volumes-$region.log" 2>&1
    test -s "0/$region/C"
    test -s "0/$region/Vc"
done
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
sha256sum -c generator-sha256.txt >> inputs-verified.log
admit
foamMultiRun > solver.log 2>&1
grep -Eq '^[[:space:]]*End[[:space:]]*$' solver.log
test -s 900/hot/T
test -s 900/cold/phi
test -s 900/plate/T
