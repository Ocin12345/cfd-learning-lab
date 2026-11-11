#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
case=${1:?Provide a new absolute case directory}
[[ $case == /* ]]
root=$(dirname "$case")
projector=$(dirname "$(realpath "$0")")/project-cone-wedge-v002.py
src=/opt/openfoam14/tutorials/shockFluid/biconic25-55Run35
admit() {
    [[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
    [[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|datToFoam|collapseEdges|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "Conflicting native process" >&2
        exit 3
    fi
}
admit
mkdir "$case"
cd "$case"
cp -a "$src/0" "$src/constant" "$src/system" .
diff -qr 0 "$src/0"
diff -qr constant "$src/constant"
diff -qr system "$src/system"
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > tutorial-input-sha256.txt
foamDictionary system/controlDict -entry writeFormat -set ascii
foamDictionary system/controlDict -entry runTimeModifiable -set false
[[ $(foamDictionary system/controlDict -entry endTime -value | xargs) == 0.001 ]]
[[ $(foamDictionary system/controlDict -entry deltaT -value | xargs) == 1e-07 ]]
set -o noclobber
blockMesh > block-mesh.log 2>&1
datToFoam constant/grid256.dat > convert.log 2>&1
points=constant/polyMesh/points
sed -n '1,17{s/binary/ascii/;p;}' "$points" > block-points-header.txt
cat block-points-header.txt constant/points.tmp > native-grid-points.txt
cp native-grid-points.txt "$points"
collapseEdges > collapse.log 2>&1
sed -i '/wedge/,/nFaces/{s/patch/wedge/}' constant/polyMesh/boundary
python3 "$projector" > wedge-projection.json
sha256sum "$projector" preprojection-wedge-points wedge-projection.json > wedge-projection-sha256.txt
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
[[ $(awk '$1=="cells:" {print $2; exit}' check-mesh.log) == 32385 ]]
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
sha256sum -c wedge-projection-sha256.txt >> inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -Eq '^[[:space:]]*End[[:space:]]*$' solver.log
test -s 0.001/U
test -s 0.001/p
test -s 0.001/T
