#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
phase=${1:?Specify local, steady, or blended}
parent=${2:?Specify original 600-iteration case}
source_case=${3:?Specify checkpoint case}
case=${4:?Specify fresh absolute case path}
audit=${5:?Specify original audit JSON}
[[ $case == /* && $case != "$parent" && $case != "$source_case" ]]
root=$(dirname "$case")
script_dir=$(dirname "$(realpath "$0")")
admit() {
    [[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
    [[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|snappyHexMesh|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo 'Conflicting native process' >&2
        exit 3
    fi
}
admit
python3 "$script_dir/prepare-flap-restart-v002.py" "$phase" "$parent" "$source_case" "$case" "$audit"
cd "$case"
set -o noclobber
sha256sum -c input-sha256.txt > inputs-verified.log
checkMesh > mesh-readback.log 2>&1
grep -Fq 'Mesh OK.' mesh-readback.log
sha256sum -c input-sha256.txt >> inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -Eq '^[[:space:]]*End[[:space:]]*$' solver.log
end=1800
if [[ $phase == local ]]; then
    end=1500
    grep -Fq 'Using LTS' solver.log
elif [[ $phase == blended ]]; then
    end=1200
fi
test -s "$end/U"
test -s "$end/phi"
