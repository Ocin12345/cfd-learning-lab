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
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|foamPostProcess)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}

admit
mkdir "$case"
cd "$case"
mkdir 0 constant system
src=/opt/openfoam14/tutorials/incompressibleFluid/venturiTube
cp "$src/system/"{fvSchemes,fvSolution} system/
cp "$src/constant/momentumTransport" constant/
grep -Eq '^simulationType[[:space:]]+laminar;' constant/momentumTransport
grep -Fq 'default         steadyState;' system/fvSchemes
cat > system/controlDict <<'EOF'
FoamFile { format ascii; class dictionary; object controlDict; }
solver incompressibleFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 1000;
deltaT 1;
writeControl timeStep;
writeInterval 25;
purgeWrite 0;
writeFormat ascii;
writePrecision 12;
writeCompression off;
runTimeModifiable no;
EOF
cat > system/blockMeshDict <<'EOF'
FoamFile { format ascii; class dictionary; object blockMeshDict; }
vertices
(
    (-0.25 0 -0.05) (0 0 -0.05) (2 0 -0.05)
    (-0.25 0.6 -0.05) (0 0.6 -0.05) (2 0.6 -0.05)
    (-0.25 0 0.05) (0 0 0.05) (2 0 0.05)
    (-0.25 0.6 0.05) (0 0.6 0.05) (2 0.6 0.05)
);
blocks
(
    hex (0 1 4 3 6 7 10 9) (20 96 1) simpleGrading (1 15 1)
    hex (1 2 5 4 7 8 11 10) (160 96 1) simpleGrading (1 15 1)
);
boundary
(
    inlet { type patch; faces ((0 6 9 3)); }
    outlet { type patch; faces ((2 5 11 8)); }
    farfield { type patch; faces ((3 9 10 4) (4 10 11 5)); }
    upstream { type patch; faces ((0 1 7 6)); }
    plate { type wall; faces ((1 2 8 7)); }
    frontAndBack { type empty; faces ((0 3 4 1) (1 4 5 2) (6 7 10 9) (7 8 11 10)); }
);
EOF
cat > constant/physicalProperties <<'EOF'
FoamFile { format ascii; class dictionary; object physicalProperties; }
viscosityModel constant;
nu 0.00025;
EOF
cat > 0/U <<'EOF'
FoamFile { format ascii; class volVectorField; location "0"; object U; }
dimensions [velocity];
internalField uniform (1 0 0);
boundaryField
{
    inlet { type fixedValue; value uniform (1 0 0); }
    outlet { type zeroGradient; }
    upstream { type slip; }
    farfield { type slip; }
    plate { type noSlip; }
    frontAndBack { type empty; }
}
EOF
cat > 0/p <<'EOF'
FoamFile { format ascii; class volScalarField; location "0"; object p; }
dimensions [kinematicPressure];
internalField uniform 0;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    upstream { type zeroGradient; }
    farfield { type zeroGradient; }
    plate { type zeroGradient; }
    frontAndBack { type empty; }
}
EOF
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
[[ $(awk '$1=="cells:" {print $2; exit}' check-mesh.log) == 17280 ]]
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
