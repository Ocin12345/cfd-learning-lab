#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
case=${1:?Pass an absolute path to a new case directory}
[[ $case == /* ]]
root=$(dirname "$case")
projector=$(dirname "$(realpath "$0")")/project-sphere-points.py

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
geometry
{
    sphereSurface
    {
        type sphere;
        centre (0 0 0);
        radius 1;
    }
}
vertices
(
    (-0.5773502691896258 -0.5773502691896258 -0.5773502691896258)
    (0.5773502691896258 -0.5773502691896258 -0.5773502691896258)
    (0.5773502691896258 0.5773502691896258 -0.5773502691896258)
    (-0.5773502691896258 0.5773502691896258 -0.5773502691896258)
    (-0.5773502691896258 -0.5773502691896258 0.5773502691896258)
    (0.5773502691896258 -0.5773502691896258 0.5773502691896258)
    (0.5773502691896258 0.5773502691896258 0.5773502691896258)
    (-0.5773502691896258 0.5773502691896258 0.5773502691896258)
    (-18 -18 -18)
    (18 -18 -18)
    (18 18 -18)
    (-18 18 -18)
    (-18 -18 18)
    (18 -18 18)
    (18 18 18)
    (-18 18 18)
);
blocks
(
    hex (0 4 7 3 8 12 15 11) (16 16 24) simpleGrading (1 1 80)
    hex (1 2 6 5 9 10 14 13) (16 16 24) simpleGrading (1 1 80)
    hex (0 1 5 4 8 9 13 12) (16 16 24) simpleGrading (1 1 80)
    hex (3 7 6 2 11 15 14 10) (16 16 24) simpleGrading (1 1 80)
    hex (0 3 2 1 8 11 10 9) (16 16 24) simpleGrading (1 1 80)
    hex (4 5 6 7 12 13 14 15) (16 16 24) simpleGrading (1 1 80)
);
edges
(
    project 0 4 (sphereSurface)
    project 4 7 (sphereSurface)
    project 3 7 (sphereSurface)
    project 0 3 (sphereSurface)
    project 1 2 (sphereSurface)
    project 2 6 (sphereSurface)
    project 5 6 (sphereSurface)
    project 1 5 (sphereSurface)
    project 0 1 (sphereSurface)
    project 4 5 (sphereSurface)
    project 6 7 (sphereSurface)
    project 2 3 (sphereSurface)
);
faces
(
    project (0 4 7 3) sphereSurface
    project (1 2 6 5) sphereSurface
    project (0 1 5 4) sphereSurface
    project (3 7 6 2) sphereSurface
    project (0 3 2 1) sphereSurface
    project (4 5 6 7) sphereSurface
);
boundary
(
    inlet { type patch; faces ((8 12 15 11)); }
    outlet { type patch; faces ((9 10 14 13)); }
    farfield { type patch; faces ((8 9 13 12) (11 15 14 10) (8 11 10 9) (12 13 14 15)); }
    sphere { type wall; faces ((3 7 4 0) (5 6 2 1) (4 5 1 0) (2 6 7 3) (1 2 3 0) (7 6 5 4)); }
);
EOF
cat > constant/physicalProperties <<'EOF'
FoamFile { format ascii; class dictionary; object physicalProperties; }
viscosityModel constant;
nu 1;
EOF
cat > 0/U <<'EOF'
FoamFile { format ascii; class volVectorField; location "0"; object U; }
dimensions [velocity];
internalField uniform (0.05 0 0);
boundaryField
{
    inlet { type fixedValue; value uniform (0.05 0 0); }
    outlet { type zeroGradient; }
    farfield { type fixedValue; value uniform (0.05 0 0); }
    sphere { type noSlip; }
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
    farfield { type zeroGradient; }
    sphere { type zeroGradient; }
}
EOF
cat > system/functions <<'EOF'
FoamFile { format ascii; class dictionary; object functions; }
sphereForce
{
    type forces;
    libs ("libforces.so");
    patches (sphere);
    rho rhoInf;
    rhoInf 1;
    CofR (0 0 0);
    pitchAxis (0 1 0);
    writeControl timeStep;
    writeInterval 1;
}
EOF
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > pre-mesh-input-sha256.txt
sha256sum -c pre-mesh-input-sha256.txt > pre-mesh-inputs-verified.log
set -o noclobber
blockMesh > mesh.log 2>&1
python3 "$projector" . --write > projection.log
checkMesh > check-mesh.log 2>&1
grep -Fq 'Mesh OK.' check-mesh.log
[[ $(awk '$1=="cells:" {print $2; exit}' check-mesh.log) == 36864 ]]
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
sha256sum preprojection-points sphere-projection.json projection.log > projection-sha256.txt
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
