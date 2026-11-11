#!/usr/bin/env bash
set -eo pipefail
source /opt/openfoam14/etc/bashrc
[[ $WM_PROJECT_VERSION == 14 ]]
root=${1:?Pass an absolute path to an existing workspace directory}
[[ $root == /* && -d $root ]]
mesh="$root/spiral-pipe-mesh-v002"
case="$root/spiral-pipe-v001"

admit() {
    [[ $(awk '/^MemAvailable:/ {print $2}' /proc/meminfo) -ge 4194304 ]]
    [[ $(df -Pk "$root" | awk 'NR==2 {print $4}') -ge 10485760 ]]
    if ps -eo comm= | awk '$0 ~ /^(foamRun|foamMultiRun|mpirun|mpiexec|reconstructPar|checkMesh|blockMesh|extrudeMesh|snappyHexMesh|potentialFoam|foamPostProcess|foamToVTK)$/ {busy=1} END {exit !busy}'; then
        echo "A conflicting solver or mesh process is active" >&2
        exit 3
    fi
}
admit
(
    cd "$mesh"
    sha256sum -c mesh-sha256.txt
    grep -Fq 'Mesh OK.' check-mesh.log
    grep -Fq 'Number of regions: 1 (OK).' check-mesh.log
)
mkdir "$case"
cd "$case"
mkdir 0 system
cp -a "$mesh/constant" "$mesh/VTK" .
cp "$mesh/"{check-mesh.log,extrude.log,mesh.log,patch-export.log,mesh-sha256.txt} .
for file in boundary faces neighbour owner points; do
    cmp "constant/polyMesh/$file" "$mesh/constant/polyMesh/$file"
done
cp /opt/openfoam14/tutorials/incompressibleFluid/venturiTube/system/{fvSchemes,fvSolution} system/
cp /opt/openfoam14/tutorials/incompressibleFluid/venturiTube/constant/momentumTransport constant/
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
cat > constant/physicalProperties <<'EOF'
FoamFile { format ascii; class dictionary; object physicalProperties; }
viscosityModel constant;
nu 0.0025;
EOF
cat > 0/U <<'EOF'
FoamFile { format ascii; class volVectorField; location "0"; object U; }
dimensions [velocity];
internalField uniform (0 0 0);
boundaryField
{
    inlet
    {
        type flowRateInletVelocity;
        volumetricFlowRate constant 0.04;
        value uniform (0 0 0);
    }
    outlet
    {
        type pressureInletOutletVelocity;
        value uniform (0 0 0);
    }
    walls { type noSlip; }
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
    walls { type zeroGradient; }
}
EOF
cat > system/functions <<'EOF'
FoamFile { format ascii; class dictionary; object functions; }
#includeFunc patchFlowRate(patch=inlet)
#includeFunc patchFlowRate(patch=outlet)
#includeFunc patchAverage(patch=inlet, fields=(p))
#includeFunc patchAverage(patch=outlet, fields=(p))
EOF
grep -Eq '^simulationType[[:space:]]+laminar;' constant/momentumTransport
grep -Fq 'default         steadyState;' system/fvSchemes
set -o noclobber
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
foamPostProcess -func writeCellVolumes -time 0 > volumes.log 2>&1
test -s 0/C
test -s 0/Vc
find 0 constant system VTK -type f -print0 | sort -z | xargs -0 sha256sum > input-sha256.txt
sha256sum -c input-sha256.txt > inputs-verified.log
admit
foamRun > solver.log 2>&1
grep -qx End solver.log
