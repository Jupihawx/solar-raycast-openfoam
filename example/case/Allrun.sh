#!/bin/bash
# TU Delft example, solar ray casting (+ optional buoyant run), OpenFOAM (openfoam.com, tested v2412/v2606):
#   geo/Buildings.stl (OSM buildings) -> blockMesh + snappyHexMesh -> [z0/albedo map from z0-worldcover-openfoam]
#   -> solar_illumination.py -> [buoyantSimpleFoam with the sol-air BC]
#
# Usage:  bash Allrun.sh            (source OpenFOAM's etc/bashrc first)
#   PYTHON=/path/to/venv/bin/python   Python with requirements.txt installed (default: python3)
#   Z0_REPO=/path/to/z0-worldcover-openfoam   ground albedo, sensible fraction and z0 per face from
#                                     ESA WorldCover; without it: uniform ground (albedo 0.2, z0 0.1 m)
#   BUOYANT=1                         also run buoyantSimpleFoam (optional, ~10-30 min)
#   NPROCS=4                          MPI ranks for the solver
set -e
cd "$(dirname "$0")"
REPO=$(cd ../.. && pwd)
PY=${PYTHON:-python3}
nprocs=${NPROCS:-4}
samples=256        # hemisphere rays per face; RAM ~ wallFaces x samples x 4 B
mkdir -p logs

# Example day (see 0/include/ABLConditions): 2023-06-21 14:00 UTC, clear sky, illustrative values
SUN="--sun-azimuth 235.3 --sun-elevation 51.5 --dni 750 --dhi 130"
SOLAIR="--sol-air 10 298.15"       # H = h in the 0/T codedMixed code, TAIR = Tair in ABLConditions

# 1. Geometry: geo/Buildings.stl and its projection origin geo/site.json (shipped with the example)
LAT0=$($PY -c "import json; print(json.load(open('geo/site.json'))['lat0'])")
LON0=$($PY -c "import json; print(json.load(open('geo/site.json'))['lon0'])")

# 2. Mesh: box in the wind frame (blockMesh), rotated about z so x points downwind, + buildings (snappy)
#    flowDir = (-sin dirFrom, -cos dirFrom): rotation angle = atan2(-cos dirFrom, -sin dirFrom)
if [ ! -f constant/polyMesh/owner ]; then
    DIRFROM=$(foamDictionary -entry dirFrom -value 0/include/ABLConditions)
    ANGLE=$($PY -c "import math; d=math.radians($DIRFROM); print(round(math.degrees(math.atan2(-math.cos(d), -math.sin(d))), 4))")
    mkdir -p constant/triSurface
    cp geo/Buildings.stl constant/triSurface/
    surfaceFeatureExtract > logs/surfaceFeatureExtract.log
    blockMesh > logs/blockMesh.log
    transformPoints -rotate-angle "((0 0 1) $ANGLE)" > logs/transformPoints.log
    snappyHexMesh -overwrite > logs/snappyHexMesh.log
    checkMesh > logs/checkMesh.log || echo "checkMesh reported problems, see logs/checkMesh.log"
fi

# 3. Ground materials (optional): albedo, sensibleFraction, z0 from ESA WorldCover
MATERIALS="--albedo Terrain=0.20 --sensible-fraction Terrain=0.7"
if [ -n "$Z0_REPO" ]; then
    $PY "$Z0_REPO/make_z0.py" --case . --lat0 "$LAT0" --lon0 "$LON0" --dx 10 --radius 2400 | tee logs/make_z0.log
    MATERIALS="--material-map Terrain=constant/boundaryData/Terrain"
fi

# 4. Ray casting -> solarIllumination/illumination.vtp, T_solar/<patch>.Tsolair (+ Terrain.z0)
$PY "$REPO/solar_illumination.py" . --patches Terrain Buildings $SUN $MATERIALS \
    --albedo Buildings=0.30 --sensible-fraction Buildings=1.0 $SOLAIR --samples $samples \
    | tee logs/solar.log
# without the material map there is no Terrain.z0 list: the wall functions get a uniform value
[ -f solarIllumination/T_solar/Terrain.z0 ] || echo "uniform 0.1" > solarIllumination/T_solar/Terrain.z0

# 5. Buoyant run (optional)
[ "$BUOYANT" = 1 ] || { echo "Done (ray casting). BUOYANT=1 bash Allrun.sh to run buoyantSimpleFoam."; exit 0; }
sed -i "s/^numberOfSubdomains.*/numberOfSubdomains $nprocs;/" system/decomposeParDict
rm -rf processor*
decomposePar -force > logs/decomposePar.log
# decomposePar (v2412) copies the z0 list of the atm wall functions unsplit: write each processor's slice
if [ -n "$Z0_REPO" ]; then
    $PY "$Z0_REPO/tools/decompose_z0.py" . Terrain --list "solarIllumination/T_solar/{patch}.z0" \
        | tee logs/decompose_z0.log
fi
mpirun -np $nprocs buoyantSimpleFoam -parallel > logs/buoyantSimpleFoam.log
reconstructPar -latestTime > logs/reconstructPar.log
echo "Done. Open the case in ParaView (touch case.foam) and solarIllumination/illumination.vtp"
