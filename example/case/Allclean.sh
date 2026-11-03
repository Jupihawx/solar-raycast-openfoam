#!/bin/bash
# Remove everything Allrun.sh generated (keeps geo/Buildings.stl and geo/site.json: delete them to re-download)
cd "$(dirname "$0")"
rm -rf constant/polyMesh constant/triSurface constant/extendedFeatureEdgeMesh constant/boundaryData \
       constant/wallLists processor* [1-9]* postProcessing logs dynamicCode solarIllumination \
       geo/*.png geo/*.csv *.foam check_*.vtp
