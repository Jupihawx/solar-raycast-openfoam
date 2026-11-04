# Example: TU Delft campus (solar ray casting + optional buoyant run)

Same site as [z0-worldcover-openfoam](https://github.com/Jupihawx/z0-worldcover-openfoam)'s example: the
TU Delft campus and the housing to its east, ~1.3 km × 1.3 km of buildings, origin 52.000825 N, 4.377 E.

| | |
|---|---|
| Geometry | `geo/Buildings.stl`, shipped: OSM buildings, cleaned (small sheds and broken triangles removed). 406 solids, tallest 90 m (EWI) |
| Domain | box aligned with the wind: buildings within 1000 m of the origin, plus 300 m upstream and to the sides, 900 m downstream, 300 m high. **Smaller than the COST 732 guideline** to keep the example short |
| Mesh | blockMesh 20 m cells, rotated to `dirFrom` with `transformPoints`. snappyHexMesh: the whole ground at 10 m (no refinement boxes), 5 m within 10 m of buildings, 2.5 m on building surfaces. 855 k cells |
| Sun | 2023-06-21 14:00 UTC: azimuth 235.3°, elevation 51.5°; clear sky DNI 750, DHI 130 W/m² (illustrative) |
| Inflow | 5 m/s at 100 m from 225°, `z0In` 0.1 m, T2m 25 °C (illustrative) |
| Materials | Buildings albedo 0.30, sensible fraction 1.0. Ground: per face from WorldCover with `Z0_REPO`, else albedo 0.20, sensible fraction 0.7, z0 0.1 m |
| Solver (optional) | `buoyantSimpleFoam`, Boussinesq EoS, `buoyantKEpsilon`, rough atm wall functions, sol-air T BC (h = 10 W/m²K) |

## Step by step

Linux or WSL2 (Ubuntu). Commands start in the repository root.

**1. OpenFOAM (openfoam.com; tested with v2412 and v2606).** Install it from <https://develop.openfoam.com/Development/openfoam/-/wikis/precompiled>
(Ubuntu: `apt` package `openfoam2412-default`), then load it in every new shell:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
buoyantSimpleFoam -help > /dev/null && echo OK    # check
```

**2. Python environment.** open3d needs Python ≤ 3.12; uv fetches that version for you:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
~/.local/bin/uv venv --python 3.12 ~/venv-solar
~/.local/bin/uv pip install --python ~/venv-solar -r requirements.txt
sudo apt install libusb-1.0-0                    # WSL/Ubuntu: open3d import fails without it
~/venv-solar/bin/python -c "import open3d, pyvista; print('OK')"
```

**3. (Optional) the z0 repository**, for ground albedo, sensible fraction and z0 per face from ESA WorldCover:

```bash
git clone https://github.com/Jupihawx/z0-worldcover-openfoam ~/z0-worldcover-openfoam
```

Its Python packages are already in this repository's `requirements.txt`.

**4. Run the example** (with `Z0_REPO`, internet is needed on the first run for one 55 MB WorldCover tile):

```bash
cd example/case
PYTHON=~/venv-solar/bin/python bash Allrun.sh                                   # A: ray casting, uniform ground
PYTHON=~/venv-solar/bin/python Z0_REPO=~/z0-worldcover-openfoam bash Allrun.sh   # B: + ground from WorldCover
PYTHON=~/venv-solar/bin/python Z0_REPO=~/z0-worldcover-openfoam BUOYANT=1 NPROCS=4 bash Allrun.sh   # C: + buoyant run
```

Pick one; the mesh is reused if it exists. `Allrun.sh` does, with a log per step in `logs/`:

| # | Step | Writes |
|---|---|---|
| 1 | read the origin from `geo/site.json` | |
| 2 | `blockMesh`, `transformPoints` (box turned to `dirFrom`), `snappyHexMesh`, `checkMesh` | `constant/polyMesh` |
| 3 | B, C: `make_z0.py` of the z0 repository | `constant/boundaryData/Terrain` |
| 4 | `solar_illumination.py` with `--sol-air 10 298.15` | `solarIllumination/` (see below) |
| 5 | C: `decomposePar`, `decompose_z0.py`, `buoyantSimpleFoam` (1000 iterations), `reconstructPar` | `1000/`, `postProcessing/` |

**5. Look at the results**

```bash
touch case.foam          # ParaView: latest time, patches Terrain and Buildings
```

- `solarIllumination/illumination.vtp`: open directly in ParaView, colour by `qIncident`, `sunlit`,
  `skyViewFactor`, `albedo`.

See [Check](#check) below. `bash Allclean.sh` removes everything generated (keeps `geo/Buildings.stl` and `geo/site.json`).

To use your own site, see "Your own site" in the z0 repository's `example/README.md` (geometry,
`site.json`, domain size).

## Runtime (measured)

16-core WSL2 machine, 8 MPI ranks, OpenFOAM v2606, WorldCover tile cached, `Z0_REPO` and `BUOYANT=1` set:

| Step | Time |
|---|---|
| blockMesh + rotate + snappyHexMesh (serial) | ~2.5 min, 855 k cells |
| z0 map (`make_z0.py`) | ~3 s |
| Ray casting, 243 k wall faces × 256 samples (62 M rays) | ~20 s, hit table 0.25 GB of RAM; `illumination.vtp` 13 MB on disk |
| buoyantSimpleFoam, 1000 iterations | ~15 min (0.88 s per iteration) |
| **Total** | **~18 min** (1066 s); without `BUOYANT=1` about 3 min |

Result of that run: 77 % of wall faces sunlit, mean incident 520 W/m² over all walls (reference GHI on
open ground 717). Buoyant run stable: `limitTemperature` inactive (0 cells), max air temperature 311 K.
The upstream 10 m probe settles by ~400 iterations; the 100 m probe above the site still drifts by
1–3 % at 1000.

## Check

- `logs/solar.log`: the reference GHI (DNI·sin(el) + DHI) against `incident` on open ground; share of
  sunlit faces.
- `solarIllumination/illumination.vtp` in ParaView: shadows fall to the north-east (sun from the
  south-west); `skyViewFactor` drops in courtyards; with `Z0_REPO`, `albedo` on the ground follows land cover.
- Buoyant run: `limitTemperature` should report 0 cells once settled; probes in `postProcessing/probes/`;
  `wallHeatFlux` per write time.

## Licence of the geometry

`geo/Buildings.stl` is derived from OpenStreetMap: © OpenStreetMap contributors, available under the
Open Database Licence (ODbL 1.0), <https://www.openstreetmap.org/copyright>. It is a derived database
under the ODbL, separate from the MIT-licensed code; if you share or adapt it, keep this attribution and
the ODbL terms.