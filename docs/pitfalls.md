# Pitfalls

Learned on OpenFOAM v2412 with `buoyantSimpleFoam` (34 M cells, 9.3 M wall faces, 32 cores) and on the
example.

## `h` in 0/T must match `--sol-air H`

`Tsolair = Tair + qSensible/H` is written by the tool; `valueFraction` uses `h` in `0/T`. If they differ,
the energy balance is wrong without any error message.

## Smooth wall functions on buildings

`nutkWallFunction` gave `nut = 0` at low y+ on building walls, which killed wall heat transport. Use the
rough `atmNutkWallFunction` / `atmEpsilonWallFunction` on `Buildings` too (z0 0.05 m), with
`compressible::alphatWallFunction` (Prt 0.85) and `boundNut true`.

## `prghPressure` with a uniform p

On the outer boundary, `prghPressure` with a uniform `p` imposed no hydrostatic gradient and produced a
ring of very high velocities along the boundary. Use `fixedValue` for `p_rgh` at the outlet (with `p`
calculated) and `fixedFluxPressure` at the walls and inlet.

## z0 lists not split by decomposePar (v2412)

`decomposePar` splits `refValue` of the T BC correctly but not the `z0` lists of the atm wall functions:
woven pattern in `alphat`, `nut`, wall T and wall heat flux on the ground. Run `tools/decompose_z0.py`
of the z0 repository after `decomposePar`:

    python3 <z0-repo>/tools/decompose_z0.py . Terrain --list "solarIllumination/T_solar/{patch}.z0"

`tools/check_wall_lists.py` (same repository) shows the lists the solver actually used, face by face,
next to `illumination.vtp`.

## Lists belong to one mesh

`T_solar/*` are in polyMesh face order. After remeshing, rerun the tool; an old list either fails on the
face count or, worse, has the right count from a similar mesh and silently puts values on wrong faces.

## Which mesh is read

The tool takes the latest `<time>/polyMesh` (snappyHexMesh without `-overwrite` writes `1/`, `2/`), else
`constant/polyMesh`. Force one with `--time`. OpenFOAM itself reads `constant/polyMesh` unless told
otherwise: make sure both use the same mesh.

## Wall results depend on the first-cell height

Wall temperature and heat flux go through the wall functions (`alphat ∝ y/ln((y+z0)/z0)`,
`epsilon ∝ 1/y`), so the edge of a refinement box shows up in them. Interpret wall results where the
first cell height is uniform, or refine the ground to one level.

## RAM

The hit table is `faces × samples × 4 B`: 9.3 M faces × 256 samples ≈ 10 GB. Lower `--samples` or
`--chunk-rays` if memory is short; `--chunk-rays` limits the rays per batch, not the table.

## Shell variables in Allrun

`$samples` is defined inside `Allrun.sh`. Pasting the `solar_illumination.py` line alone into a shell
gives `--samples: expected one argument`.

## open3d import error on WSL

`import open3d` fails on a fresh Ubuntu/WSL because libusb is missing → `sudo apt install libusb-1.0-0`.

## Convergence

Residuals alone are not enough in buoyant runs. Watch a near-ground probe and the `wallHeatFlux`
integral; on the large case they needed 1000–1500 iterations.
