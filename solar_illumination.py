#!/usr/bin/env python3
"""
solar_illumination.py — per-face solar irradiance on OpenFOAM wall patches.

Pipeline
--------
1. Read wall patches from an OpenFOAM case (PyVista / vtkOpenFOAMReader).
2. Build an Embree ray-casting scene (Open3D) from the triangulated walls.
3. Direct pass:  one shadow ray per face toward the sun  -> q_direct.
4. Hemisphere pass: N cosine-weighted rays per face, cast ONCE; every ray
   records whether it escaped (sees sky) or which wall face it hit.
5. Sky view factor -> q_sky.  Diffuse bounces 1 and 2 are then pure table
   lookups on the cached hit table (no new rays).
6. Output: illumination.vtp for ParaView, and per-face lists in polyMesh face
   order for #include in 0/*: T_solar/<patch>.Tsolair (sol-air temperature,
   with --sol-air) and T_solar/<patch>.z0 (with a --material-map holding z0).

Conventions
-----------
* z is up.
* Sun azimuth is compass style: 0 deg = +y (north), 90 deg = +x (east).
  Elevation is degrees above the horizon.
* DNI = direct normal irradiance [W/m2], DHI = diffuse horizontal irradiance.
* albedo = diffuse reflectivity of a patch; absorbed flux = (1-albedo)*incident.
* OpenFOAM boundary-face normals point OUT of the fluid domain; the tool
  flips them to get the illuminated (fluid-facing) side.

Example
-------
python solar_illumination.py <case> --sun-azimuth 135 --sun-elevation 35 \
    --dni 800 --dhi 120 --albedo Terrain=0.2 --albedo Buildings=0.35 --sol-air 10 298
"""

import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pyvista as pv
import open3d as o3d
from scipy.spatial import cKDTree

INVALID_ID = 0xFFFFFFFF  # Open3D "no hit" primitive id


# --------------------------------------------------------------------------- #
#  Case / geometry loading
# --------------------------------------------------------------------------- #

def find_foam_file(case_dir: Path) -> Path:
    if case_dir.is_file() and case_dir.suffix == ".foam":
        return case_dir
    hits = sorted(case_dir.glob("*.foam"))
    if hits:
        return hits[0]
    stub = case_dir / f"{case_dir.name}.foam"
    stub.touch()
    return stub


def read_foam_list(path: Path, ncomp=1) -> np.ndarray:
    """ASCII OpenFOAM list (scalarField / vectorField, e.g. boundaryData) -> array."""
    text = path.read_text()
    body = text[text.index("(", text.index("}") if "FoamFile" in text else 0) + 1:text.rindex(")")]
    vals = np.array(body.replace("(", " ").replace(")", " ").split(), float)
    return vals.reshape(-1, ncomp) if ncomp > 1 else vals


def mesh_time(case_dir: Path, requested):
    """Time folder holding the mesh: --time if given, else the latest <time>/polyMesh
    (snappyHexMesh without -overwrite), else constant. Returns (time value or None, polyMesh dir)."""
    if requested is not None:
        d = case_dir / requested / "polyMesh"
        return (float(requested), d) if d.is_dir() else (None, case_dir / "constant" / "polyMesh")
    times = []
    for d in case_dir.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        if (d / "polyMesh" / "boundary").exists():
            times.append((t, d / "polyMesh"))
    return max(times) if times else (None, case_dir / "constant" / "polyMesh")


def parse_boundary_patches(poly_dir: Path) -> dict:
    """Return {patchName: patchType} from <polyMesh>/boundary
    (its dictionary part is ASCII even in binary cases)."""
    bfile = poly_dir / "boundary"
    text = bfile.read_text(errors="ignore")
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//.*", "", text)
    patches = {}
    for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", text):
        name, body = m.group(1), m.group(2)
        if name == "FoamFile":
            continue
        t = re.search(r"\btype\s+(\w+)\s*;", body)
        if t:
            patches[name] = t.group(1)
    return patches


def load_wall_surface(foam_file: Path, patch_names: list, time_value=None) -> pv.PolyData:
    """Read the requested patches and return one merged PolyData whose cell
    order is patch-by-patch. Adds cell arrays: patchId, gfid."""
    reader = pv.OpenFOAMReader(str(foam_file))
    if time_value is not None:
        reader.set_active_time_value(time_value)
    reader.disable_all_cell_arrays()    # geometry only: never parse 0/ (#include, #calc, mappedFile)
    available = {n.split("/", 1)[1] for n in reader.patch_array_names
                 if n.startswith("patch/")}
    missing = [p for p in patch_names if p not in available]
    if missing:
        sys.exit(f"Patches not found in case: {missing}. Available: {sorted(available)}")
    reader.disable_all_patch_arrays()
    for p in patch_names:
        reader.enable_patch_array(f"patch/{p}")
    reader.cell_to_point_creation = False
    bnd = reader.read()["boundary"]

    surfs = []
    for i, p in enumerate(patch_names):
        try:
            s = bnd[p].extract_surface(algorithm=None)   # newer pyvista
        except TypeError:
            s = bnd[p].extract_surface()                 # older pyvista: no 'algorithm'
        s.cell_data.clear()
        s.point_data.clear()
        s.cell_data["patchId"] = np.full(s.n_cells, i, np.int32)
        surfs.append(s)
    merged = surfs[0] if len(surfs) == 1 else pv.merge(surfs, merge_points=False)
    # gfid must be assigned after merging (merge may reorder cells; cell data
    # stays aligned, so patchId is safe, but downstream indexing uses position)
    merged.cell_data["gfid"] = np.arange(merged.n_cells, dtype=np.int64)
    return merged


# --------------------------------------------------------------------------- #
#  Sun vector
# --------------------------------------------------------------------------- #

def sun_vector(args) -> np.ndarray:
    """Unit vector pointing from the scene TOWARD the sun."""
    if args.sun_dir is not None:
        v = np.asarray(args.sun_dir, float)
        v /= np.linalg.norm(v)
        if v[2] < 0:  # user probably gave the light travel direction
            v = -v
            print("note: --sun-dir pointed below the horizon; flipped it to point at the sun.")
        return v
    az = np.radians(args.sun_azimuth)
    el = np.radians(args.sun_elevation)
    return np.array([np.sin(az) * np.cos(el),
                     np.cos(az) * np.cos(el),
                     np.sin(el)])


# --------------------------------------------------------------------------- #
#  Ray casting passes
# --------------------------------------------------------------------------- #

def build_scene(mesh: pv.PolyData):
    """Triangulate and load into Embree. Returns (scene, tri_to_face)."""
    tri = mesh.triangulate()
    tri_to_face = np.asarray(tri.cell_data["gfid"], np.int64)
    scene = o3d.t.geometry.RaycastingScene()
    scene.add_triangles(
        o3d.core.Tensor(np.asarray(tri.points, np.float32)),
        o3d.core.Tensor(tri.regular_faces.astype(np.uint32)),
    )
    return scene, tri_to_face


def direct_pass(scene, origins, normals, sun, dni):
    cos = normals @ sun
    lit_candidate = cos > 1e-6
    idx = np.flatnonzero(lit_candidate)
    occluded = np.zeros(len(origins), bool)
    if idx.size:
        rays = np.empty((idx.size, 6), np.float32)
        rays[:, :3] = origins[idx]
        rays[:, 3:] = sun.astype(np.float32)
        occluded[idx] = scene.test_occlusions(o3d.core.Tensor(rays)).numpy()
    sunlit = lit_candidate & ~occluded
    q_direct = np.where(sunlit, dni * np.maximum(cos, 0.0), 0.0)
    return q_direct.astype(np.float32), sunlit


def cosine_hemisphere_dirs(normals, n_samples, rng):
    """(F,3) normals -> (F,S,3) cosine-weighted directions about each normal."""
    F = len(normals)
    u1 = rng.random((F, n_samples), np.float32)
    u2 = rng.random((F, n_samples), np.float32)
    r = np.sqrt(u1)
    phi = (2.0 * np.pi) * u2
    lx, ly = r * np.cos(phi), r * np.sin(phi)
    lz = np.sqrt(np.maximum(1.0 - u1, 0.0))
    # per-face tangent frame
    n = normals
    a = np.where(np.abs(n[:, 2:3]) < 0.9, [[0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]])
    t1 = np.cross(a, n)
    t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    t2 = np.cross(n, t1)
    return (lx[..., None] * t1[:, None, :]
            + ly[..., None] * t2[:, None, :]
            + lz[..., None] * n[:, None, :]).astype(np.float32)


def hemisphere_pass(scene, origins, normals, tri_to_face, n_samples, chunk_rays, seed):
    """Cast S rays per face once. Returns:
       hit_table (F,S) int32: gfid of the face hit, -1 sky escape, -2 escape
       below horizon (dark);  and sky view factor (F,)."""
    F = len(origins)
    hit_table = np.empty((F, n_samples), np.int32)
    rng = np.random.default_rng(seed)
    faces_per_chunk = max(1, chunk_rays // n_samples)
    t0 = time.time()
    for start in range(0, F, faces_per_chunk):
        end = min(start + faces_per_chunk, F)
        dirs = cosine_hemisphere_dirs(normals[start:end], n_samples, rng)
        nrays = (end - start) * n_samples
        rays = np.empty((nrays, 6), np.float32)
        rays[:, :3] = np.repeat(origins[start:end], n_samples, axis=0)
        rays[:, 3:] = dirs.reshape(-1, 3)
        ans = scene.cast_rays(o3d.core.Tensor(rays))
        prim = ans["primitive_ids"].numpy()
        escaped = prim == INVALID_ID
        gfid = np.where(escaped, -1, tri_to_face[np.where(escaped, 0, prim)])
        gfid = gfid.astype(np.int32)
        up = rays[:, 5] > 0.0
        gfid[escaped & ~up] = -2  # escaped below horizon: neither sky nor surface
        hit_table[start:end] = gfid.reshape(end - start, n_samples)
        done = end * n_samples
        print(f"  hemisphere pass: {done/1e6:.1f} M rays "
              f"({end}/{F} faces, {time.time()-t0:.1f}s)", end="\r")
    print()
    sky_vf = (hit_table == -1).mean(axis=1).astype(np.float32)
    return hit_table, sky_vf


def bounce_pass(hit_table, albedo, e_prev):
    """One diffuse bounce: mean over rays of (albedo * E_prev) of the hit face."""
    radiosity = (albedo * e_prev).astype(np.float32)
    lut = np.concatenate([radiosity, [0.0, 0.0]]).astype(np.float32)  # -1,-2 -> 0
    return lut[hit_table].mean(axis=1)


# --------------------------------------------------------------------------- #
#  OpenFOAM output
# --------------------------------------------------------------------------- #

def write_face_list(path: Path, vals):
    """Bare per-face list, in polyMesh face order, for  <entry> #include "<path>";  in 0/*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(f"nonuniform List<scalar> \n{len(vals)}\n(\n")
        f.write("\n".join(f"{v:.6g}" for v in vals))
        f.write("\n)\n")


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("case", type=Path, help="OpenFOAM case directory (or .foam file)")
    ap.add_argument("--patches", nargs="+", default=None,
                    help="wall patches to illuminate (default: all type-wall patches)")
    sun = ap.add_argument_group("sun")
    sun.add_argument("--sun-azimuth", type=float, default=135.0,
                     help="compass azimuth of the sun [deg], 0=N(+y), 90=E(+x)")
    sun.add_argument("--sun-elevation", type=float, default=35.0,
                     help="sun elevation above horizon [deg]")
    sun.add_argument("--sun-dir", nargs=3, type=float, metavar=("X", "Y", "Z"),
                     help="explicit vector toward the sun (overrides azimuth/elevation)")
    rad = ap.add_argument_group("radiation")
    rad.add_argument("--dni", type=float, default=800.0, help="direct normal irradiance [W/m2]")
    rad.add_argument("--dhi", type=float, default=120.0, help="diffuse horizontal irradiance [W/m2]")
    rad.add_argument("--albedo", action="append", default=[], metavar="PATCH=VALUE",
                     help="per-patch albedo, repeatable")
    rad.add_argument("--default-albedo", type=float, default=0.2)
    rad.add_argument("--sensible-fraction", action="append", default=[], metavar="PATCH=VALUE",
                     help="fraction of ABSORBED solar delivered to the air as sensible heat "
                          "(the rest goes to evaporation / ground storage); per-patch, "
                          "repeatable. asphalt~1.0, grass~0.5, water~0.2")
    rad.add_argument("--default-sensible-fraction", type=float, default=1.0)
    rad.add_argument("--material-map", action="append", default=[], metavar="PATCH=DIR",
                     help="per-face albedo/sensibleFraction (and z0 if present) from a boundaryData "
                          "folder (points + 0/albedo, 0/sensibleFraction, 0/z0, e.g. "
                          "constant/boundaryData/Terrain from z0-worldcover-openfoam); nearest point; "
                          "overrides --albedo/--sensible-fraction on that patch")
    rad.add_argument("--bounces", type=int, default=2, choices=(0, 1, 2))
    mc = ap.add_argument_group("sampling")
    mc.add_argument("--samples", type=int, default=256, help="hemisphere rays per face")
    mc.add_argument("--chunk-rays", type=int, default=4_000_000, help="rays per Embree batch")
    mc.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-flip-normals", action="store_true",
                    help="do NOT flip patch normals (default flips, since OpenFOAM "
                         "boundary normals point out of the fluid)")
    ap.add_argument("--time", default=None,
                    help="time folder holding the mesh (default: latest <time>/polyMesh, else constant)")
    rad.add_argument("--sol-air", nargs=2, type=float, metavar=("H", "TAIR"), default=None,
                     help="write T_solar/<patch>.Tsolair = TAIR + qSensible/H [K], the sol-air "
                          "temperature for a mixed BC (H [W/m2/K]: surface loss coefficient, "
                          "longwave + storage; TAIR [K]: air temperature)")
    ap.add_argument("--out", type=Path, default=None,
                    help="output directory (default: <case>/solarIllumination)")
    args = ap.parse_args()

    args.case = args.case.resolve()  # so "." works when run from inside a case
    case_dir = args.case if args.case.is_dir() else args.case.parent
    out_dir = args.out or case_dir / "solarIllumination"
    out_dir.mkdir(parents=True, exist_ok=True)

    time_value, poly_dir = mesh_time(case_dir, args.time)
    print(f"mesh: {poly_dir.relative_to(case_dir)}")
    boundary = parse_boundary_patches(poly_dir)
    patch_names = args.patches or [n for n, t in boundary.items() if t == "wall"]
    if not patch_names:
        sys.exit("No wall patches found; specify --patches explicitly.")
    print(f"patches: {patch_names}")

    albedo_map = {p: args.default_albedo for p in patch_names}
    for spec in args.albedo:
        name, _, val = spec.partition("=")
        if name not in albedo_map:
            sys.exit(f"--albedo {spec}: patch '{name}' is not in {patch_names}")
        albedo_map[name] = float(val)
    print(f"albedo: {albedo_map}")

    sens_map = {p: args.default_sensible_fraction for p in patch_names}
    for spec in args.sensible_fraction:
        name, _, val = spec.partition("=")
        if name not in sens_map:
            sys.exit(f"--sensible-fraction {spec}: patch '{name}' is not in {patch_names}")
        sens_map[name] = float(val)
    print(f"sensible fraction: {sens_map}")

    sun = sun_vector(args)
    print(f"sun vector (toward sun): {np.round(sun, 4)}")
    if sun[2] <= 0:
        print("warning: sun at or below horizon; direct pass will be zero.")

    # ---- geometry ----
    t0 = time.time()
    mesh = load_wall_surface(find_foam_file(args.case), patch_names, time_value)
    mesh = mesh.compute_normals(cell_normals=True, point_normals=False,
                                split_vertices=False, consistent_normals=False,
                                auto_orient_normals=False)
    normals = np.asarray(mesh.cell_data["Normals"], np.float64)
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    if not args.no_flip_normals:
        normals = -normals  # fluid-facing side
    centroids = mesh.cell_centers().points.astype(np.float64)
    areas = np.asarray(mesh.compute_cell_sizes(length=False, volume=False)
                       .cell_data["Area"], np.float64)
    patch_id = np.asarray(mesh.cell_data["patchId"])
    albedo = np.array([albedo_map[p] for p in patch_names], np.float32)[patch_id]
    sensible = np.array([sens_map[p] for p in patch_names], np.float32)[patch_id]
    z0_faces = {}
    for spec in args.material_map:
        name, _, folder = spec.partition("=")
        if name not in patch_names:
            sys.exit(f"--material-map {spec}: patch '{name}' is not in {patch_names}")
        pts = read_foam_list(Path(folder) / "points", 3)
        sel = patch_id == patch_names.index(name)
        near = cKDTree(pts).query(centroids[sel])[1]
        albedo[sel] = read_foam_list(Path(folder) / "0" / "albedo")[near]
        sensible[sel] = read_foam_list(Path(folder) / "0" / "sensibleFraction")[near]
        if (Path(folder) / "0" / "z0").exists():   # per-face z0 for the wall functions
            z0_faces[name] = read_foam_list(Path(folder) / "0" / "z0")[near]
        print(f"material map on {name}: {len(pts)} points, albedo {albedo[sel].min():.2f}-"
              f"{albedo[sel].max():.2f}, sensible {sensible[sel].min():.2f}-{sensible[sel].max():.2f}")
    F = mesh.n_cells
    print(f"loaded {F} wall faces in {time.time()-t0:.1f}s "
          f"(mean upward-normal fraction: {(normals[:, 2] > 0).mean():.2f})")

    scene, tri_to_face = build_scene(mesh)
    eps = 0.02 * np.sqrt(areas) + 1e-6
    origins = (centroids + eps[:, None] * normals).astype(np.float32)
    normals32 = normals.astype(np.float32)

    # ---- direct sun ----
    t0 = time.time()
    q_direct, sunlit = direct_pass(scene, origins, normals, sun, args.dni)
    print(f"direct pass: {sunlit.mean()*100:.1f}% of faces sunlit ({time.time()-t0:.1f}s)")

    # ---- hemisphere: sky + bounce visibility, cast once ----
    hit_table, sky_vf = hemisphere_pass(scene, origins, normals32, tri_to_face,
                                        args.samples, args.chunk_rays, args.seed)
    q_sky = (args.dhi * sky_vf).astype(np.float32)

    # ---- bounces (table lookups only) ----
    e0 = q_direct + q_sky
    q_b1 = bounce_pass(hit_table, albedo, e0) if args.bounces >= 1 else np.zeros(F, np.float32)
    q_b2 = bounce_pass(hit_table, albedo, q_b1) if args.bounces >= 2 else np.zeros(F, np.float32)
    q_inc = e0 + q_b1 + q_b2
    q_abs = (1.0 - albedo) * q_inc          # physical absorbed shortwave
    q_sens = sensible * q_abs               # part delivered to the air (BC value)

    ghi = args.dni * max(sun[2], 0.0) + args.dhi
    print(f"\nresults (W/m2):        mean      max")
    for name, q in (("direct", q_direct), ("sky", q_sky), ("bounce1", q_b1),
                    ("bounce2", q_b2), ("incident", q_inc), ("absorbed", q_abs),
                    ("sensible", q_sens)):
        print(f"  {name:10s} {q.mean():10.1f} {q.max():9.1f}")
    print(f"  (reference GHI on open flat ground: {ghi:.1f})")

    # ---- outputs ----
    for name, arr in (("qDirect", q_direct), ("qSky", q_sky), ("qBounce1", q_b1),
                      ("qBounce2", q_b2), ("qIncident", q_inc), ("qAbsorbed", q_abs),
                      ("qSensible", q_sens), ("skyViewFactor", sky_vf),
                      ("albedo", albedo), ("sensibleFraction", sensible)):
        mesh.cell_data[name] = arr
    mesh.cell_data["sunlit"] = sunlit.astype(np.uint8)
    mesh.field_data["patchNames"] = np.array(patch_names)
    vtp = out_dir / "illumination.vtp"
    mesh.save(vtp)
    print(f"\nwrote {vtp}")

    for i, p in enumerate(patch_names):
        sel = patch_id == i
        if args.sol_air:                                 # sol-air temperature, polyMesh face order
            h_sa, t_air = args.sol_air
            write_face_list(out_dir / "T_solar" / f"{p}.Tsolair", t_air + q_sens[sel] / h_sa)
            print(f"wrote {out_dir / 'T_solar' / (p + '.Tsolair')}")
        if p in z0_faces:                                # z0 from --material-map, same face order
            write_face_list(out_dir / "T_solar" / f"{p}.z0", z0_faces[p])
            print(f"wrote {out_dir / 'T_solar' / (p + '.z0')}")


if __name__ == "__main__":
    main()
