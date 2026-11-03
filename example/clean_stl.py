#!/usr/bin/env python3
"""Clean a building STL exported from Blender + Blosm for snappyHexMesh.

    python clean_stl.py <in.stl> <out.stl> [--min-area 25]

Removes, per solid (connected part):
  - flat sheets (no height),
  - solids whose footprint (area of the upward faces) is below --min-area m2: sheds, small parts that
    are smaller than the mesh resolution and only make slivers,
  - parts that start above the ground with nothing under them (unsupported floating parts).
    Stacked building parts (an upper part resting on a lower roof) are kept.
Then snaps vertices to 1 mm and drops triangles that collapse or are duplicated: OpenFOAM merges close
points when reading an STL, and one such triangle made surfaceFeatureExtract stop with a floating
point exception.

The example's geo/Buildings.stl was made with this script (--min-area 25) from a Blosm export of
OSM buildings (lon 4.36745-4.38655, lat 51.99439-52.00726), after removing some buildings by hand.
"""
import argparse

import numpy as np
import trimesh

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("inp")
ap.add_argument("out")
ap.add_argument("--min-area", type=float, default=25.0, help="smallest footprint kept, m2")
a = ap.parse_args()

m = trimesh.load(a.inp, force="mesh")
b = m.split(only_watertight=False)
lo = np.array([x.bounds[0] for x in b])
hi = np.array([x.bounds[1] for x in b])
fp = np.array([x.area_faces[x.face_normals[:, 2] > 0.5].sum() for x in b])
supported = np.ones(len(b), bool)
for i in np.flatnonzero(lo[:, 2] > 0.1):          # starts above the ground: anything underneath?
    ov = (lo[:, 0] < hi[i, 0]) & (hi[:, 0] > lo[i, 0]) & (lo[:, 1] < hi[i, 1]) & (hi[:, 1] > lo[i, 1])
    ov[i] = False
    supported[i] = (ov & (hi[:, 2] >= lo[i, 2] - 0.5) & (lo[:, 2] < lo[i, 2])).any()
flat = hi[:, 2] - lo[:, 2] < 0.01
small = fp < a.min_area
keep = ~flat & ~small & supported
print(f"{len(b)} solids: flat {flat.sum()}, footprint < {a.min_area:g} m2 {small.sum()}, "
      f"unsupported {(~supported).sum()} -> kept {keep.sum()}")

out = trimesh.util.concatenate([b[i] for i in np.flatnonzero(keep)])
n0 = len(out.faces)
out.vertices = np.round(out.vertices, 3)
out.merge_vertices()
out.update_faces(out.nondegenerate_faces(height=1e-3))
out.update_faces(out.unique_faces())
out.remove_unreferenced_vertices()
out.export(a.out)
print(f"repair: removed {n0 - len(out.faces)} degenerate/duplicate triangles; "
      f"{len(out.faces)} triangles, height 0-{out.bounds[1, 2]:.1f} m -> {a.out}")
