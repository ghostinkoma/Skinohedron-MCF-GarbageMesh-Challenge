"""
run_all.py
==========
Runs every per-section verification, prints each section's conclusion (結句), and
writes results/results.json consumed by viewer/viewer.html.

The JSON carries:
  * convergence curves (regular + irregular, pointwise + spectral)  [from s6]
  * the per-section scalar diagnostics                              [s2,s3,s4,s7,s10]
  * mesh geometry + spherical-harmonic vertex fields for a few levels (for the
    GLSL sphere in the viewer).
"""
from __future__ import annotations
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import s2_layered_complex as s2
import s3_dec as s3
import s4_trace_free as s4
import s6_laplacian_consistency as s6
import s7_nondelaunay as s7
import s10_sparameter as s10
from ksf import mesh, sph


def _mesh_payload(levels=(2, 3, 4)):
    """Geometry + harmonic fields for the 3D viewer."""
    payload = []
    for lvl in levels:
        V, F = mesh.icosphere(lvl)
        fields = {name: sph.harmonic(V, name)[0].tolist()
                  for name in ("Y10", "Y2xy", "Y3xyz")}
        payload.append(dict(
            level=int(lvl),
            h=mesh.mean_edge_length(V, F),
            n=len(V),
            vertices=[[round(float(c), 6) for c in p] for p in V],
            faces=[[int(i) for i in f] for f in F],
            normals=[[round(float(c), 6) for c in p]
                     for p in mesh.vertex_normals(V, F)],
            fields=fields,
        ))
    return payload


def main():
    print("#" * 70)
    print("#  Kosaka Skin-o-hedron Model -- full verification run")
    print("#" * 70, "\n")

    sec2 = s2.run();   s2.main();  print()
    sec3 = s3.run();   s3.main();  print()
    sec4 = s4.run();   s4.main();  print()
    sec6 = s6.run();   s6.main();  print()
    sec7 = s7.run();   s7.main();  print()
    sec10 = s10.run(); s10.main(); print()

    results = {
        "version": "1.0.0",
        "s2_layered_complex": sec2,
        "s3_dec": sec3,
        "s4_trace_free": sec4,
        "s6_consistency": sec6,
        "s7_nondelaunay": sec7,
        "s10_sparameter": sec10,
        "viewer_mesh": _mesh_payload(),
    }
    out = os.path.join(HERE, "..", "results", "results.json")
    out = os.path.abspath(out)
    with open(out, "w") as fh:
        json.dump(results, fh)
    size = os.path.getsize(out) / 1024
    print("#" * 70)
    print(f"#  wrote {out}  ({size:.0f} KB)")
    print("#" * 70)


if __name__ == "__main__":
    main()
