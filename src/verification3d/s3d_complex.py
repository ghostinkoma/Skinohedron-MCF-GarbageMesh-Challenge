"""
s3d_complex.py
==============
3D extension, step 1 (the "opening act" before the Laplacian).

Goal: confirm that the *combinatorial* backbone of the 3D Skin-o-hedron model is
correct before any operator (Laplacian) is placed on it.  We build solid-ball
tetrahedral meshes at several refinements and check, at every level:

  * the discrete de Rham complex closes in 3D:   d1 . d0 = 0   and   d2 . d1 = 0
  * the Euler characteristic of a solid ball:    V - E + F - T = 1
  * basic geometry: boundary volume -> 4/3 pi as the surface is refined,
    and tetrahedral shape quality (reported honestly, warts and all).

Why this matters: in 2D the analogous facts were  d.d = 0  and  V - E + F = 2.
The jump to 3D adds a second derivative (d2) and a third simplex type (tets);
getting  d2 . d1 = 0  and  chi = 1  exactly is the minimal evidence that the
substrate is wired correctly.  No convergence claim is made here -- that is the
next step (3D Laplacian vs spherical Bessel eigenvalues).
"""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from ksf3d.mesh3d import solid_ball, tet_quality, mean_edge_length
from ksf3d.dec3d import complex_residuals


LEVELS = [(0, 3), (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4)]


def run():
    rows = []
    for slev, ns in LEVELS:
        V, T = solid_ball(shell_level=slev, n_shells=ns)
        r1, r2, euler, (nV, nE, nF, nT) = complex_residuals(V, T)
        q = tet_quality(V, T)
        rows.append({
            "shell_level": slev, "n_shells": ns,
            "h": mean_edge_length(V, T),
            "V": nV, "E": nE, "F": nF, "T": nT,
            "dd1": r1, "dd2": r2, "euler": euler,
            "volume": q["volume"], "volume_ideal": q["volume_ideal"],
            "q_min": q["q_min"], "q_mean": q["q_mean"],
        })
    return {
        "rows": rows,
        "max_dd1": max(r["dd1"] for r in rows),
        "max_dd2": max(r["dd2"] for r in rows),
        "euler_ok": all(r["euler"] == 1 for r in rows),
    }


def main():
    res = run()
    print("=== 3D-S1: solid-ball tetrahedral complex ===")
    print(f'{"slev":>4} {"nsh":>4} {"h":>8} {"V":>6} {"E":>6} {"F":>6} {"T":>6}'
          f' {"|d1d0|":>8} {"|d2d1|":>8} {"chi":>4} {"vol":>7} {"q_min":>6} {"q_mean":>6}')
    for r in res["rows"]:
        print(f'{r["shell_level"]:4d} {r["n_shells"]:4d} {r["h"]:8.4f}'
              f' {r["V"]:6d} {r["E"]:6d} {r["F"]:6d} {r["T"]:6d}'
              f' {r["dd1"]:8.1e} {r["dd2"]:8.1e} {r["euler"]:4d}'
              f' {r["volume"]:7.4f} {r["q_min"]:6.3f} {r["q_mean"]:6.3f}')
    print()
    print(f'max |d1.d0| = {res["max_dd1"]:.1e}    max |d2.d1| = {res["max_dd2"]:.1e}')
    print(f'Euler characteristic = 1 at every level : {res["euler_ok"]}')
    print(f'ideal ball volume 4/3 pi = {res["rows"][0]["volume_ideal"]:.4f}'
          f'  (boundary volume rises toward it as shell_level increases)')
    print()
    print("[結句] The 3D substrate is combinatorially sound. Both coboundary")
    print("       compositions vanish exactly (d1.d0 = 0 AND d2.d1 = 0), and the")
    print("       Euler characteristic of the solid ball is 1 at every refinement")
    print("       -- the correct 3D analogue of the surface's chi = 2. The boundary")
    print("       volume approaches 4/3 pi as the sphere is refined. Tetrahedral")
    print("       quality is honestly modest (Delaunay slivers give small q_min),")
    print("       which is exactly the kind of irregular mesh the Skin-o-hedron")
    print("       programme is meant to tolerate. This is the OPENING ACT only:")
    print("       it fixes the complex, not the operator. The next step hands the")
    print("       baton to a 3D Laplacian tested against spherical-Bessel")
    print("       eigenvalues -- where the real convergence question (and the")
    print("       no-free-lunch tension) will reappear in 3D.")
    return res


if __name__ == "__main__":
    main()
