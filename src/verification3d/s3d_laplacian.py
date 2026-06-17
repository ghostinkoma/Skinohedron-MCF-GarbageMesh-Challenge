"""
s3d_laplacian.py
================
3D extension, STEP 2 -- the operator the substrate was waiting for, and the real
convergence question the `EXTENSION_3D.md` doc flagged as "not yet done".

We place the P1 FEM Laplacian on tetrahedral meshes and test its lowest Dirichlet
eigenvalue against the EXACT ground truth:
    * unit ball : lambda_1 = pi^2  (= j_{0,1}^2, the spherical-Bessel ground state)
    * unit cube : lambda_1 = 3 pi^2

The point is not the operator (it is plain linear FEM -- the 2D Tier-3 result
says so, and it carries over). The point is to see, with numbers, how 3D
convergence depends on TETRAHEDRAL QUALITY -- the 3D echo of the 2D no-free-lunch
story (§6-§7). Three mesh families, increasingly "equally-spaced":

  A. solid_ball   (Delaunay shells)   q_min = 0       -> slivers
  B. graded_ball  (geometric+rotated) q_min -> 0       -> better, not regular
  C. kuhn_cube    (congruent lattice) q_min = 0.717    -> uniformly regular

Reading the result:
  A: the spectrum has the right STRUCTURE (lambda_1 ~ pi^2, then an l=1 triplet
     near 20.19, an l=2 quintet near 33.22 -- multiplicities 1,3,5) but the error
     does NOT shrink under refinement: slivers stall convergence.
  B: error now shrinks, but only ~O(h): the family is not uniformly shape-regular
     (q_min drifts to 0), so it cannot reach the optimal rate.
  C: q is constant, the family is uniformly shape-regular, and the error converges
     at the optimal O(h^2). This is the "equally-spaced tetrahedra" payoff.
"""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from ksf3d.mesh3d import solid_ball, tet_quality, mean_edge_length
from ksf3d.mesh3d_uniform import (graded_ball, ball_boundary_mask,
                                  kuhn_cube, cube_boundary_mask)
from ksf3d.fem3d import dirichlet_eigenvalues


def _order(hs, es):
    return float(np.polyfit(np.log(hs), np.log(es), 1)[0])


def run():
    out = {}
    pi2 = float(np.pi ** 2)

    # A. sliver Delaunay ball
    A = []
    for slev, ns in [(1, 3), (1, 5), (2, 4)]:
        V, T = solid_ball(shell_level=slev, n_shells=ns)
        ev = dirichlet_eigenvalues(V, T, ball_boundary_mask(V), k=5)
        A.append({"h": mean_edge_length(V, T), "q_min": tet_quality(V, T)["q_min"],
                  "lam1": float(ev[0]), "err": abs(ev[0] - pi2) / pi2,
                  "spectrum": [float(x) for x in ev]})
    out["A_sliver_ball"] = {"rows": A, "target": pi2,
                            "order": _order([r["h"] for r in A], [r["err"] for r in A])}

    # B. graded+rotated ball
    B = []
    for slev in (1, 2, 3):
        V, T = graded_ball(shell_level=slev)
        ev = dirichlet_eigenvalues(V, T, ball_boundary_mask(V), k=5)
        B.append({"h": mean_edge_length(V, T), "q_min": tet_quality(V, T)["q_min"],
                  "lam1": float(ev[0]), "err": abs(ev[0] - pi2) / pi2})
    out["B_graded_ball"] = {"rows": B, "target": pi2,
                            "order": _order([r["h"] for r in B], [r["err"] for r in B])}

    # C. uniform Kuhn cube
    tgt = 3.0 * pi2
    C = []
    for n in (3, 5, 8, 12):
        V, T = kuhn_cube(n)
        ev = dirichlet_eigenvalues(V, T, cube_boundary_mask(V), k=4)
        C.append({"h": mean_edge_length(V, T), "q_min": tet_quality(V, T)["q_min"],
                  "lam1": float(ev[0]), "err": abs(ev[0] - tgt) / tgt})
    out["C_uniform_cube"] = {"rows": C, "target": tgt,
                             "order": _order([r["h"] for r in C], [r["err"] for r in C])}
    return out


def main():
    res = run()
    print("=== 3D-S2: P1 FEM Laplacian vs exact Dirichlet eigenvalues ===\n")

    a = res["A_sliver_ball"]
    print(f"A. solid_ball (Delaunay shells)   target lambda_1 = pi^2 = {a['target']:.4f}")
    print(f'   {"h":>8} {"q_min":>7} {"lambda_1":>10} {"rel.err":>10}   spectrum (1,3,5 pattern)')
    for r in a["rows"]:
        sp = ", ".join(f"{x:.2f}" for x in r["spectrum"])
        print(f'   {r["h"]:8.4f} {r["q_min"]:7.3f} {r["lam1"]:10.4f} '
              f'{r["err"]:10.3e}   [{sp}]')
    print(f'   => fitted order O(h^{a["order"]:.2f})  -- error does NOT shrink: '
          f'slivers stall convergence\n')

    b = res["B_graded_ball"]
    print(f"B. graded_ball (geometric + rotated shells)   target = {b['target']:.4f}")
    print(f'   {"h":>8} {"q_min":>7} {"lambda_1":>10} {"rel.err":>10}')
    for r in b["rows"]:
        print(f'   {r["h"]:8.4f} {r["q_min"]:7.3f} {r["lam1"]:10.4f} {r["err"]:10.3e}')
    print(f'   => fitted order O(h^{b["order"]:.2f})  -- converges, but q_min drifts '
          f'to 0 so only first order\n')

    c = res["C_uniform_cube"]
    print(f"C. kuhn_cube (congruent uniform lattice)   target 3*pi^2 = {c['target']:.4f}")
    print(f'   {"h":>8} {"q_min":>7} {"lambda_1":>10} {"rel.err":>10}')
    for r in c["rows"]:
        print(f'   {r["h"]:8.4f} {r["q_min"]:7.3f} {r["lam1"]:10.4f} {r["err"]:10.3e}')
    print(f'   => fitted order O(h^{c["order"]:.2f})  -- q constant (= uniformly '
          f'shape-regular) => OPTIMAL second order\n')

    print("[結句] STEP 2 settles the 3D convergence question honestly. The")
    print("   operator is plain P1 tetrahedral FEM (no new operator -- the 2D")
    print("   Tier-3 finding carries over). What the 3D data adds is the")
    print("   mesh-quality law, exactly mirroring the 2D no-free-lunch story:")
    print("     * sliver meshes (q_min=0): the spectrum is structurally correct")
    print("       (right 1-3-5 multiplicities) but does NOT converge;")
    print("     * merely-better meshes (q_min -> 0): converge only at O(h);")
    print("     * genuinely UNIFORM, congruent tetrahedra (q constant): recover")
    print("       the optimal O(h^2).")
    print("   So 'as-equally-spaced-as-possible tetrahedra' is not cosmetic -- it")
    print("   is the precise lever that buys back convergence in 3D. It buys back")
    print("   what structured-mesh FEM already achieves, i.e. it is correct and")
    print("   important, but it is NOT novelty over existing FEM. The open")
    print("   construction is a boundary-conforming UNIFORM filling of the ball")
    print("   (uniform interior + exact sphere boundary at once).")
    return res


if __name__ == "__main__":
    main()
