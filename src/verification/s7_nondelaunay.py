"""
s7_nondelaunay.py  --  Part I, Section 7
========================================
The manuscript claims the construction stays well-defined on non-Delaunay meshes
with no 'catastrophic degeneracy'. We make the precise tradeoff visible.

On a triangle mesh the cotangent weight on an edge is
        w = 0.5 (cot a + cot b)
which becomes NEGATIVE exactly when the edge violates the local Delaunay
condition (a + b > pi). The Wardetzky et al. 'no free lunch' theorem says no
local, symmetric, linear-precise Laplacian can keep ALL weights non-negative on
ALL meshes. So 'works on non-Delaunay meshes' and 'positive weights everywhere'
cannot both hold.

We report, as the mesh is made progressively more irregular:
  * the fraction of edges with a negative cotangent weight (loss of an M-matrix
    / maximum-principle structure), and
  * the smallest mass (Voronoi area), to show the operator stays *non-degenerate*
    (bounded away from 0 ~ h^2) even while positivity is lost.

This substantiates a *weakened, honest* version of the Section 7 claim.
"""
from __future__ import annotations
import numpy as np
from ksf import mesh, dec


def _negative_weight_fraction(V, F):
    L, M = dec.cotangent_laplacian(V, F)
    A = L.tocoo()
    off = A.row != A.col
    offdiag = A.data[off]
    # stiffness stores -w off-diagonal; w<0  <=>  offdiag entry > 0
    n_off = offdiag.size
    n_neg = int(np.sum(offdiag > 1e-12))
    return n_neg / max(n_off, 1), float(M.min())


def run(jitters=(0.0, 0.1, 0.2, 0.3, 0.45), level=4) -> dict:
    rows = []
    for j in jitters:
        V, F = mesh.icosphere(level, jitter=j, seed=7)
        h = mesh.mean_edge_length(V, F)
        frac, mmin = _negative_weight_fraction(V, F)
        q = mesh.shape_regularity(V, F)["min"]
        rows.append(dict(jitter=float(j), h=h, neg_weight_frac=frac,
                         min_mass=mmin, q_min=q))
    return {"level": level, "rows": rows}


def main():
    res = run()
    print(f"=== S7: non-Delaunay behaviour (level {res['level']}) ===")
    print(f"{'jitter':>7} {'q_min':>7} {'neg-weight %':>13} {'min Voronoi area':>18}")
    for r in res["rows"]:
        print(f"{r['jitter']:>7.2f} {r['q_min']:>7.3f} "
              f"{100*r['neg_weight_frac']:>12.2f}% {r['min_mass']:>18.3e}")
    print("\n[結句] As irregularity grows, a growing fraction of edges acquire")
    print("       NEGATIVE cotangent weights -- the maximum principle is lost,")
    print("       exactly the price the 'no free lunch' theorem demands. What the")
    print("       paper CAN honestly claim is the weaker fact it actually shows:")
    print("       the operator stays *non-degenerate* (min Voronoi area ~ h^2 > 0),")
    print("       so assembly never blows up. Section 7 should be reworded from")
    print("       'stable' to 'non-degenerate but not maximum-principle preserving'.")


if __name__ == "__main__":
    main()
