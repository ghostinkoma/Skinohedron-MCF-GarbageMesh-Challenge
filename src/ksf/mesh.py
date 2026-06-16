"""
ksf.mesh
========
Construction of the Kosaka Skin-o-hedron sequence S = {S_k} (paper Part I, Def. 2.1)
as a concrete, reproducible icosphere refinement family on S^2.

Each refinement level halves the mesh size h_k, giving the sequence h_k -> 0 the
paper requires. We also expose a *tangential jitter* option that produces a
persistently-irregular (non shape-improving) mesh family, used in s7 to probe the
"non-Delaunay stability" claim.

All meshes live exactly on the unit sphere so the ambient metric g is known in
closed form, which lets every other module compare discrete operators against an
exact ground truth (the Laplace-Beltrami spectrum).
"""
from __future__ import annotations
import numpy as np


# --------------------------------------------------------------------------- #
#  base icosahedron                                                           #
# --------------------------------------------------------------------------- #
def icosahedron() -> tuple[np.ndarray, np.ndarray]:
    """Return (V, F) of a unit icosahedron: 12 vertices, 20 triangular faces."""
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    V = np.array(
        [[-1, phi, 0], [1, phi, 0], [-1, -phi, 0], [1, -phi, 0],
         [0, -1, phi], [0, 1, phi], [0, -1, -phi], [0, 1, -phi],
         [phi, 0, -1], [phi, 0, 1], [-phi, 0, -1], [-phi, 0, 1]], float)
    V /= np.linalg.norm(V, axis=1, keepdims=True)
    F = np.array(
        [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
         [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
         [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
         [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]], int)
    return V, F


def _subdivide(V: np.ndarray, F: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One step of 1->4 midpoint subdivision with re-projection onto the sphere.

    Implements the refinement map phi_k of Def. 2.1: every face subdivides
    locally, incidence is preserved, and the new vertices are pushed to S^2.
    """
    verts = list(map(tuple, V))
    index = {v: i for i, v in enumerate(verts)}
    cache: dict[tuple[int, int], int] = {}

    def midpoint(a: int, b: int) -> int:
        key = (min(a, b), max(a, b))
        if key in cache:
            return cache[key]
        p = (np.asarray(verts[a]) + np.asarray(verts[b])) * 0.5
        p /= np.linalg.norm(p)
        t = tuple(p)
        if t not in index:
            index[t] = len(verts)
            verts.append(t)
        cache[key] = index[t]
        return index[t]

    new_faces = []
    for a, b, c in F:
        ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
        new_faces += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
    return np.asarray(verts, float), np.asarray(new_faces, int)


# --------------------------------------------------------------------------- #
#  public constructor                                                         #
# --------------------------------------------------------------------------- #
def icosphere(level: int, jitter: float = 0.0, seed: int = 0
              ) -> tuple[np.ndarray, np.ndarray]:
    """Build S_level.

    Parameters
    ----------
    level  : number of 1->4 subdivisions of the base icosahedron.
    jitter : if > 0, add a tangential perturbation of size  jitter * h  to every
             vertex then re-project to the sphere. Because the perturbation
             scales with h, the family stays *persistently* irregular under
             refinement (the relative distortion does not vanish) -- exactly the
             stress case for the non-Delaunay claim.
    seed   : RNG seed for the jitter.
    """
    V, F = icosahedron()
    for _ in range(level):
        V, F = _subdivide(V, F)

    if jitter > 0.0:
        rng = np.random.default_rng(seed)
        n = V / np.linalg.norm(V, axis=1, keepdims=True)
        r = rng.normal(size=V.shape)
        tangential = r - np.sum(r * n, axis=1, keepdims=True) * n
        h = mean_edge_length(V, F)
        V = V + jitter * h * tangential
        V /= np.linalg.norm(V, axis=1, keepdims=True)
    return V, F


# --------------------------------------------------------------------------- #
#  mesh quality / metric-fidelity diagnostics (Def. 2.1 axioms 1 & 3)         #
# --------------------------------------------------------------------------- #
def edges(F: np.ndarray) -> np.ndarray:
    """Unique undirected edges as an (E,2) int array."""
    s = set()
    for a, b, c in F:
        for u, v in ((a, b), (b, c), (c, a)):
            s.add((min(u, v), max(u, v)))
    return np.asarray(sorted(s), int)


def mean_edge_length(V: np.ndarray, F: np.ndarray) -> float:
    e = edges(F)
    return float(np.linalg.norm(V[e[:, 0]] - V[e[:, 1]], axis=1).mean())


def shape_regularity(V: np.ndarray, F: np.ndarray) -> dict:
    """Return min/mean/max of the per-triangle quality  q = 4*sqrt(3)*A / sum(l^2).

    q = 1 for an equilateral triangle, q -> 0 for a degenerate sliver. A family
    is shape-regular (Def. 2.1) iff inf_k min(q) > 0.
    """
    qs = []
    for a, b, c in F:
        va, vb, vc = V[a], V[b], V[c]
        area = 0.5 * np.linalg.norm(np.cross(vb - va, vc - va))
        l2 = (np.sum((vb - va) ** 2) + np.sum((vc - vb) ** 2)
              + np.sum((va - vc) ** 2))
        qs.append(4.0 * np.sqrt(3.0) * area / l2)
    qs = np.asarray(qs)
    return {"min": float(qs.min()), "mean": float(qs.mean()),
            "max": float(qs.max())}


def vertex_normals(V: np.ndarray, F: np.ndarray) -> np.ndarray:
    """Angle-weighted vertex normals (handy for the viewer / shading)."""
    n = np.zeros_like(V)
    for a, b, c in F:
        fn = np.cross(V[b] - V[a], V[c] - V[a])
        nn = np.linalg.norm(fn)
        if nn > 0:
            fn = fn / nn
        n[a] += fn; n[b] += fn; n[c] += fn
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    ln[ln == 0] = 1.0
    return n / ln
