"""
ksf3d.dec3d
===========
Discrete Exterior Calculus on a tetrahedral 3-complex (3D extension, Part I §3).

The de Rham complex in 3D has three exterior derivatives

    Omega^0 --d0--> Omega^1 --d1--> Omega^2 --d2--> Omega^3
   (vertices)     (edges)        (faces)        (tets)

We build them as the transposes of the simplicial boundary operators
d_p = (partial_{p+1})^T, using the *sorted-simplex with alternating signs*
convention:

    partial[v0,...,vp] = sum_i (-1)^i [v0,..,v^i,..,vp]    (vertices sorted)

With this convention the identity  partial . partial = 0  holds *exactly* in
integer arithmetic, hence  d1 . d0 = 0  and  d2 . d1 = 0  are machine-zero by
construction.  Verifying it numerically is still worthwhile: it confirms the
face/edge enumeration and indexing are wired correctly.
"""
from __future__ import annotations
import numpy as np
from scipy.sparse import coo_matrix
from itertools import combinations


def _enumerate_simplices(T: np.ndarray):
    """From tets T (M,4) build sorted, de-duplicated edge/face/tet index maps."""
    edges, faces, tets = {}, {}, {}
    for tet in T:
        t = tuple(sorted(int(x) for x in tet))
        tets.setdefault(t, len(tets))
        for f in combinations(t, 3):
            faces.setdefault(f, len(faces))
        for e in combinations(t, 2):
            edges.setdefault(e, len(edges))
    return edges, faces, tets


def boundary_operators(V: np.ndarray, T: np.ndarray):
    """Return dict with d0,d1,d2 (cochain coboundary) and the index maps.

    d0 : (E x V)   d1 : (F x E)   d2 : (Tt x F)
    """
    edges, faces, tets = _enumerate_simplices(T)
    nV, nE, nF, nT = len(V), len(edges), len(faces), len(tets)

    # partial_1 : edges -> vertices    [a,b] -> b - a
    r, c, val = [], [], []
    for (a, b), ei in edges.items():
        r += [a, b]; c += [ei, ei]; val += [-1.0, +1.0]
    d0 = coo_matrix((val, (c, r)), shape=(nE, nV)).tocsr()   # = partial_1^T

    # partial_2 : faces -> edges    [a,b,c] -> [b,c] - [a,c] + [a,b]
    r, c, val = [], [], []
    for (a, b, cc), fi in faces.items():
        for sgn, e in ((+1.0, (b, cc)), (-1.0, (a, cc)), (+1.0, (a, b))):
            r.append(fi); c.append(edges[e]); val.append(sgn)
    # coboundary d1 : Omega^1 -> Omega^2 is the face<-edge incidence (F x E),
    # which is exactly partial_2 with this sign convention.
    d1 = coo_matrix((val, (r, c)), shape=(nF, nE)).tocsr()

    # partial_3 : tets -> faces   [a,b,c,d] -> [b,c,d]-[a,c,d]+[a,b,d]-[a,b,c]
    r, c, val = [], [], []
    for (a, b, cc, dd), ti in tets.items():
        for sgn, f in ((+1.0, (b, cc, dd)), (-1.0, (a, cc, dd)),
                       (+1.0, (a, b, dd)), (-1.0, (a, b, cc))):
            r.append(ti); c.append(faces[f]); val.append(sgn)
    d2 = coo_matrix((val, (r, c)), shape=(nT, nF)).tocsr()

    return {"d0": d0, "d1": d1, "d2": d2,
            "edges": edges, "faces": faces, "tets": tets,
            "counts": (nV, nE, nF, nT)}


def complex_residuals(V: np.ndarray, T: np.ndarray):
    """Return (max|d1.d0|, max|d2.d1|, euler) ; both residuals should be 0."""
    ops = boundary_operators(V, T)
    d0, d1, d2 = ops["d0"], ops["d1"], ops["d2"]
    r1 = float(np.abs((d1 @ d0).toarray()).max())
    r2 = float(np.abs((d2 @ d1).toarray()).max())
    nV, nE, nF, nT = ops["counts"]
    euler = nV - nE + nF - nT
    return r1, r2, euler, ops["counts"]
