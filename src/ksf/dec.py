"""
ksf.dec
=======
Discrete Exterior Calculus operators on a triangulated surface (paper Part I, §3).

We provide:
  * the coboundary (discrete exterior derivative) d0 : Omega^0 -> Omega^1 and
    d1 : Omega^1 -> Omega^2, and verify the complex property d1 . d0 = 0;
  * the cotangent stiffness matrix L (the canonical primal/dual Hodge Laplacian,
    L = d0^T * star1 * d0 with the cotan Hodge-star), and the lumped Voronoi mass
    matrix M. The pair (L, M) is exactly the linear-FEM Laplace-Beltrami
    discretisation, which is the reference operator we hold the paper's claims to.

Sign convention: L is symmetric positive semidefinite and approximates the
*positive* operator  -Delta. Hence for an eigenfunction with  Delta f = -lam f
we have, weakly,  L f ~ lam * M f, i.e.  M^{-1} L f ~ lam f.
"""
from __future__ import annotations
import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
from .mesh import edges


# --------------------------------------------------------------------------- #
#  coboundary operators                                                       #
# --------------------------------------------------------------------------- #
def coboundary(V: np.ndarray, F: np.ndarray):
    """Return (d0, d1, edge_list, edge_index).

    d0 : (E x N) signed incidence vertices->edges.
    d1 : (T x E) signed incidence edges->faces (consistent orientation).
    """
    E = edges(F)
    eidx = {(a, b): i for i, (a, b) in enumerate(E)}
    n, ne, nt = len(V), len(E), len(F)

    # d0
    rows, cols, vals = [], [], []
    for i, (a, b) in enumerate(E):           # oriented a->b
        rows += [i, i]; cols += [b, a]; vals += [1.0, -1.0]
    d0 = coo_matrix((vals, (rows, cols)), shape=(ne, n)).tocsr()

    # d1
    rows, cols, vals = [], [], []
    for t, (a, b, c) in enumerate(F):
        for u, v in ((a, b), (b, c), (c, a)):
            key = (min(u, v), max(u, v))
            sign = 1.0 if (u, v) == key else -1.0
            rows.append(t); cols.append(eidx[key]); vals.append(sign)
    d1 = coo_matrix((vals, (rows, cols)), shape=(nt, ne)).tocsr()
    return d0, d1, E, eidx


def complex_residual(V: np.ndarray, F: np.ndarray) -> float:
    """max |d1 . d0|  -- should be ~0 (the de Rham complex closes)."""
    d0, d1, _, _ = coboundary(V, F)
    return float(np.abs((d1 @ d0).toarray()).max())


# --------------------------------------------------------------------------- #
#  cotangent Laplacian + mass                                                 #
# --------------------------------------------------------------------------- #
def _cot(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """cot of the interior angle at vertex a in triangle (a,b,c)."""
    u, v = b - a, c - a
    return float(np.dot(u, v) / np.linalg.norm(np.cross(u, v)))


def cotangent_laplacian(V: np.ndarray, F: np.ndarray):
    """Return (L, M) : L = sparse cotan stiffness (SPSD ~ -Delta),
    M = lumped Voronoi (barycentric) mass as a dense 1-D array of areas."""
    n = len(V)
    rows, cols, vals = [], [], []
    M = np.zeros(n)
    for i, j, k in F:
        vi, vj, vk = V[i], V[j], V[k]
        ci = _cot(vi, vj, vk)   # angle at i, opposite edge (j,k)
        cj = _cot(vj, vk, vi)   # angle at j, opposite edge (k,i)
        ck = _cot(vk, vi, vj)   # angle at k, opposite edge (i,j)
        for a, b, w in ((j, k, ci), (k, i, cj), (i, j, ck)):
            hw = 0.5 * w
            rows += [a, b, a, b]
            cols += [b, a, a, b]
            vals += [-hw, -hw, hw, hw]
        area = 0.5 * np.linalg.norm(np.cross(vj - vi, vk - vi))
        M[i] += area / 3.0; M[j] += area / 3.0; M[k] += area / 3.0
    L = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    return L, M


def mass_matrix_sparse(M: np.ndarray) -> csr_matrix:
    n = len(M)
    return coo_matrix((M, (range(n), range(n))), shape=(n, n)).tocsr()
