"""
ksf3d.mesh3d
============
Solid-ball tetrahedral meshes for the 3D extension of the Skin-o-hedron model
(Part I, 3D).  This is the 3D analogue of the 2D `ksf.mesh` icosphere.

Construction
------------
We fill the closed unit ball B^3 = { |x| <= 1 } with tetrahedra by:

  1. taking `L` concentric spherical shells at radii  r_l = l / L   (l = 1..L),
     each populated by the vertices of a level-`s` icosphere scaled to r_l;
  2. adding the single centre vertex at the origin;
  3. computing the Delaunay tetrahedralisation of that point set.

Refinement is controlled by (s, L): increasing either lowers the mesh size h.
This is deliberately simple and dependency-light (scipy.spatial.Delaunay only),
so the result is easy to reproduce and audit.  It is NOT claimed to be an
optimal-quality mesh -- quality diagnostics are reported honestly by
`tet_quality()` so the reader can see exactly what they are getting.

Outputs are plain numpy arrays:
  V : (N,3) float64 vertex positions
  T : (M,4) int     tetrahedra (each row = 4 vertex indices)
"""
from __future__ import annotations
import numpy as np
from scipy.spatial import Delaunay

# reuse the 2D icosphere surface for the shell vertices
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from ksf.mesh import icosphere


# --------------------------------------------------------------------------- #
#  solid ball mesh                                                            #
# --------------------------------------------------------------------------- #
def solid_ball(shell_level: int = 1, n_shells: int = 3,
               dedup_tol: float = 1e-9):
    """Return (V, T) for a tetrahedralised solid unit ball.

    shell_level : icosphere subdivision level used on each shell (more => denser
                  angular resolution).
    n_shells    : number of radial shells (more => denser radial resolution).
    """
    Sv, _ = icosphere(shell_level)                 # unit-sphere surface verts
    Sv = Sv / np.linalg.norm(Sv, axis=1, keepdims=True)

    pts = [np.zeros((1, 3))]                        # centre
    for l in range(1, n_shells + 1):
        r = l / n_shells
        pts.append(Sv * r)
    P = np.vstack(pts)

    # deduplicate (the centre is unique; shells don't collide, but be safe)
    P = _dedup(P, dedup_tol)

    tri = Delaunay(P)                              # 3-D Delaunay => tetrahedra
    T = tri.simplices.astype(np.int64)

    # orient every tet positively (det >= 0) for a consistent convention
    T = _orient_positive(P, T)
    return P, T


def _dedup(P: np.ndarray, tol: float) -> np.ndarray:
    keep, seen = [], {}
    for p in P:
        key = tuple(np.round(p / tol).astype(np.int64))
        if key not in seen:
            seen[key] = True
            keep.append(p)
    return np.asarray(keep, dtype=np.float64)


def _orient_positive(V: np.ndarray, T: np.ndarray) -> np.ndarray:
    out = T.copy()
    for m, (a, b, c, d) in enumerate(T):
        mat = np.array([V[b] - V[a], V[c] - V[a], V[d] - V[a]])
        if np.linalg.det(mat) < 0:
            out[m, 2], out[m, 3] = T[m, 3], T[m, 2]   # swap to flip sign
    return out


# --------------------------------------------------------------------------- #
#  geometry diagnostics                                                       #
# --------------------------------------------------------------------------- #
def tet_volume(V: np.ndarray, tet) -> float:
    a, b, c, d = (V[i] for i in tet)
    return abs(np.linalg.det(np.array([b - a, c - a, d - a]))) / 6.0


def mean_edge_length(V: np.ndarray, T: np.ndarray) -> float:
    s = set()
    for a, b, c, d in T:
        for u, v in ((a, b), (a, c), (a, d), (b, c), (b, d), (c, d)):
            s.add((min(u, v), max(u, v)))
    e = np.array(list(s))
    return float(np.linalg.norm(V[e[:, 0]] - V[e[:, 1]], axis=1).mean())


def tet_quality(V: np.ndarray, T: np.ndarray) -> dict:
    """Radius-ratio quality rho = 3 r_in / r_circ in [0,1] (1 = regular tet).

    Returns min / mean over all tets, plus total volume and the ideal 4*pi/3.
    """
    q = np.empty(len(T))
    vol_tot = 0.0
    for m, tet in enumerate(T):
        a, b, c, d = (V[i] for i in tet)
        vol = tet_volume(V, tet)
        vol_tot += vol
        # face areas
        faces = [(a, b, c), (a, b, d), (a, c, d), (b, c, d)]
        A = sum(0.5 * np.linalg.norm(np.cross(y - x, z - x))
                for x, y, z in faces)
        # inradius r_in = 3V / A_total
        r_in = 3.0 * vol / A if A > 0 else 0.0
        # circumradius via the standard formula
        r_circ = _circumradius(a, b, c, d, vol)
        q[m] = (3.0 * r_in / r_circ) if r_circ > 0 else 0.0
    return {
        "q_min": float(q.min()),
        "q_mean": float(q.mean()),
        "volume": float(vol_tot),
        "volume_ideal": float(4.0 / 3.0 * np.pi),
    }


def _circumradius(a, b, c, d, vol) -> float:
    if vol <= 0:
        return 0.0
    # |product of opposite edge-pair sums| / (24 V)  -- standard tet circumradius
    A = np.linalg.norm(b - a); a_ = np.linalg.norm(d - c)
    B = np.linalg.norm(c - a); b_ = np.linalg.norm(d - b)
    C = np.linalg.norm(d - a); c_ = np.linalg.norm(c - b)
    p = A * a_; q = B * b_; r = C * c_
    s = (p + q + r) * (-p + q + r) * (p - q + r) * (p + q - r)
    if s <= 0:
        return 0.0
    return float(np.sqrt(s) / (24.0 * vol))
