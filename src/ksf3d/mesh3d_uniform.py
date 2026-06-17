"""
ksf3d.mesh3d_uniform
====================
"As equally-spaced as possible" tetrahedral fillings -- the meshes that matter
for *convergence* (not just for the combinatorial substrate).

The original `mesh3d.solid_ball` Delaunay-tetrahedralises concentric icosphere
shells. It is fine for the topological checks (d.d=0, chi=1) but its tetrahedral
quality is q_min = 0 (Delaunay slivers between nearby shells), and -- as
`s3d_laplacian` shows -- that destroys eigenvalue convergence. This module
provides two better fillings:

  graded_ball : keep the exact spherical boundary, but (1) space the shells
                GEOMETRICALLY so radial and angular spacing match (isotropic
                cells), and (2) rotate each shell so vertices are not radially
                aligned. Quality improves a lot (q_mean ~ 0.8), but q_min still
                drifts toward 0 under refinement, so the family is NOT uniformly
                shape-regular and convergence is only ~O(h).

  kuhn_cube   : the unit cube split into congruent Kuhn (Freudenthal)
                tetrahedra -- 6 per sub-cube, every tet identical up to
                isometry. Quality is CONSTANT (q = 0.717) at every refinement,
                i.e. a uniformly shape-regular family. This is the genuinely
                "equally-spaced" tetrahedral lattice, and it restores clean
                O(h^2) eigenvalue convergence. Its boundary conforms exactly to
                the cube, so it isolates the mesh-quality effect with no
                domain-approximation confound.

Boundary-conforming UNIFORM filling of the *ball* (uniform interior + exact
sphere boundary at once) is the remaining construction; it needs a conforming
mesher (isosurface-stuffing / BCC-snap), noted as future work.
"""
from __future__ import annotations
import os, sys, itertools
import numpy as np
from scipy.spatial import Delaunay
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from ksf.mesh import icosphere, mean_edge_length as _surf_edge


def _rot(seed):
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def graded_ball(shell_level=1):
    """Solid unit ball: geometric radial shells matched to angular spacing, each
    shell rotated to avoid radial alignment. Exact spherical boundary."""
    Sv, _ = icosphere(shell_level)
    Sv = Sv / np.linalg.norm(Sv, axis=1, keepdims=True)
    hs = _surf_edge(*icosphere(shell_level))           # unit-sphere edge length
    radii = [1.0]
    while radii[-1] / (1.0 + hs) > 0.5 * hs:
        radii.append(radii[-1] / (1.0 + hs))
    radii = sorted(radii)
    pts = [np.zeros((1, 3))]
    for k, r in enumerate(radii):
        pts.append((Sv * r) @ _rot(k + 1).T)
    P = np.vstack(pts)
    T = Delaunay(P).simplices.astype(np.int64)
    return P, T


def ball_boundary_mask(V, tol=1e-6):
    return np.linalg.norm(V, axis=1) >= 1.0 - tol


def kuhn_cube(n=4):
    """Unit cube [0,1]^3 as n^3 sub-cubes, each split into 6 congruent Kuhn
    tetrahedra (all sharing the main diagonal). Every tet is identical up to
    isometry -> uniformly shape-regular at all n."""
    def idx(i, j, k):
        return (i * (n + 1) + j) * (n + 1) + k
    V = np.array([[i / n, j / n, k / n]
                  for i in range(n + 1) for j in range(n + 1)
                  for k in range(n + 1)], float)
    e = np.eye(3, dtype=int)
    T = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                base = np.array([i, j, k])
                for sg in itertools.permutations(range(3)):
                    cur = base.copy(); path = [cur.copy()]
                    for d in sg:
                        cur = cur + e[d]; path.append(cur.copy())
                    T.append([idx(*p) for p in path])
    return V, np.array(T, int)


def cube_boundary_mask(V):
    return (np.isclose(V, 0.0) | np.isclose(V, 1.0)).any(axis=1)


def kuhn_box(nx, ny, nz, Lx=2.0):
    """A rectangular box [−Lx/2, Lx/2] × [0, wy] × [0, wz] meshed by congruent
    Kuhn tetrahedra of edge h = Lx/nx (wy = ny·h, wz = nz·h). Used for the
    quasi-1D plane-wave channel of theory/01c. Same congruent cells as kuhn_cube;
    only the aspect ratio differs (long in x, thin in y,z)."""
    h = Lx / nx

    def idx(i, j, k):
        return (i * (ny + 1) + j) * (nz + 1) + k

    V = np.array([[i * h, j * h, k * h]
                  for i in range(nx + 1) for j in range(ny + 1)
                  for k in range(nz + 1)], float)
    e = np.eye(3, dtype=int)
    T = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                base = np.array([i, j, k])
                for sg in itertools.permutations(range(3)):
                    cur = base.copy(); path = [cur.copy()]
                    for d in sg:
                        cur = cur + e[d]; path.append(cur.copy())
                    T.append([idx(*p) for p in path])
    V[:, 0] -= Lx / 2.0     # centre x at 0 so the interface sits at x=0
    return V, np.array(T, int)
