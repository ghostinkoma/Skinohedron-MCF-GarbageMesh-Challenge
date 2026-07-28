"""
ksf3d.fem3d
===========
The 3D operator the substrate was missing. `dec3d` fixed the *complex*
(d0,d1,d2 with d.d = 0); this module places the lowest-order **P1 finite-element
Laplacian** on a tetrahedral mesh and provides the *exact* ground truth needed to
test it.

Honesty up front (the 2D Tier-3 result carries over): this is plain linear
tetrahedral FEM. It is NOT a new operator. Its only job here is to let us measure
*convergence* against known eigenvalues, so we can see -- with numbers -- how the
3D no-free-lunch / mesh-quality tension behaves.

Ground truth
------------
Dirichlet `-Delta u = lambda u`, `u|boundary = 0`:

  * unit ball  B^3 :  lambda_{l,k} = j_{l,k}^2 , the squared zeros of the
    spherical Bessel function j_l ; multiplicity 2l+1.
  * unit cube  [0,1]^3 :  lambda = pi^2 (a^2 + b^2 + c^2) , integers a,b,c >= 1.
"""
from __future__ import annotations
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.special import spherical_jn
from scipy.optimize import brentq


# --------------------------------------------------------------------------- #
#  P1 FEM assembly on tetrahedra                                              #
# --------------------------------------------------------------------------- #
def fem_laplacian(V: np.ndarray, T: np.ndarray):
    """Return (K, M): K = sparse P1 stiffness (SPSD ~ -Delta),
    M = lumped (barycentric) mass as a dense 1-D array of vertex volumes."""
    n = len(V)
    rows, cols, vals = [], [], []
    M = np.zeros(n)
    n_degenerate = 0
    for tet in T:
        ix = [int(i) for i in tet]
        P = V[ix]
        vol = abs(np.linalg.det(
            np.array([P[1] - P[0], P[2] - P[0], P[3] - P[0]]))) / 6.0
        if vol <= 1e-15:
            n_degenerate += 1
            continue
        # P1 gradients: rows 1..3 of inv([1 x y z]) are the constant gradients
        C = np.linalg.inv(np.column_stack([np.ones(4), P]))
        G = C[1:4, :]                       # (3,4) gradient of each hat
        Ke = vol * (G.T @ G)
        for a in range(4):
            M[ix[a]] += vol / 4.0
            for b in range(4):
                rows.append(ix[a]); cols.append(ix[b]); vals.append(Ke[a, b])
    if n_degenerate:
        import warnings
        warnings.warn(
            f"fem_laplacian: {n_degenerate} degenerate tet(s) skipped (vol <= 1e-15); "
            "assembled operator may be incomplete.",
            stacklevel=2,
        )
    K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    return K, M


def dirichlet_eigenvalues(V, T, boundary_mask, k=6):
    """Lowest k Dirichlet eigenvalues of K u = lambda M u with u=0 on the
    boundary vertices flagged by `boundary_mask` (bool array over vertices)."""
    K, M = fem_laplacian(V, T)
    I = np.where(~boundary_mask)[0]
    Ki = K[I][:, I].tocsr()
    Mi = coo_matrix((M[I], (range(len(I)), range(len(I))))).tocsr()
    vals = eigsh(Ki, k=k, M=Mi, sigma=0.0, which="LM", return_eigenvectors=False)
    return np.sort(vals)


# --------------------------------------------------------------------------- #
#  exact ground-truth eigenvalues                                             #
# --------------------------------------------------------------------------- #
def _spherical_bessel_zeros(l, n, xmax=80.0, samples=8000):
    xs = np.linspace(1e-4, xmax, samples)
    f = spherical_jn(l, xs)
    zeros = []
    for i in range(len(xs) - 1):
        if f[i] == 0.0 or f[i] * f[i + 1] < 0.0:
            zeros.append(brentq(lambda x: spherical_jn(l, x), xs[i], xs[i + 1]))
            if len(zeros) >= n:
                break
    return zeros


def ball_dirichlet_spectrum(n_levels=8):
    """Lowest distinct Dirichlet eigenvalues of the unit ball, with multiplicity:
    returns sorted list of (lambda, l, multiplicity)."""
    out = []
    for l in range(0, 6):
        for z in _spherical_bessel_zeros(l, 4):
            out.append((z * z, l, 2 * l + 1))
    out.sort()
    return out[:n_levels]


def cube_dirichlet_spectrum(n_levels=8):
    """Lowest Dirichlet eigenvalues of the unit cube [0,1]^3: pi^2 (a^2+b^2+c^2)."""
    vals = {}
    for a in range(1, 6):
        for b in range(1, 6):
            for c in range(1, 6):
                lam = np.pi ** 2 * (a * a + b * b + c * c)
                vals[round(lam, 6)] = vals.get(round(lam, 6), 0) + 1
    items = sorted(vals.items())
    return [(lam, m) for lam, m in items][:n_levels]
