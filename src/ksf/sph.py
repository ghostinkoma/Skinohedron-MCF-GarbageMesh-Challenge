"""
ksf.sph
=======
Exact ground truth for every convergence test: real spherical harmonics on S^2,
which are the eigenfunctions of the Laplace-Beltrami operator with

        Delta_{S^2} Y_l^m = - l(l+1) Y_l^m .

We expose a small set of low-degree real harmonics as closed-form polynomials in
the Cartesian coordinates of points already lying on the unit sphere, together
with their eigenvalues lam = l(l+1) (the positive eigenvalue of -Delta).
"""
from __future__ import annotations
import numpy as np

# name -> (callable(V)->values, eigenvalue l(l+1))
_HARMONICS = {
    # l = 1  (eigenvalue 2)
    "Y10": (lambda V: V[:, 2], 2.0),                         # ~ z
    "Y11": (lambda V: V[:, 0], 2.0),                         # ~ x
    # l = 2  (eigenvalue 6)
    "Y2xy": (lambda V: V[:, 0] * V[:, 1], 6.0),              # ~ xy
    "Y2z2": (lambda V: 3.0 * V[:, 2] ** 2 - 1.0, 6.0),       # ~ 3z^2 - 1
    # l = 3  (eigenvalue 12)
    "Y3xyz": (lambda V: V[:, 0] * V[:, 1] * V[:, 2], 12.0),  # ~ xyz
}


def harmonic(V: np.ndarray, name: str = "Y2xy") -> tuple[np.ndarray, float]:
    """Return (f_at_vertices, eigenvalue l(l+1)) for a named real harmonic."""
    fn, lam = _HARMONICS[name]
    f = fn(V).astype(float)
    return f, lam


def available() -> list[str]:
    return list(_HARMONICS)


def l2_norm(values: np.ndarray, mass: np.ndarray) -> float:
    """Discrete L^2 norm with the lumped mass weights."""
    return float(np.sqrt(np.sum(mass * values ** 2)))
