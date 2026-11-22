"""
fluid_highorder_accuracy_verify.py  --  where the Taylor-Green error lives, and 1-D enrichment (theory/06f)

RECONSTRUCTED (2026-11). The original script was cited by theory/06f and the viewers
but never committed. This version reproduces every number of the 06f tables with the
repository's own operators (to the printed digits; P8 differs at the 1e-13 level).

  CHECK A  Taylor-Green mode u = sin(2 pi x) cos(2 pi y) on the PERIODIC Kuhn cube
           (the Box of fluid_navier_stokes_verify.py), n = 8, 10, 12, 16:
             stiffness  u^T K u  vs  2 k^2 E   (k = 2 pi, E = 1/4)   -> O(h^2) error
             mass       u^T M u  vs  E                           -> exact
           The mass exactness is a property of this test field (trapezoidal rule on a
           whole-period sinusoid, see mass_quadrature_verify.py), not of lumped mass
           in general.
  CHECK B  1-D enrichment: periodic [0, 1], Lagrange P_p elements on Gauss-Lobatto
           nodes, consistent mass, ne = round(12 / p) elements (~12 DOF); relative
           error of the eigenvalue (2 pi)^2 of -d^2/dx^2 for p = 1..6, 8.
           Spectral decay to machine precision -- in 1-D ONLY (no 2-D/3-D high-order
           assembly exists in this repository).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from numpy.polynomial import legendre as Leg
import numpy.polynomial.polynomial as Poly
from scipy.linalg import eigh
from fluid_navier_stokes_verify import Box

K_WAVE = 2 * np.pi
# values printed in theory/06f (for comparison)
DOC_TG = {8: 18.745, 10: 19.098, 12: 19.292, 16: 19.487}
DOC_1D = {1: 2.30e-2, 2: 1.58e-3, 3: 1.37e-4, 4: 1.29e-5, 5: 2.64e-7, 6: 3.43e-9, 8: 2.35e-13}


def check_A():
    print("\n[A] Taylor-Green mode on the periodic Kuhn cube: stiffness vs mass")
    target = 2 * K_WAVE**2 * 0.25
    errs = []
    for n in (8, 10, 12, 16):
        box = Box(n)
        ux = np.sin(K_WAVE * box.Vr[:, 0]) * np.cos(K_WAVE * box.Vr[:, 1])
        uKu = float(ux @ (box.Kp @ ux)); uMu = float(np.sum(box.Mp * ux * ux))
        e = uKu / target - 1.0
        errs.append(abs(e))
        print(f"    n={n:>2}: u^T K u = {uKu:.3f} (06f: {DOC_TG[n]:.3f}, err {e*100:+.2f}%)"
              f"   u^T M u = {uMu:.4f} (err {abs(uMu-0.25)/0.25:.1e})")
        assert abs(uKu - DOC_TG[n]) < 1e-3 and abs(uMu - 0.25) < 1e-12
    h = 1.0 / np.array([8, 10, 12, 16], float)
    order = float(np.polyfit(np.log(h), np.log(errs), 1)[0])
    print(f"    stiffness error order ~ h^{order:.2f}")
    assert 1.8 < order < 2.2
    print("    PASS  -- stiffness O(h^2); mass exact for THIS field (06f table reproduced)")


def _gll(p):
    if p == 1:
        return np.array([-1.0, 1.0])
    c = np.zeros(p + 1); c[-1] = 1.0
    inner = np.sort(np.real(Leg.legroots(Leg.legder(c))))
    return np.concatenate([[-1.0], inner, [1.0]])


def _element(p):
    x = _gll(p); q, w = Leg.leggauss(p + 4)
    coef = np.linalg.inv(np.vander(x, p + 1, increasing=True))
    phi = np.array([Poly.polyval(q, coef[:, i]) for i in range(p + 1)])
    dphi = np.array([Poly.polyval(q, Poly.polyder(coef[:, i])) for i in range(p + 1)])
    return (dphi * w) @ dphi.T, (phi * w) @ phi.T          # on the reference [-1, 1]


def eig_error_1d(p, ne):
    Ke, Me = _element(p); h = 1.0 / ne; N = ne * p
    K = np.zeros((N, N)); M = np.zeros((N, N))
    for e in range(ne):
        idx = [(e * p + i) % N for i in range(p + 1)]     # periodic
        for a in range(p + 1):
            for b in range(p + 1):
                K[idx[a], idx[b]] += Ke[a, b] * (2.0 / h)
                M[idx[a], idx[b]] += Me[a, b] * (h / 2.0)
    lam = np.sort(eigh(K, M, eigvals_only=True))
    return abs(lam[1] - K_WAVE**2) / K_WAVE**2, N


def check_B():
    print("\n[B] 1-D local enrichment P1 -> P8 (periodic, ~12 DOF)")
    for p in (1, 2, 3, 4, 5, 6, 8):
        ne = max(1, round(12 / p))
        err, N = eig_error_1d(p, ne)
        print(f"    P{p}: {ne:>2} elements, {N:>2} DOF, eigenvalue error {err:.2e}   (06f: {DOC_1D[p]:.2e})")
        if p <= 6:
            assert abs(err - DOC_1D[p]) / DOC_1D[p] < 0.01, "does not reproduce 06f"
        else:
            assert err < 1e-12
    print("    PASS  -- spectral decay to machine precision, in 1-D (06f table reproduced;")
    print("            P8 agrees at the 1e-13 round-off level)")


if __name__ == "__main__":
    print("06f: Taylor-Green error split (stiffness vs mass) and 1-D high-order enrichment")
    check_A()
    check_B()
    print("\nConclusion: for this benchmark the error lives in the stiffness (O(h^2)); the")
    print("mass form is exact for this whole-period sinusoid. Local enrichment reaches")
    print("machine precision in 1-D; 2-D/3-D enrichment is not implemented here.")
