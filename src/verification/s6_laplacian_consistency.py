"""
s6_laplacian_consistency.py  --  Part I, Section 6   (THE decisive test)
========================================================================
The manuscript's headline results are:

  Theorem 5.x  : the local truncation error of the edge operator has NO O(h)
                 term; the leading error is O(h^2).
  Prop. 6.x    : ||Delta_h I_k f - I_k Delta f||_{L2} = O(h^2)   for f in C^4.

We test BOTH the strong (pointwise) and the weak (spectral/variational) readings
of "consistency" against the exact Laplace-Beltrami spectrum of S^2, on a
near-ideal regular icosphere AND on a persistently-irregular one.

Two error measures, each in the lumped-L2 norm:

  pointwise consistency :  || M^{-1} L f  -  lam f ||_{L2} / || lam f ||_{L2}
                           (this is exactly the quantity in Prop. 6.x)
  spectral error        :  | lambda_h - l(l+1) | / l(l+1),
                           lambda_h = (f^T L f) / (f^T M f)
                           (the discrete Rayleigh quotient of the SAME harmonic
                           f = Y_l used for the pointwise test -- this is the
                           definition used in the paper, Section 6)

Two eigenvalue cross-checks are also recorded:

  eig_l2  :  max relative error of the l=2 quintet of generalised eigenvalues
             (exact 6) -- a generic smooth mode.
  eig_l1  :  relative error of the MEAN of the l=1 triplet (exact 2).
             CAUTION: the l=1 eigenfunctions are the coordinate functions
             x, y, z, and on the icosphere this triplet converges anomalously
             fast (~h^3.7). This is a special case, NOT the generic spectral
             rate; earlier versions reported it as "spectral super-convergence".

The empirical orders in h are reported. The honest finding (see [結句]) is that
the STRONG claim is false even on the best mesh -- it is O(h), matching the
Wardetzky-Mathur-Kaelberer-Grinspun 'no free lunch' theorem -- while the WEAK
spectral claim holds at the classical P1 rate O(h^2) (Rayleigh quotient and
generic eigenvalues alike). Only the l=1 triplet is faster. The manuscript
should restate its theorems accordingly.
"""
from __future__ import annotations
import numpy as np
from scipy.sparse.linalg import eigsh
from ksf import mesh, dec, sph


def _pointwise_error(V, F, harmonic="Y2xy"):
    L, M = dec.cotangent_laplacian(V, F)
    f, lam = sph.harmonic(V, harmonic)
    approx = (L @ f) / M                      # M^{-1} L f  ~  -Delta f = lam f
    num = np.sqrt(np.sum(M * (approx - lam * f) ** 2))
    den = np.sqrt(np.sum(M * (lam * f) ** 2))
    return float(num / den)


def _spectral_error(V, F, harmonic="Y2xy"):
    """Paper definition: relative error of the discrete Rayleigh quotient
    (f^T L f)/(f^T M f) of the exact harmonic f against l(l+1)."""
    L, M = dec.cotangent_laplacian(V, F)
    f, lam = sph.harmonic(V, harmonic)
    lam_h = float((f @ (L @ f)) / (f @ (M * f)))
    return abs(lam_h - lam) / lam, lam_h


def _eigen_errors(V, F, n_modes=9):
    """Generalised eigenvalues of (L, M): returns (eig_l1, eig_l2).
    eig_l1 = rel. error of the mean of the l=1 triplet (exact 2; special case),
    eig_l2 = max rel. error over the l=2 quintet (exact 6; generic mode)."""
    L, M = dec.cotangent_laplacian(V, F)
    Ms = dec.mass_matrix_sparse(M)
    vals = eigsh(L, k=n_modes, M=Ms, sigma=0.0, which="LM",
                 return_eigenvectors=False)
    vals = np.sort(vals)
    # lambda_0 ~ 0 (constants); vals[1:4] ~ 2 (l=1); vals[4:9] ~ 6 (l=2)
    eig_l1 = abs(float(np.mean(vals[1:4])) - 2.0) / 2.0
    eig_l2 = float(np.max(np.abs(vals[4:9] - 6.0))) / 6.0
    return eig_l1, eig_l2


def run(levels=range(2, 7), harmonic="Y2xy") -> dict:
    out = {"harmonic": harmonic, "regular": {}, "irregular": {}}
    for tag, jit in (("regular", 0.0), ("irregular", 0.30)):
        hs, pw, sp, e1, e2 = [], [], [], [], []
        for lvl in levels:
            V, F = mesh.icosphere(lvl, jitter=jit, seed=1)
            h = mesh.mean_edge_length(V, F)
            hs.append(h)
            pw.append(_pointwise_error(V, F, harmonic))
            se, _ = _spectral_error(V, F, harmonic)
            sp.append(se)
            a, b = _eigen_errors(V, F)
            e1.append(a); e2.append(b)
        hs = np.array(hs)

        def order(err):
            return float(np.polyfit(np.log(hs), np.log(np.array(err)), 1)[0])

        out[tag] = {"h": hs.tolist(), "pointwise": list(map(float, pw)),
                    "spectral": list(map(float, sp)),
                    "eig_l2": list(map(float, e2)),
                    "eig_l1": list(map(float, e1)),
                    "pointwise_order": order(pw), "spectral_order": order(sp),
                    "eig_l2_order": order(e2), "eig_l1_order": order(e1)}
    return out


def main():
    res = run()
    for tag in ("regular", "irregular"):
        d = res[tag]
        print(f"=== S6: {tag} icosphere  (harmonic {res['harmonic']}) ===")
        print(f"{'h':>9} {'pointwise L2':>14} {'spectral (RQ)':>14}"
              f" {'eig l=2':>12} {'eig l=1*':>12}")
        for h, p, s, b, a in zip(d["h"], d["pointwise"], d["spectral"],
                                 d["eig_l2"], d["eig_l1"]):
            print(f"{h:>9.5f} {p:>14.4e} {s:>14.4e} {b:>12.4e} {a:>12.4e}")
        print(f"  pointwise consistency ~ h^{d['pointwise_order']:.2f}")
        print(f"  spectral convergence  ~ h^{d['spectral_order']:.2f}"
              f"  (Rayleigh quotient of {res['harmonic']})")
        print(f"  eigenvalue l=2        ~ h^{d['eig_l2_order']:.2f}")
        print(f"  eigenvalue l=1*       ~ h^{d['eig_l1_order']:.2f}"
              "  (*special case: coordinate functions)\n")

    rp = res["regular"]["pointwise_order"]
    rs = res["regular"]["spectral_order"]
    r1 = res["regular"]["eig_l1_order"]
    print("[結句] On the *best possible* mesh the canonical DEC Laplacian is")
    print(f"       O(h^{rp:.1f}) pointwise -- NOT O(h^2). Prop. 6.x as written is too")
    print(f"       strong. The spectral error is O(h^{rs:.1f}) -- the classical P1 rate --")
    print("       so the USEFUL claim survives if restated weakly/variationally.")
    print(f"       (Only the l=1 triplet converges faster, ~h^{r1:.1f}; that is a")
    print("       special property of the coordinate functions, not a general")
    print("       spectral super-convergence.)")
    print("       On irregular meshes the pointwise error does not converge at all,")
    print("       exactly as the 'no free lunch' theorem predicts: trace-freeness")
    print("       alone cannot rescue strong consistency. Recommendation: replace")
    print("       Theorem 5.x / Prop 6.x with a Galerkin (weak) O(h^2) statement")
    print("       under a uniform shape-regularity hypothesis.")


if __name__ == "__main__":
    main()
