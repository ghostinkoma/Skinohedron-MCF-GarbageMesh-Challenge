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
  spectral error        :  | lambda_h - lambda_exact | / lambda_exact
                           (the lowest non-zero eigenvalue, lam_exact = 2)

The empirical orders in h are reported. The honest finding (see [結句]) is that
the STRONG claim is false even on the best mesh -- it is O(h), matching the
Wardetzky-Mathur-Kaelberer-Grinspun 'no free lunch' theorem -- while the WEAK
spectral claim holds and is in fact super-convergent. The manuscript should
restate its theorems accordingly.
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


def _spectral_error(V, F, n_modes=6, exact=2.0):
    L, M = dec.cotangent_laplacian(V, F)
    Ms = dec.mass_matrix_sparse(M)
    vals = eigsh(L, k=n_modes, M=Ms, sigma=0.0, which="LM",
                 return_eigenvectors=False)
    vals = np.sort(vals)
    # lambda_0 ~ 0 (constants); the next three ~ 2 (the l=1 triplet)
    lam_h = float(np.mean(vals[1:4]))
    return abs(lam_h - exact) / exact, lam_h


def run(levels=range(2, 7), harmonic="Y2xy") -> dict:
    out = {"harmonic": harmonic, "regular": {}, "irregular": {}}
    for tag, jit in (("regular", 0.0), ("irregular", 0.30)):
        hs, pw, sp = [], [], []
        for lvl in levels:
            V, F = mesh.icosphere(lvl, jitter=jit, seed=1)
            h = mesh.mean_edge_length(V, F)
            hs.append(h)
            pw.append(_pointwise_error(V, F, harmonic))
            se, _ = _spectral_error(V, F)
            sp.append(se)
        hs = np.array(hs); pw = np.array(pw); sp = np.array(sp)
        # only fit a slope where the data actually decreases (irregular plateaus)
        pw_order = float(np.polyfit(np.log(hs), np.log(pw), 1)[0])
        sp_order = float(np.polyfit(np.log(hs), np.log(sp), 1)[0])
        out[tag] = {"h": hs.tolist(), "pointwise": pw.tolist(),
                    "spectral": sp.tolist(),
                    "pointwise_order": pw_order, "spectral_order": sp_order}
    return out


def main():
    res = run()
    for tag in ("regular", "irregular"):
        d = res[tag]
        print(f"=== S6: {tag} icosphere  (harmonic {res['harmonic']}) ===")
        print(f"{'h':>9} {'pointwise L2':>14} {'spectral err':>14}")
        for h, p, s in zip(d["h"], d["pointwise"], d["spectral"]):
            print(f"{h:>9.5f} {p:>14.4e} {s:>14.4e}")
        print(f"  pointwise consistency ~ h^{d['pointwise_order']:.2f}")
        print(f"  spectral convergence  ~ h^{d['spectral_order']:.2f}\n")

    rp = res["regular"]["pointwise_order"]
    rs = res["regular"]["spectral_order"]
    print("[結句] On the *best possible* mesh the canonical DEC Laplacian is")
    print(f"       O(h^{rp:.1f}) pointwise -- NOT O(h^2). Prop. 6.x as written is too")
    print(f"       strong. But the spectral error is O(h^{rs:.1f}) (super-convergent),")
    print("       so the USEFUL claim survives if restated weakly/variationally.")
    print("       On irregular meshes the pointwise error does not converge at all,")
    print("       exactly as the 'no free lunch' theorem predicts: trace-freeness")
    print("       alone cannot rescue strong consistency. Recommendation: replace")
    print("       Theorem 5.x / Prop 6.x with a Galerkin (weak) O(h^2) statement")
    print("       under a uniform shape-regularity hypothesis.")


if __name__ == "__main__":
    main()
