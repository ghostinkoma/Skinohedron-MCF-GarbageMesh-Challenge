r"""
fix02_geometric_anisotropy.py
=============================
Continuation of the review response -- the Tier-2 item REVISION_PROPOSAL.md left
open:

    "Tie E to real geometry ... Then L(G) would model anisotropic diffusion on a
     surface, a real and useful PDE class. Needed: a definition mapping
     geometry -> E, and a convergence test against a problem with a known
     anisotropic answer."

fix03 has just conceded Tier 3: the anisotropic operator is, identically,
anisotropic P1 FEM (no new operator). So this script does NOT try to claim
novelty. It answers the strictly weaker, still-worth-settling question:

    With a GEOMETRICALLY DEFINED anisotropy E, is L(G) a *consistent*
    discretisation of a genuine anisotropic surface energy -- i.e. does the
    discrete energy converge to the exact continuous one at O(h^2)?

Two findings, both with numbers:

  PART 1 (a clean negative).  The most natural geometric source of anisotropy is
  curvature (the shape operator / second fundamental form). On the unit sphere
  every point is UMBILIC: principal curvatures are equal (k1 = k2 = 1), so the
  *deviatoric* (trace-free) part of the shape operator is identically zero. Hence
  curvature CANNOT supply the trace-free anisotropy on S^2 -- E must come from an
  external director field. We confirm the umbilicity numerically with a discrete
  shape-operator estimate (k1 - k2 -> 0).

  PART 2 (a clean positive, but modest).  Take an external tangential director
  field d (here the tangential projection of the global x-axis), build the
  geometric trace-free tensor E = d d^T - (1/2)(I - n n^T), set
  G = (I - n n^T) + eps E, and test the discrete anisotropic Dirichlet energy
  u^T L(G) u against the EXACT continuous energy a(u) = \int_S grad u . G . grad u
  computed by high-order spherical quadrature, for u = x*y. The discrete energy
  converges to the exact value at O(h^2).

Honest scope: O(h^2) *consistency* is exactly what anisotropic P1 FEM already
delivers (fix03). PART 2 therefore confirms the operator is a sound anisotropic-
diffusion discretisation -- it does NOT show any advantage over anisotropic FEM.
That remains open (Tier 3 burden, now conceded as unmet).
r"""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))

from ksf.mesh import icosphere, mean_edge_length, vertex_normals
from fix03_aniso_is_fem import route_A, E_director


# --------------------------------------------------------------------------- #
#  PART 1 : the sphere is umbilic -> deviatoric shape operator = 0             #
# --------------------------------------------------------------------------- #
def deviatoric_shape_ratio(V, F, exact_normals):
    """Mean over vertices of  |k1 - k2| / |k1 + k2|  (dimensionless anisotropy of
    the shape operator), estimated by least-squares fit of the Gauss map over the
    one-ring.

    exact_normals=True uses the exact unit-sphere Gauss map N(p)=p (so the fit
    sees the true shape operator); =False uses angle-weighted discrete normals (a
    generic estimator that does not know it is on a sphere). The first isolates
    the true geometry; the second exposes only estimator noise.
    """
    if exact_normals:
        N = V / np.linalg.norm(V, axis=1, keepdims=True)
    else:
        N = vertex_normals(V, F)
    nbr = [set() for _ in range(len(V))]
    for a, b, c in F:
        for u, v in ((a, b), (b, c), (c, a)):
            nbr[u].add(v); nbr[v].add(u)
    ratios = []
    for i in range(len(V)):
        n = N[i]
        a = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
        t1 = np.cross(n, a); t1 /= np.linalg.norm(t1); t2 = np.cross(n, t1)
        T = np.column_stack([t1, t2])
        X, Y = [], []
        for j in nbr[i]:
            dx = V[j] - V[i]; X.append(T.T @ (dx - np.dot(dx, n) * n))
            dn = N[j] - N[i]; Y.append(T.T @ (dn - np.dot(dn, n) * n))
        if len(X) < 3:
            continue
        X = np.array(X); Y = np.array(Y)
        S, *_ = np.linalg.lstsq(X, Y, rcond=None); S = 0.5 * (S + S.T)
        ev = np.linalg.eigvalsh(S); tr = ev.sum()
        if abs(tr) > 1e-12:
            ratios.append(abs(ev[0] - ev[1]) / abs(tr))
    return float(np.mean(ratios))


# --------------------------------------------------------------------------- #
#  PART 2 : exact continuous anisotropic energy by spherical quadrature        #
# --------------------------------------------------------------------------- #
def exact_energy(eps, n_theta=200, n_phi=400):
    r"""a(u) = \int_{S^2} (grad_S u) . G . (grad_S u) dA   for u = x*y,
    with G = (I - n n^T) + eps E, E from the same global-x director as E_director.
    Gauss-Legendre in cos(theta), trapezoid in phi (spectral for this smooth,
    low-degree integrand) -> reference accurate to ~machine precision."""
    xg, wg = np.polynomial.legendre.leggauss(n_theta)          # nodes in [-1,1]=cos th
    phs = np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False)
    dphi = 2 * np.pi / n_phi
    a_axis = np.array([1.0, 0.0, 0.0])
    total = 0.0
    for ct, wt in zip(xg, wg):
        st = np.sqrt(max(0.0, 1.0 - ct * ct))
        for ph in phs:
            p = np.array([st * np.cos(ph), st * np.sin(ph), ct])
            n = p                                              # unit sphere normal
            P = np.eye(3) - np.outer(n, n)
            d = P @ a_axis; nd = np.linalg.norm(d)
            E = np.zeros((3, 3)) if nd < 1e-12 else \
                (np.outer(d / nd, d / nd) - 0.5 * P)
            G = P + eps * E
            g = P @ np.array([p[1], p[0], 0.0])                # grad_S(xy)
            f = float(g @ G @ g)
            total += wt * dphi * f                             # dA = d(cos th) dphi
    return total


def run():
    out = {}
    # PART 1
    p1 = []
    for lvl in (2, 3, 4, 5):
        V, F = icosphere(lvl)
        p1.append({"level": lvl, "h": mean_edge_length(V, F),
                   "dev_exact": deviatoric_shape_ratio(V, F, True),
                   "dev_approx": deviatoric_shape_ratio(V, F, False)})
    out["umbilic"] = p1

    # PART 2 : exact references
    a_exact = {eps: exact_energy(eps) for eps in (0.0, 0.3)}
    out["a_exact"] = a_exact
    out["a_exact_iso_analytic"] = 8.0 * np.pi / 5.0            # = 6 * \int (xy)^2 dA

    p2 = {}
    for eps in (0.0, 0.3):
        rows = []; prev = None
        for lvl in (2, 3, 4, 5):
            V, F = icosphere(lvl)
            L = route_A(V, F, E_director, eps)
            u = (V[:, 0] * V[:, 1]).astype(float)
            Eh = float(u @ (L @ u))
            err = abs(Eh - a_exact[eps]) / abs(a_exact[eps])
            h = mean_edge_length(V, F)
            rows.append({"level": lvl, "h": h, "E_h": Eh, "rel_err": err})
        # fitted order
        hs = np.array([r["h"] for r in rows]); es = np.array([r["rel_err"] for r in rows])
        order = float(np.polyfit(np.log(hs), np.log(es), 1)[0])
        p2[eps] = {"rows": rows, "order": order}
    out["energy"] = p2
    return out


def main():
    res = run()
    print("=== FIX 02: geometric anisotropy -- can E come from geometry, and is "
          "L(G) a consistent anisotropic operator? ===\n")

    print("PART 1 -- the sphere is UMBILIC, so curvature gives NO trace-free part")
    print("   anisotropy of the shape operator,  |k1 - k2| / |k1 + k2|  :")
    print(f'   {"level":>5} {"h":>8} {"exact Gauss map":>18} {"approx normals":>16}')
    for r in res["umbilic"]:
        print(f'   {r["level"]:5d} {r["h"]:8.4f} {r["dev_exact"]:18.3e} '
              f'{r["dev_approx"]:16.3e}')
    print("   exact column ~ 1e-16 at every level: the true shape operator of S^2")
    print("      is the identity (k1=k2=1) -> deviatoric part is EXACTLY zero.")
    print("   approx column is just the estimator's own discretisation noise and")
    print("      shrinks under refinement -- not a real geometric anisotropy.")
    print("   => curvature supplies no trace-free anisotropy on S^2; E MUST come")
    print("      from an external director field.\n")

    print("PART 2 -- discrete anisotropic energy u^T L(G) u  vs  exact  a(u)")
    print(f'   exact a(u), eps=0   (analytic 8*pi/5 = {res["a_exact_iso_analytic"]:.6f})'
          f' : quadrature {res["a_exact"][0.0]:.6f}')
    print(f'   exact a(u), eps=0.3 : quadrature {res["a_exact"][0.3]:.6f}\n')
    for eps in (0.0, 0.3):
        p = res["energy"][eps]
        print(f'   eps = {eps}')
        print(f'      {"level":>5} {"h":>8} {"E_h":>12} {"rel.err":>12}')
        for r in p["rows"]:
            print(f'      {r["level"]:5d} {r["h"]:8.4f} {r["E_h"]:12.6f} '
                  f'{r["rel_err"]:12.3e}')
        print(f'      fitted convergence order ~ O(h^{p["order"]:.2f})\n')

    print("[結句] PART 1 is a clean, correct negative: the natural geometric")
    print("   anisotropy (curvature) is identically degenerate on the sphere, so")
    print("   the trace-free tensor's content here is an EXTERNAL director, not an")
    print("   intrinsic geometric invariant of S^2. PART 2 is a clean positive but")
    print("   deliberately modest: with such a director, L(G) is a CONSISTENT")
    print("   discretisation of a real anisotropic surface energy, converging at")
    print("   O(h^2) to the exact value (and reproducing the analytic 8*pi/5 at")
    print("   eps=0). Combined with fix03, the honest standing of the whole")
    print("   programme is now: the construction yields a genuine, well-posed,")
    print("   O(h^2)-consistent anisotropic-diffusion operator -- which is exactly")
    print("   anisotropic P1 FEM, no more. Any remaining claim (that the layered /")
    print("   trace-free *formalism* organises this better than FEEC / material")
    print("   Hodge stars) is unproven and must be demonstrated, not asserted.")
    return res


if __name__ == "__main__":
    main()
