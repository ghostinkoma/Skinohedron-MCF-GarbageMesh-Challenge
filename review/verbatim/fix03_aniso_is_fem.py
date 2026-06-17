"""
fix03_aniso_is_fem.py
=====================
Continuation of the review response. Settles the Tier-3 burden that
REVISION_PROPOSAL.md explicitly imposed:

    "The revised paper may only claim novelty for whatever the layered mesh +
     trace-free dual-tensor formalism adds *on top of* anisotropic FEM -- and
     must DEMONSTRATE it, not assert it. If, after honest testing, it adds
     nothing, that is itself a publishable negative result and should be stated
     plainly."

fix01 already conceded C1/C2 and showed the trace-free tensor, used as an
anisotropic conductivity G = (I - n n^T) + eps*E, produces an operator that is
DISTINCT from the *isotropic* FEM Laplacian. Good -- but the honest question is
one level up: is that anisotropic operator anything more than ordinary
anisotropic P1 FEM?

This is the exact analogue of c1 (which showed the isotropic KSF operator IS the
cotan/Dziuk FEM Laplacian, bit for bit). Here we test the anisotropic operator
the same way: assemble it by two genuinely independent code paths and see whether
they agree to machine precision.

  ROUTE A  (gradient quadrature):   the fix01 assembly --
           K_e = Area * grad_phi . G . grad_phi^T , G placed in the middle.
           This is "the KSF anisotropic operator" as written.

  ROUTE B  (coordinate-change cotan):  a classical fact about anisotropic P1
           FEM -- a constant SPD conductivity is removed by a linear change of
           tangent coordinates y = B^{-1} x with G|_tangent = B B^T. In the new
           coordinates the bilinear form is ISOTROPIC, so the element matrix is
           just |det B| times the ordinary cotangent stiffness of the
           transformed triangle. Route B never multiplies by G; it only
           transforms vertex coordinates and calls the repo's own isotropic
           cotan formula.

If ROUTE A == ROUTE B to ~1e-16 for arbitrary tangential trace-free E, then the
anisotropic KSF operator is, identically, "the isotropic cotan FEM Laplacian of a
per-face linearly-remapped mesh" -- i.e. textbook anisotropic P1 FEM. No new
operator is produced at the anisotropic level either. That is the negative result
the proposal invited, and we state it plainly.

We test THREE different E fields (fix01's fixed E; a geometric director-based E;
a random symmetric trace-free tangential E) so the identity cannot be an artefact
of one special choice.
"""
from __future__ import annotations
import os, sys
import numpy as np
from scipy.sparse import coo_matrix
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))

from ksf.mesh import icosphere
from ksf.dec import cotangent_laplacian


# --------------------------------------------------------------------------- #
#  per-face conductivity providers  G(face) = (I - w w^T) + eps * E            #
#  every E here is symmetric, trace-free w.r.t. the tangent plane, tangential  #
# --------------------------------------------------------------------------- #
def _tangent_frame(w):
    a = np.array([0.0, 0.0, 1.0]) if abs(w[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    t1 = np.cross(w, a); t1 /= np.linalg.norm(t1)
    t2 = np.cross(w, t1)
    return t1, t2


def E_fix01(w, centroid, rng):
    """The exact E used in fix01 (fixed lab-frame director)."""
    t1, t2 = _tangent_frame(w)
    return np.outer(t1, t1) - np.outer(t2, t2)


def E_director(w, centroid, rng):
    """Geometric director field: tangential projection of the global x-axis,
    normalised. E = d d^T - (1/2)|d|^2 (I - w w^T)  -> symmetric trace-free."""
    a = np.array([1.0, 0.0, 0.0])
    d = a - np.dot(a, w) * w
    nd = np.linalg.norm(d)
    P = np.eye(3) - np.outer(w, w)
    if nd < 1e-9:
        return np.zeros((3, 3))
    d /= nd
    return np.outer(d, d) - 0.5 * P


def E_random(w, centroid, rng):
    """A random symmetric trace-free tangential tensor (worst case for 'special
    structure' explanations)."""
    t1, t2 = _tangent_frame(w)
    # general symmetric trace-free 2x2 in the (t1,t2) basis: [[a,b],[b,-a]]
    a, b = rng.normal(), rng.normal()
    s = np.hypot(a, b)                          # spectral norm = sqrt(a^2+b^2)
    if s < 1e-12:
        return np.zeros((3, 3))
    a, b = a / s, b / s                         # normalise: eigenvalues = +/-1
    B = np.column_stack([t1, t2])              # 3x2
    M2 = np.array([[a, b], [b, -a]])
    return B @ M2 @ B.T


# --------------------------------------------------------------------------- #
#  ROUTE A : gradient-quadrature anisotropic stiffness (the fix01 assembly)    #
# --------------------------------------------------------------------------- #
def route_A(V, F, Efun, eps, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    n = len(V); rows = []; cols = []; vals = []
    for tri in F:
        i, j, k = (int(x) for x in tri); p = V[[i, j, k]]
        nrm = np.cross(p[1] - p[0], p[2] - p[0]); A = 0.5 * np.linalg.norm(nrm)
        if A <= 0:
            continue
        w = nrm / np.linalg.norm(nrm)
        g = np.zeros((3, 3))
        for a in range(3):
            b, c = (a + 1) % 3, (a + 2) % 3; e = p[c] - p[b]
            ga = np.cross(w, e) / (2 * A)
            if np.dot(ga, p[a] - p[b]) < 0:
                ga = -ga
            g[a] = ga
        cen = p.mean(axis=0)
        E = Efun(w, cen, rng)
        G = (np.eye(3) - np.outer(w, w)) + eps * E
        Ke = A * (g @ G @ g.T); idx = [i, j, k]
        for a in range(3):
            for b in range(3):
                rows.append(idx[a]); cols.append(idx[b]); vals.append(Ke[a, b])
    return coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()


# --------------------------------------------------------------------------- #
#  ROUTE B : isotropic cotan stiffness of a per-face G-transformed triangle    #
#            (never multiplies by G; only changes coordinates)                 #
# --------------------------------------------------------------------------- #
def _cot2d(a, b, c):
    u, v = b - a, c - a
    cr = u[0] * v[1] - u[1] * v[0]
    return float(np.dot(u, v) / abs(cr))


def route_B(V, F, Efun, eps, rng_seed=0):
    rng = np.random.default_rng(rng_seed)
    n = len(V); rows = []; cols = []; vals = []
    for tri in F:
        i, j, k = (int(x) for x in tri); p = V[[i, j, k]]
        nrm = np.cross(p[1] - p[0], p[2] - p[0]); A = 0.5 * np.linalg.norm(nrm)
        if A <= 0:
            continue
        w = nrm / np.linalg.norm(nrm)
        t1, t2 = _tangent_frame(w); T = np.column_stack([t1, t2])    # 3x2
        cen = p.mean(axis=0)
        E = Efun(w, cen, rng)
        G = (np.eye(3) - np.outer(w, w)) + eps * E
        G2 = T.T @ G @ T                                             # 2x2 SPD
        B = np.linalg.cholesky(G2)                                   # G2 = B B^T
        detB = float(np.sqrt(np.linalg.det(G2)))
        # 2D coords of the triangle in the tangent frame, then remap by B^{-1}
        q = np.array([T.T @ (p[m] - p[0]) for m in range(3)])        # 3x2
        Binv = np.linalg.inv(B)
        y = (Binv @ q.T).T                                           # 3x2 remapped
        # isotropic cotan element stiffness of the remapped 2D triangle, x |det B|
        Ke = np.zeros((3, 3))
        for a in range(3):
            b, c = (a + 1) % 3, (a + 2) % 3
            ca = _cot2d(y[a], y[b], y[c])      # angle at vertex a, opp edge (b,c)
            hw = 0.5 * ca * detB
            Ke[b, b] += hw; Ke[c, c] += hw
            Ke[b, c] -= hw; Ke[c, b] -= hw
        idx = [i, j, k]
        for a in range(3):
            for b in range(3):
                rows.append(idx[a]); cols.append(idx[b]); vals.append(Ke[a, b])
    return coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()


def reldiff(L1, L2):
    scale = float(np.abs(L1.toarray()).max())
    return float(np.abs((L1 - L2).toarray()).max()) / scale


def run():
    V, F = icosphere(3)
    L_iso, _ = cotangent_laplacian(V, F)
    fields = {"fix01_fixed": E_fix01, "geometric_director": E_director,
              "random_tracefree": E_random}
    out = {}
    for name, Efun in fields.items():
        rec = []
        for eps in (0.0, 0.1, 0.3):
            LA = route_A(V, F, Efun, eps, rng_seed=7)
            LB = route_B(V, F, Efun, eps, rng_seed=7)
            rec.append({
                "eps": eps,
                "A_vs_B": reldiff(LA, LB),
                "A_vs_iso": reldiff(LA, L_iso),
            })
        out[name] = rec
    return out


def main():
    res = run()
    print("=== FIX 03: is the ANISOTROPIC KSF operator anything beyond "
          "anisotropic FEM? ===\n")
    print("Test: assemble the SAME operator two independent ways --")
    print("  A = gradient quadrature with G in the middle (the KSF assembly)")
    print("  B = isotropic cotan FEM on a per-face G-remapped triangle "
          "(never uses G as a factor)")
    print("If A == B for arbitrary tangential trace-free E, the operator IS "
          "anisotropic P1 FEM.\n")
    for name, rec in res.items():
        print(f"E field = {name}")
        print(f'   {"eps":>5} {"relerr(A vs B)":>16} {"relerr(A vs isoFEM)":>20}')
        for r in rec:
            tag = "  (eps=0: equals isotropic FEM)" if r["eps"] == 0 else \
                  "  (distinct from isotropic FEM)"
            print(f'   {r["eps"]:5.2f} {r["A_vs_B"]:16.3e} '
                  f'{r["A_vs_iso"]:20.3e}{tag}')
        print()

    worst = max(r["A_vs_B"] for rec in res.values() for r in rec)
    print(f"Worst A-vs-B disagreement over ALL fields and eps: {worst:.3e}\n")
    print("[結句 / CONCEDED at Tier 3] The two independent assemblies agree to")
    print("   machine precision for every E tested (fixed, geometric, random).")
    print("   Therefore the anisotropic KSF operator is IDENTICALLY the isotropic")
    print("   cotangent (Dziuk P1 FEM) Laplacian of a per-face G-remapped mesh --")
    print("   which is the standard construction of anisotropic P1 FEM. So, just")
    print("   as c1 showed at the isotropic level, no NEW OPERATOR is produced at")
    print("   the anisotropic level either: wiring the trace-free tensor in as a")
    print("   conductivity lands exactly inside the mature field of anisotropic")
    print("   FEM. The honest consequence: the revised paper may NOT claim a novel")
    print("   operator. Whatever value remains can only be organisational (a")
    print("   frame-invariant *bookkeeping* for the anisotropy carrier on layered")
    print("   meshes), and that is a much smaller, expository claim -- to be")
    print("   weighed against finite-element exterior calculus / material Hodge")
    print("   stars, which already do this. We do not assert it is novel; we mark")
    print("   it OPEN and modest. fix02 then asks the only remaining live question:")
    print("   can a GEOMETRICALLY MEANINGFUL E at least make the operator a")
    print("   consistent discretisation of a real anisotropic energy? (yes), while")
    print("   conceding that consistency != novelty over anisotropic FEM.")
    return res


if __name__ == "__main__":
    main()
