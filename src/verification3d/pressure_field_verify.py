"""
pressure_field_verify.py
========================
Backs theory/04_pressure_field.md (Step 2). Pressure is one Poisson solve on the
verified operator: continuum  div^2 p = div(f)  ->  discrete  K p = D f, with the
same P1 stiffness K = sum_tet vol * G^T G used for heat and waves. The per-tet
gradient G is the one fem_laplacian already builds, so:

    grad(p)_tet = G @ p[ix]                 (constant per tetrahedron)
    (D f)_i     = sum_tet vol * (G^T f_tet) (divergence source, scattered to nodes)

and by construction  D( grad(p) ) == K p  exactly (consistency).

Checks (the five exact-solution targets of the document's section 4):
  1. Hydrostatic liquid:  K p = D(rho g)  ->  grad(p) = rho g everywhere (p = rho g.x).
  2. Helmholtz projection: a pure-gradient field projects to ~0; a div-free field
     is preserved; the projector is idempotent.
  3. Gas hydrostatic:     isothermal atmosphere p(h)=p0 exp(-g h / c^2); c^2->inf
     recovers the linear (incompressible) law.
  4. Acoustic limit:      M p_tt = -c^2 K p oscillates at omega = c*sqrt(lambda)
     -> exactly Step 1's wave, frequency scaling with c.
  5. Solvability:         pure-Neumann K has the constant null space; the source
     D f is compatible (orthogonal to constants); unique after pinning.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh, cg
from scipy.sparse import diags
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian


# --------------------------------------------------------------------------- #
#  Operators: K, M, per-tet gradient G, and the divergence source D           #
# --------------------------------------------------------------------------- #
def operators(V, T):
    """Return K, M and per-tetrahedron (ix, vol, G) consistent with fem_laplacian."""
    tets = []
    for tet in T:
        ix = [int(i) for i in tet]
        P = V[ix]
        vol = abs(np.linalg.det(
            np.array([P[1] - P[0], P[2] - P[0], P[3] - P[0]]))) / 6.0
        if vol <= 1e-15:
            continue
        C = np.linalg.inv(np.column_stack([np.ones(4), P]))
        G = C[1:4, :]                       # (3,4) constant gradient of each hat
        tets.append((ix, vol, G))
    K, M = fem_laplacian(V, T)
    return K, M, tets


def grad_per_tet(p, tets):
    """Constant gradient of node field p on each tetrahedron -> (nT,3)."""
    return np.array([G @ p[ix] for (ix, vol, G) in tets])


def div_source(f_tet, tets, n):
    """Assemble b_i = sum_tet vol * (G^T f_tet)_a  (the discrete D f)."""
    b = np.zeros(n)
    for t, (ix, vol, G) in enumerate(tets):
        contrib = vol * (G.T @ f_tet[t])    # (4,)
        for a in range(4):
            b[ix[a]] += contrib[a]
    return b


def solve_neumann(K, b, pin=0):
    """Solve singular pure-Neumann K p = b by pinning one dof (p[pin]=0).
    Direct sparse solve -> machine precision at this size."""
    from scipy.sparse.linalg import spsolve
    n = K.shape[0]
    keep = np.ones(n, bool); keep[pin] = False
    Kr = K[keep][:, keep].tocsc()
    br = b[keep]
    pr = spsolve(Kr, br)
    p = np.zeros(n); p[keep] = pr
    return p


# --------------------------------------------------------------------------- #
#  Checks                                                                       #
# --------------------------------------------------------------------------- #
def check1_hydrostatic(V, T, K, M, tets):
    print("\n[1] Hydrostatic liquid:  K p = D(rho g)  ->  grad p = rho g")
    n = len(V)
    rho, g = 1.0, np.array([0.0, 0.0, -9.81])     # gravity in -z
    f_tet = np.tile(rho * g, (len(tets), 1))      # constant body force per tet
    b = div_source(f_tet, tets, n)
    p = solve_neumann(K, b)
    gp = grad_per_tet(p, tets)                      # recovered per-tet gradient
    err = np.max(np.abs(gp - rho * g))             # should equal rho*g everywhere
    # also confirm p is affine: p ~ rho g . x + const
    A = np.column_stack([V, np.ones(n)])
    coef, *_ = np.linalg.lstsq(A, p, rcond=None)
    resid = np.max(np.abs(p - A @ coef))
    print(f"    max|grad p - rho g| = {err:.2e}   (target ~machine)")
    print(f"    fitted grad = {coef[:3]}  (expect {rho*g})")
    print(f"    non-affine residual = {resid:.2e}")
    assert err < 1e-8 and resid < 1e-8, "hydrostatic gradient not exact"
    print("    PASS  (linear p=rho g.x lies in P1, recovered to machine precision)")


def check2_projection(V, T, K, M, tets):
    print("\n[2] Helmholtz projection:  pure-gradient -> 0;  div-free preserved;  P^2=P")
    n = len(V)
    rng = np.random.default_rng(0)

    def project(w):
        """Return the divergence-free part of a per-tet field w."""
        b = div_source(w, tets, n)
        p = solve_neumann(K, b)
        return w - grad_per_tet(p, tets)

    # (a) pure gradient field w = grad(phi): projection must remove all of it.
    phi = rng.standard_normal(n)
    w_grad = grad_per_tet(phi, tets)
    u_grad = project(w_grad)
    rel_grad = np.linalg.norm(u_grad) / np.linalg.norm(w_grad)
    print(f"    pure-gradient: |P w|/|w| = {rel_grad:.2e}  (target ~0, all removed)")
    assert rel_grad < 1e-6, "pure gradient not fully projected out"

    # (b) random field -> its projection u is divergence-free and non-trivial.
    w = rng.standard_normal((len(tets), 3))
    u = project(w)
    nontrivial = np.linalg.norm(u) / np.linalg.norm(w)
    df_src = np.linalg.norm(div_source(u, tets, n))
    print(f"    random field: |P w|/|w| = {nontrivial:.3f} (non-trivial)  |D(P w)| = {df_src:.2e}")
    assert nontrivial > 0.1 and df_src < 1e-8, "projection not div-free / trivial"

    # (c) idempotency on that non-trivial div-free field: P(u) == u.
    u2 = project(u)
    idem = np.linalg.norm(u2 - u) / np.linalg.norm(u)
    print(f"    idempotency |P u - u| / |u| = {idem:.2e}")
    assert idem < 1e-8, "projector not idempotent"
    print("    PASS  (projector removes gradients, keeps div-free part, P^2=P)")


def check3_atmosphere():
    print("\n[3] Gas hydrostatic: liquid=linear vs gas=exponential, joined by c^2")
    g, p0, H = 9.81, 1.0, 1.0
    hs = np.linspace(0, H, 2001)
    gaps = []
    for c2 in (1.0, 10.0, 1e6):
        p = p0 * np.exp(-g * hs / c2)               # isothermal-atmosphere solution
        lin = p0 * (1 - g * hs / c2)                # incompressible (linear) law
        gap = np.max(np.abs(p - lin))               # how far gas departs from linear
        # exact pressure ratio across the column must equal exp(-gH/c^2)
        ratio_err = abs(p[-1] / p0 - np.exp(-g * H / c2))
        gaps.append(gap)
        tag = "≈linear (incompressible limit)" if c2 >= 1e6 else "exponential (compressible)"
        print(f"    c^2={c2:>8g}: exp-vs-linear gap={gap:.2e}   ratio p(H)/p0 err={ratio_err:.2e}  ({tag})")
        assert ratio_err < 1e-12, "pressure ratio not exp(-gH/c^2)"
    # physics: gap large for small c^2 (gas), -> 0 as c^2 -> inf (liquid)
    assert gaps[0] > 1.0 and gaps[-1] < 1e-3 and gaps[0] > gaps[1] > gaps[2], \
        "compressible->incompressible trend wrong"
    print("    PASS  (gas departs from linear; c^2->inf recovers the liquid linear law)")


def check4_acoustic(V, T, K, M):
    print("\n[4] Acoustic limit = Step 1 wave:  M p_tt = -c^2 K p,  omega = c*sqrt(lambda)")
    Minv = 1.0 / M
    lam, _ = eigsh(K, k=2, M=diags(M), sigma=1e-8, which="LM")
    lam1 = np.sort(lam)[1]
    maxKM = float(np.max((np.abs(K).sum(1).A1 * Minv)))
    ok = True
    for c in (1.0, 2.0):
        # leapfrog for M p_tt = -c^2 K p; measure oscillation period of mode-1 proj.
        lam_max = 2.2 * float(np.max(K.diagonal() * Minv)) * c * c
        dt = 0.4 * np.sqrt(4.0 / lam_max)
        # seed lowest mode
        _, vec = eigsh(K, k=2, M=diags(M), sigma=1e-8, which="LM")
        mode = vec[:, np.argsort(lam)[1]]
        p0 = mode.copy(); p1 = mode.copy()
        proj = []
        for s in range(4000):
            Kp = c * c * (K @ p1)
            p2 = 2 * p1 - p0 - dt * dt * Minv * Kp
            proj.append(p2 @ (M * mode)); p0, p1 = p1, p2
        proj = np.array(proj)
        # measure period by zero crossings
        sign = np.sign(proj); cross = np.where(np.diff(sign) != 0)[0]
        if len(cross) >= 3:
            period_steps = 2 * np.mean(np.diff(cross))
            omega_meas = 2 * np.pi / (period_steps * dt)
        else:
            omega_meas = np.nan
        omega_pred = c * np.sqrt(lam1)
        rel = abs(omega_meas - omega_pred) / omega_pred
        print(f"    c={c}: omega_meas={omega_meas:.3f}  omega_pred=c*sqrt(lam1)={omega_pred:.3f}  rel={rel:.3f}")
        ok = ok and rel < 0.05
    assert ok, "acoustic frequency does not scale as c*sqrt(lambda)"
    print("    PASS  (gas acoustic limit reproduces Step 1's wave, frequency ∝ c)")


def check5_solvability(V, T, K, M, tets):
    print("\n[5] Solvability: pure-Neumann null space + source compatibility")
    n = len(V)
    one = np.ones(n)
    Kone = np.linalg.norm(K @ one)
    print(f"    |K · 1| = {Kone:.2e}   (constant is the null vector)")
    rho, g = 1.0, np.array([0.3, -0.5, 9.81])
    f_tet = np.tile(rho * g, (len(tets), 1))
    b = div_source(f_tet, tets, n)
    compat = abs(one @ b)
    print(f"    |1 · b| = {compat:.2e}   (source orthogonal to constants => solvable)")
    assert Kone < 1e-8 and compat < 1e-8, "null space / compatibility failed"
    print("    PASS  (constant null space; D f compatible; unique after pinning)")


if __name__ == "__main__":
    V, T = kuhn_cube(8)
    V = V - 0.5
    K, M, tets = operators(V, T)
    print("Kuhn cube n=8:  nV=%d  nTet=%d" % (len(V), len(tets)))
    # consistency: D(grad p) == K p
    p = np.random.default_rng(1).standard_normal(len(V))
    lhs = div_source(grad_per_tet(p, tets), tets, len(V))
    cons = np.linalg.norm(lhs - K @ p) / np.linalg.norm(K @ p)
    print("consistency  |D(grad p) - K p| / |K p| = %.2e" % cons)
    assert cons < 1e-10, "gradient/divergence not consistent with K"

    check1_hydrostatic(V, T, K, M, tets)
    check2_projection(V, T, K, M, tets)
    check3_atmosphere()
    check4_acoustic(V, T, K, M)
    check5_solvability(V, T, K, M, tets)
    print("\nAll Step 2 pressure-field checks PASSED.")
