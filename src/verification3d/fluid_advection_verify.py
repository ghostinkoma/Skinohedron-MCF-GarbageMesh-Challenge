"""
fluid_advection_verify.py  --  V3 Stage B (scalar advection-diffusion).

Backs theory/06b_advection.md. Isolates the new convective operator C(u) by carrying
a passive scalar in a PRESCRIBED velocity, so the problem is linear in phi and C(u)
can be checked against exact solutions before the nonlinear coupling of Stage C.

Operator: on each tet the P1 gradient G gives u.grad phi = u.(G phi_tet); assembled
and lumped through M this is C(u). The skew form C_skew = (C - C^T)/2 conserves the
scalar energy 1/2 phi^T M phi (= u.grad when div u = 0), the stable baseline.

Checks:
  A  Consistency:  C . 1 = 0           (a constant advects to nothing)   [machine]
  B  Skew-symmetry: phi^T C_skew phi=0 (advection conserves energy)      [machine]
  C  Mass conservation: int phi const  under advection                   [machine]
  D  Pure advection (uniform u, periodic): centre of mass moves at u     [FE ~5%],
     energy bounded (stability of the skew form)
  E  Advection-diffusion Gaussian: variance grows at 2*kappa*t [~1-2%],
     centre of mass moves at u                                           [FE ~5%]
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import coo_matrix
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian


def per_tet(V, T):
    tets = []
    for tet in T:
        ix = [int(i) for i in tet]
        P = V[ix]
        vol = abs(np.linalg.det(np.array([P[1]-P[0], P[2]-P[0], P[3]-P[0]]))) / 6.0
        if vol <= 1e-15:
            continue
        C = np.linalg.inv(np.column_stack([np.ones(4), P]))
        tets.append((ix, vol, C[1:4, :]))
    return tets


def build_C(u0, tets, n):
    """Convective matrix for constant velocity u0: (C phi)_a = sum_tet vol/4 (u0 . G phi)."""
    rows, cols, data = [], [], []
    for (ix, vol, G) in tets:
        ug = u0 @ G                       # (4,)
        for a in range(4):
            for b in range(4):
                rows.append(ix[a]); cols.append(ix[b]); data.append(vol/4.0 * ug[b])
    return coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()


def remap_periodic(V, axes, L=1.0):
    n = len(V); remap = np.arange(n)
    for ax in axes:
        for ip in np.where(np.isclose(V[:, ax], L))[0]:
            q = V[ip].copy(); q[ax] = 0.0
            j = np.where(np.all(np.isclose(V, q), axis=1))[0]
            if len(j):
                remap[ip] = j[0]
    for i in range(n):
        r = i
        while remap[r] != r:
            r = remap[r]
        remap[i] = r
    _, inv = np.unique(remap, return_inverse=True)
    return inv, int(inv.max()+1)


def fold(A, inv, nn):
    Ac = A.tocoo()
    return coo_matrix((Ac.data, (inv[Ac.row], inv[Ac.col])), shape=(nn, nn)).tocsr()


# --------------------------------------------------------------------------- #
def checkA_consistency(C, n):
    print("\n[A] Consistency: C . 1 = 0 (a constant advects to nothing)")
    r = np.linalg.norm(C @ np.ones(n))
    print(f"    ||C . 1|| = {r:.2e}")
    assert r < 1e-10, "constant is not advected to zero"
    print("    PASS")


def checkB_skew(C, M, n):
    print("\n[B] Skew-symmetry: phi^T C_skew phi = 0 (advection conserves energy)")
    Cs = (C - C.T) * 0.5
    rng = np.random.default_rng(0)
    phi = rng.standard_normal(n)
    e = phi @ (Cs @ phi)
    print(f"    phi^T C_skew phi = {e:.2e}")
    assert abs(e) < 1e-10, "skew form not energy-neutral"
    print("    PASS")


def checkC_mass(Cp, Mp, Minv, nn):
    print("\n[C] Mass conservation: int phi constant under advection (periodic)")
    Cs = (Cp - Cp.T) * 0.5
    rng = np.random.default_rng(1)
    phi = np.abs(rng.standard_normal(nn)) + 0.1
    tot0 = Mp @ phi
    dt = 2e-4
    for _ in range(500):
        phi = phi - dt*Minv*(Cs @ phi)
    tot1 = Mp @ phi
    rel = abs(tot1-tot0)/abs(tot0)
    print(f"    int phi: {tot0:.5f} -> {tot1:.5f}  rel change {rel:.2e}")
    assert rel < 1e-8, "mass not conserved on periodic domain"
    print("    PASS")


def stats_x(u, Mp, Vr):
    w = Mp*u; tot = w.sum()
    ang = 2*np.pi*Vr[:, 0]
    cx = np.arctan2((w*np.sin(ang)).sum(), (w*np.cos(ang)).sum())/(2*np.pi) % 1.0
    var = (w*(((Vr[:, 0]-cx+0.5) % 1.0-0.5)**2)).sum()/tot
    return cx, var, tot


def checkD_pure_advection(Cp, Mp, Minv, Vr, nn, u=1.0):
    print("\n[D] Pure advection (uniform u, periodic): centre moves at u; energy bounded")
    x0, s0 = 0.3, 0.10
    phi = np.exp(-(Vr[:, 0]-x0)**2/(2*s0**2))
    Cs = (Cp - Cp.T)*0.5
    dt, tend = 4e-4, 0.2
    e0 = np.sqrt(phi @ (Mp*phi)); emax = e0
    cx0, _, _ = stats_x(phi, Mp, Vr)
    u_f = phi.copy()
    for _ in range(int(tend/dt)):
        u_f = u_f - dt*Minv*(Cs @ u_f)
        emax = max(emax, np.sqrt(u_f @ (Mp*u_f)))
    cx, _, _ = stats_x(u_f, Mp, Vr)
    speed = ((cx - cx0 + 0.5) % 1.0 - 0.5)/tend
    print(f"    centre-of-mass speed = {speed:.3f}  vs u = {u:.3f}  (FE dispersion)")
    print(f"    energy ratio max/init = {emax/e0:.3f}  (bounded => stable skew form)")
    assert abs(speed-u)/u < 0.08 and emax/e0 < 1.05, "advection speed/stability off"
    print("    PASS")


def checkE_gaussian_addiff(Cp, Kp, Mp, Minv, Vr, nn, u=1.0, kappa=0.003):
    print("\n[E] Advection-diffusion Gaussian: variance ~ 2*kappa*t; centre at u")
    x0, s0 = 0.3, 0.08
    phi = np.exp(-(Vr[:, 0]-x0)**2/(2*s0**2))
    Cs = (Cp - Cp.T)*0.5
    dt, tend = 4e-4, 0.2
    cx0, var0, _ = stats_x(phi, Mp, Vr)
    u_f = phi.copy()
    for _ in range(int(tend/dt)):
        u_f = u_f - dt*Minv*(Cs @ u_f) + dt*kappa*Minv*(-(Kp @ u_f))
    cx, var, _ = stats_x(u_f, Mp, Vr)
    var_pred = s0**2 + 2*kappa*tend
    speed = ((cx - cx0 + 0.5) % 1.0 - 0.5)/tend
    print(f"    variance: {var:.5f} vs exact s0^2+2kt = {var_pred:.5f}  err {abs(var-var_pred):.2e}")
    print(f"    centre speed = {speed:.3f} vs u = {u:.3f}")
    assert abs(var-var_pred) < 0.02*var_pred and abs(speed-u)/u < 0.08, "advection-diffusion off"
    print("    PASS")


if __name__ == "__main__":
    n = 12
    V, T = kuhn_cube(n)
    K, M = fem_laplacian(V, T)
    tets = per_tet(V, T)
    nV = len(V)
    u0 = np.array([1.0, 0.0, 0.0])
    C = build_C(u0, tets, nV)
    print("Advection box n=%d: nV=%d nTet=%d" % (n, nV, len(T)))

    # algebraic operator checks (A,B use the raw assembled C)
    Minv = 1.0/M
    checkA_consistency(C, nV)
    checkB_skew(C, M, nV)

    # periodic transport checks (C,D,E — closed domain, no boundary flux)
    inv, nn = remap_periodic(V, axes=(0, 1, 2))
    Cp = fold(C, inv, nn); Kp = fold(K, inv, nn)
    Mp = np.zeros(nn)
    for i, m in enumerate(M):
        Mp[inv[i]] += m
    Minv_p = 1.0/Mp
    Vr = np.zeros((nn, 3))
    for i, p in enumerate(V):
        Vr[inv[i]] = np.where(np.isclose(p, 1.0), 0.0, p)
    checkC_mass(Cp, Mp, Minv_p, nn)
    checkD_pure_advection(Cp, Mp, Minv_p, Vr, nn)
    checkE_gaussian_addiff(Cp, Kp, Mp, Minv_p, Vr, nn)
    print("\nAll V3 Stage B (advection-diffusion) checks PASSED.")
