"""
fluid_stokes_verify.py  --  V3 Stage A (Stokes flow, linear).

Backs theory/06a_stokes_flow.md. Stokes drops the nonlinear advection term, so the
problem is fully linear and verifiable against exact profiles at machine precision,
exactly as the linear era was. This fixes the plumbing (wall Dirichlet boundaries,
vector momentum, pressure projection) before Stage B's advection.

Channel: unit cube, y in [0,1] is wall-normal (Dirichlet), x,z homogeneous.
  Couette    (wall-driven)     : nu u'' = 0,  u(0)=0, u(1)=U   -> u = U y   (linear)
  Poiseuille (pressure-driven) : nu u'' = G,  u(0)=u(1)=0      -> u = (-G/2nu) y(1-y)

Checks:
  A  Couette steady     -> machine precision vs U y
  B  Poiseuille steady  -> machine precision vs parabola, centreline = -G/(8 nu)
  C  Couette transient  -> converges to steady; slow-mode decay rate ~ nu pi^2
  D  Poiseuille transient -> converges to the parabolic steady state
  E  Incompressibility under the verified projection (the coupling anchor):
     ||D u|| -> machine precision, and idempotency P^2 = P
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian


# --------------------------------------------------------------------------- #
def operators(V, T):
    """K, M and per-tet (ix, vol, G) consistent with fem_laplacian (for D)."""
    tets = []
    for tet in T:
        ix = [int(i) for i in tet]
        P = V[ix]
        vol = abs(np.linalg.det(np.array([P[1]-P[0], P[2]-P[0], P[3]-P[0]]))) / 6.0
        if vol <= 1e-15:
            continue
        C = np.linalg.inv(np.column_stack([np.ones(4), P]))
        tets.append((ix, vol, C[1:4, :]))
    K, M = fem_laplacian(V, T)
    return K, M, tets


def grad_per_tet(p, tets):
    return np.array([G @ p[ix] for (ix, vol, G) in tets])


def div_source(f_tet, tets, n):
    b = np.zeros(n)
    for t, (ix, vol, G) in enumerate(tets):
        c = vol * (G.T @ f_tet[t])
        for a in range(4):
            b[ix[a]] += c[a]
    return b


def solve_dirichlet(K, rhs, bc_val, fixed, n):
    free = ~fixed
    u = np.zeros(n); u[fixed] = bc_val[fixed]
    Kff = K[free][:, free].tocsc()
    b = rhs[free] - (K[free][:, fixed] @ u[fixed])
    u[free] = spsolve(Kff, b)
    return u


def solve_neumann(K, b, n, pin=0):
    keep = np.ones(n, bool); keep[pin] = False
    p = np.zeros(n)
    p[keep] = spsolve(K[keep][:, keep].tocsc(), b[keep])
    return p


# --------------------------------------------------------------------------- #
def checkA_couette_steady(V, K, M, n, wall, y, U=1.0):
    print("\n[A] Couette steady (wall-driven):  u = U y")
    bcv = np.where(np.isclose(y, 1.0), U, 0.0)
    u = solve_dirichlet(K, np.zeros(n), bcv, wall, n)
    err = np.max(np.abs(u - U*y))
    print(f"    max|u - U y| = {err:.2e}   (linear -> P1 exact)")
    assert err < 1e-12, "Couette profile not exact"
    print("    PASS")


def checkB_poiseuille_steady(V, K, M, n, wall, y, nu=1.0, G=-1.0):
    print("\n[B] Poiseuille steady (pressure-driven):  u = (-G/2nu) y(1-y)")
    rhs = M * (-(G/nu)*np.ones(n))            # -lap u = -G/nu
    u = solve_dirichlet(K, rhs, np.zeros(n), wall, n)
    exact = (-G/(2*nu)) * y*(1-y)
    err = np.max(np.abs(u - exact))
    cl = u.max(); cl_exact = -G/(8*nu)
    print(f"    max|u - parabola| = {err:.2e}   centreline {cl:.5f} vs {cl_exact:.5f}")
    assert err < 1e-12 and abs(cl-cl_exact) < 1e-10, "Poiseuille profile not exact"
    print("    PASS")


def checkC_couette_transient(V, K, M, n, wall, y, U=1.0, nu=1.0):
    print("\n[C] Couette transient -> steady; slow-mode decay rate ~ nu pi^2")
    Minv = 1.0/M
    free = ~wall
    bcv = np.where(np.isclose(y, 1.0), U, 0.0)
    u = np.zeros(n); u[wall] = bcv[wall]       # start from rest + wall BC
    steady = U*y
    dt = 0.2/float(np.max(K.diagonal()*Minv))
    norms = []
    for s in range(20000):
        rhs = -nu*(K@u)
        u[free] = u[free] + dt*Minv[free]*rhs[free]
        if s % 50 == 0:
            norms.append(np.sqrt(((u-steady)[free]**2 @ M[free])))
        if norms and norms[-1] < 1e-9:
            break
    conv = norms[-1]
    # fit decay rate from the late-time tail: ||d|| ~ exp(-sigma t)
    tail = np.array(norms[-8:])
    ts = np.arange(len(tail))*50*dt
    sigma = -np.polyfit(ts, np.log(tail), 1)[0]
    sigma_pred = nu*np.pi**2
    rel = abs(sigma-sigma_pred)/sigma_pred
    print(f"    converged ||u-steady|| = {conv:.2e}")
    print(f"    decay rate sigma = {sigma:.3f}  vs nu*pi^2 = {sigma_pred:.3f}  (rel {rel*100:.1f}%)")
    assert conv < 1e-6 and rel < 0.08, "Couette transient did not match"
    print("    PASS")


def checkD_poiseuille_transient(V, K, M, n, wall, y, nu=1.0, G=-1.0):
    print("\n[D] Poiseuille transient -> parabolic steady state")
    Minv = 1.0/M
    free = ~wall
    u = np.zeros(n)                            # rest, walls 0
    steady = (-G/(2*nu)) * y*(1-y)
    body = M*(-(G/nu)*np.ones(n))             # constant forcing (= -lap u_steady)
    dt = 0.2/float(np.max(K.diagonal()*Minv))
    last = 1e9
    for s in range(40000):
        rhs = -nu*(K@u) + nu*body             # nu*(-K u) + forcing; forcing scaled
        u[free] = u[free] + dt*Minv[free]*rhs[free]
        if s % 200 == 0:
            last = np.sqrt(((u-steady)[free]**2 @ M[free]))
            if last < 1e-8:
                break
    print(f"    converged ||u-steady|| = {last:.2e}  centreline {u.max():.5f} vs {-G/(8*nu):.5f}")
    assert last < 1e-5, "Poiseuille transient did not converge to parabola"
    print("    PASS")


def checkE_incompressibility(V, K, M, tets, n, rho=1.0, dt=0.1):
    print("\n[E] Incompressibility under the verified projection (coupling anchor)")
    rng = np.random.default_rng(0)
    # build a non-trivial tentative per-tet velocity with nonzero divergence
    phi = rng.standard_normal(n)
    w = grad_per_tet(phi, tets) + rng.standard_normal((len(tets), 3))*0.3
    # project: K p = (rho/dt) D w ;  u = w - (dt/rho) grad p
    b = (rho/dt) * div_source(w, tets, n)
    p = solve_neumann(K, b, n)
    u = w - (dt/rho)*grad_per_tet(p, tets)
    div_u = np.linalg.norm(div_source(u, tets, n))
    # idempotency: projecting u again changes nothing
    b2 = (rho/dt)*div_source(u, tets, n)
    p2 = solve_neumann(K, b2, n)
    u2 = u - (dt/rho)*grad_per_tet(p2, tets)
    idem = np.linalg.norm(u2 - u)/np.linalg.norm(u)
    print(f"    ||D u|| after projection = {div_u:.2e}   (target ~machine)")
    print(f"    idempotency ||P u - u||/||u|| = {idem:.2e}   (P^2=P)")
    assert div_u < 1e-10 and idem < 1e-10, "projection coupling failed"
    print("    PASS  (vector pressure coupling is the verified projection of 04)")


if __name__ == "__main__":
    V, T = kuhn_cube(8)                        # [0,1]^3
    K, M, tets = operators(V, T)
    n = len(V); y = V[:, 1]
    wall = np.isclose(y, 0.0) | np.isclose(y, 1.0)
    print("Stokes channel n=8:  nV=%d  nTet=%d  wall nodes=%d" % (n, len(T), int(wall.sum())))
    checkA_couette_steady(V, K, M, n, wall, y)
    checkB_poiseuille_steady(V, K, M, n, wall, y)
    checkC_couette_transient(V, K, M, n, wall, y)
    checkD_poiseuille_transient(V, K, M, n, wall, y)
    checkE_incompressibility(V, K, M, tets, n)
    print("\nAll V3 Stage A (Stokes) checks PASSED.")
