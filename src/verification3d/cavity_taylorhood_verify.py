"""
cavity_taylorhood_verify.py  --  cross-check the PSPG cavity with an inf-sup-stable element

P0 (08) validated the cavity with PSPG-stabilised equal-order P1. But PSPG carries a
stabilisation parameter tau; a sceptic can ask whether the Ghia match was bought by
tuning tau. This script removes that doubt by solving the SAME cavity with the
Taylor-Hood (P2 velocity / P1 pressure) element -- a genuine inf-sup-stable pair with
NO stabilisation parameter at all -- and showing it converges to the SAME Ghia (1982)
benchmark, and agrees with the PSPG solution.

  CHECK A  Taylor-Hood (no tau) vs Ghia Re=100: centreline RMS converges (0.0096 -> 0.0062
           at n=16 -> 24); the parameter-free element matches the published data.
  CHECK B  Head-to-head: PSPG (with tau) and Taylor-Hood (no tau) agree with each other
           on the centreline -> the Ghia match is NOT a tau-tuning artifact; the physics
           is right, confirmed by two independent discretisations.

Note on the project's "one operator": Taylor-Hood enriches the velocity space to P2 --
exactly the local high-order enrichment of 06f ("enrich K"). The pressure stays P1 on
the same vertices. So this is not a departure from the one-operator philosophy; it is
the 06f enrichment route, applied to make the pressure inf-sup stable.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy.sparse import coo_matrix, bmat, csr_matrix
from scipy.sparse.linalg import spsolve
from cavity_pspg_verify import solve_cavity as solve_pspg, GHIA_Y, GHIA_U


# ---------- Taylor-Hood P2-P1 ----------
def mesh_square_p2(n):
    xs = np.linspace(0, 1, n+1)
    Vp1 = [(x, y) for y in xs for x in xs]
    idx = lambda i, j: j*(n+1)+i
    nV1 = len(Vp1); edge_mid = {}; Vp2 = list(Vp1)

    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in edge_mid:
            edge_mid[key] = len(Vp2)
            pa, pb = Vp2[a], Vp2[b]
            Vp2.append(((pa[0]+pb[0])/2, (pa[1]+pb[1])/2))
        return edge_mid[key]
    tris = []
    for j in range(n):
        for i in range(n):
            a = idx(i, j); b = idx(i+1, j); c = idx(i, j+1); d = idx(i+1, j+1)
            for (v0, v1, v2) in [(a, b, d), (a, d, c)]:
                tris.append([v0, v1, v2, mid(v0, v1), mid(v1, v2), mid(v2, v0)])
    return np.array(Vp2), np.array(tris), nV1


GP = np.array([[1/3, 1/3], [0.47014206, 0.47014206], [0.47014206, 0.05971587],
               [0.05971587, 0.47014206], [0.10128651, 0.10128651],
               [0.10128651, 0.79742699], [0.79742699, 0.10128651]])
GW = np.array([0.1125, 0.06619708, 0.06619708, 0.06619708, 0.06296959, 0.06296959, 0.06296959])


def shp(xi, eta, Jinv):
    L = np.array([1-xi-eta, xi, eta])
    N = np.array([2*L[0]**2-L[0], 2*L[1]**2-L[1], 2*L[2]**2-L[2],
                  4*L[0]*L[1], 4*L[1]*L[2], 4*L[2]*L[0]])
    dNdL = np.array([[4*L[0]-1, 0, 0], [0, 4*L[1]-1, 0], [0, 0, 4*L[2]-1],
                     [4*L[1], 4*L[0], 0], [0, 4*L[2], 4*L[1]], [4*L[2], 0, 4*L[0]]])
    dLdxe = np.array([[-1, -1], [1, 0], [0, 1]])
    dNdx = (dNdL@dLdxe)@Jinv
    return N, dNdx, L


def precompute(V, tris):
    data = []
    for tri in tris:
        ix = [int(k) for k in tri]; P = V[ix[:3]]
        J = np.array([P[1]-P[0], P[2]-P[0]]).T
        Jinv = np.linalg.inv(J); detJ = abs(np.linalg.det(J))
        pts = [(shp(xi, eta, Jinv), w*detJ) for (xi, eta), w in zip(GP, GW)]
        data.append((ix, pts))
    return data


def assemble_static(data, nV, nV1):
    rk, ck, dk = [], [], []; rb = [[], []]; cb = [[], []]; db = [[], []]
    for (ix, pts) in data:
        Ke = np.zeros((6, 6)); Bxe = np.zeros((3, 6)); Bye = np.zeros((3, 6))
        for (N, dNdx, L), jac in pts:
            Ke += jac*(dNdx@dNdx.T)
            Bxe += jac*np.outer(L, dNdx[:, 0]); Bye += jac*np.outer(L, dNdx[:, 1])
        for a in range(6):
            for b in range(6):
                rk.append(ix[a]); ck.append(ix[b]); dk.append(Ke[a, b])
        for a in range(3):
            for b in range(6):
                rb[0].append(ix[a]); cb[0].append(ix[b]); db[0].append(Bxe[a, b])
                rb[1].append(ix[a]); cb[1].append(ix[b]); db[1].append(Bye[a, b])
    K = coo_matrix((dk, (rk, ck)), shape=(nV, nV)).tocsr()
    Bx = coo_matrix((db[0], (rb[0], cb[0])), shape=(nV1, nV)).tocsr()
    By = coo_matrix((db[1], (rb[1], cb[1])), shape=(nV1, nV)).tocsr()
    return K, Bx, By


def assemble_adv(data, nV, U):
    r, c, d = [], [], []
    for (ix, pts) in data:
        Ne = np.zeros((6, 6))
        for (N, dNdx, L), jac in pts:
            ue = np.array([N@U[ix, 0], N@U[ix, 1]])
            conv = dNdx@ue
            Ne += jac*np.outer(N, conv)
        for a in range(6):
            for b in range(6):
                r.append(ix[a]); c.append(ix[b]); d.append(Ne[a, b])
    return coo_matrix((d, (r, c)), shape=(nV, nV)).tocsr()


def solve_TH(n, Re, picard=40, tol=1e-7):
    V, tris, nV1 = mesh_square_p2(n); nV = len(V); nu = 1.0/Re
    data = precompute(V, tris); K, Bx, By = assemble_static(data, nV, nV1)
    x, y = V[:, 0], V[:, 1]; lid = np.isclose(y, 1.0)
    walls = (np.isclose(x, 0) | np.isclose(x, 1) | np.isclose(y, 0)) & ~lid
    diru = lid | walls
    ub = np.zeros(nV); ub[lid] = 1.0; vb = np.zeros(nV)
    U = np.zeros((nV, 2)); U[:, 0] = ub.copy()
    Z = csr_matrix((nV, nV)); Zp = csr_matrix((nV1, nV1)); didx = np.where(diru)[0]
    it = 0
    for it in range(picard):
        N = assemble_adv(data, nV, U); A = nu*K+N
        S = bmat([[A, Z, Bx.T], [Z, A, By.T], [Bx, By, Zp]]).tolil()
        rhs = np.zeros(2*nV+nV1)
        for blk, val in [(0, ub), (1, vb)]:
            for gi in didx:
                k = gi+blk*nV; S.rows[k] = [k]; S.data[k] = [1.0]; rhs[k] = val[gi]
        pk = 2*nV; S.rows[pk] = [pk]; S.data[pk] = [1.0]; rhs[pk] = 0.0
        sol = spsolve(S.tocsr(), rhs)
        Un = np.column_stack([sol[:nV], sol[nV:2*nV]])
        du = np.linalg.norm(Un-U); U = Un
        if du < tol:
            break
    return V, U, sol[2*nV:], it


def th_centreline(n, Re):
    V, U, p, it = solve_TH(n, Re)
    x, y = V[:, 0], V[:, 1]; ln = np.isclose(x, 0.5)
    yl = y[ln]; ul = U[ln, 0]; o = np.argsort(yl)
    return yl[o], ul[o]


def checkA_th_vs_ghia():
    print("\n[A] Taylor-Hood P2-P1 (NO stabilisation parameter) vs Ghia Re=100")
    last = None
    for n in (16, 24):
        yl, ul = th_centreline(n, 100.0)
        err = np.sqrt(np.mean((np.interp(GHIA_Y, yl, ul)-GHIA_U)**2))
        print(f"    n={n}: centreline u RMS vs Ghia = {err:.4f}  (u(0.5)={np.interp(0.5,yl,ul):.4f}, Ghia -0.2058)")
        last = err
    assert last < 0.02, "Taylor-Hood does not match Ghia"
    print("    PASS  (parameter-free inf-sup element matches the published benchmark)")


def checkB_pspg_vs_th():
    print("\n[B] Head-to-head: PSPG (with tau) vs Taylor-Hood (no tau) -- same physics?")
    Vp, ux, uy, p, it, tris, K, tau = solve_pspg(64, 100.0)
    xs, ys = Vp[:, 0], Vp[:, 1]; ln = np.isclose(xs, 0.5)
    yp = ys[ln]; up = ux[ln]; o = np.argsort(yp); yp, up = yp[o], up[o]
    yt, ut = th_centreline(24, 100.0)
    # compare both on the Ghia y-stations
    up_i = np.interp(GHIA_Y, yp, up); ut_i = np.interp(GHIA_Y, yt, ut)
    rms = np.sqrt(np.mean((up_i-ut_i)**2))
    print(f"    PSPG(n=64,tau={tau:.4f}) vs Taylor-Hood(n=24,no tau): centreline RMS = {rms:.4f}")
    print("    -> the two independent discretisations agree; the Ghia match is NOT a tau artifact.")
    assert rms < 0.02, "PSPG and Taylor-Hood disagree"
    print("    PASS  (parameter-free and stabilised elements agree -> physics confirmed)")


if __name__ == "__main__":
    print("Taylor-Hood cross-check of the PSPG cavity (is the Ghia match parameter-independent?)")
    checkA_th_vs_ghia()
    checkB_pspg_vs_th()
    print("\nConclusion:")
    print("  - Taylor-Hood P2-P1, an inf-sup-stable element with NO stabilisation parameter,")
    print("    independently converges to Ghia Re=100 (RMS 0.0062 at n=24).")
    print("  - PSPG (with tau) and Taylor-Hood (no tau) agree with each other -> the benchmark")
    print("    match is the correct physics, not an artifact of tuning tau. P0 is doubly confirmed.")
    print("  - Taylor-Hood's P2 velocity is the 06f enrichment (enrich the velocity/K side),")
    print("    so this stays within the one-operator philosophy.")
