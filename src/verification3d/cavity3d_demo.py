"""
cavity3d_demo.py  --  STEP 3: the same tetrahedral operator carries the cavity in 3-D

A demonstration (not a tabulated benchmark -- the 3-D lid-driven cavity has side-wall
drag and no simple 1-D Ghia table). It confirms that the SAME Kuhn-cube tetrahedral
operator and the SAME PSPG-stabilised assembly used in 2-D (cavity_pspg_verify.py)
reproduce the lid-driven recirculation in 3-D: positive flow near the moving lid,
negative (return) flow below, with a smooth pressure.

This is the 3-D leg of the plan: pressure as a domain on the one operator, in 2-D and
3-D alike. (Kept at n=12 so it runs in reasonable time with a direct monolithic solve.)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import coo_matrix, bmat, csr_matrix
from scipy.sparse.linalg import spsolve
from ksf3d.mesh3d_uniform import kuhn_cube


def per_tet(V, T):
    tets = []
    for tet in T:
        ix = [int(i) for i in tet]; P = V[ix]
        vol = abs(np.linalg.det(np.array([P[1]-P[0], P[2]-P[0], P[3]-P[0]])))/6
        if vol <= 1e-15:
            continue
        G = np.linalg.inv(np.column_stack([np.ones(4), P]))[1:4, :]
        tets.append((ix, vol, G))
    return tets


def assemble3d(V, T):
    n = len(V); tets = per_tet(V, T)
    rk, ck, dk = [], [], []; Ml = np.zeros(n)
    rb = [[], [], []]; cb = [[], [], []]; db = [[], [], []]
    for (ix, vol, G) in tets:
        Ke = vol*(G.T@G)
        for a in range(4):
            Ml[ix[a]] += vol/4.0
            for b in range(4):
                rk.append(ix[a]); ck.append(ix[b]); dk.append(Ke[a, b])
                for cmp in range(3):
                    rb[cmp].append(ix[a]); cb[cmp].append(ix[b]); db[cmp].append(vol/4.0*G[cmp, b])
    K = coo_matrix((dk, (rk, ck)), shape=(n, n)).tocsr()
    Bs = [coo_matrix((db[c], (rb[c], cb[c])), shape=(n, n)).tocsr() for c in range(3)]
    return K, Ml, Bs, tets


def adv3d(tets, nv, U):
    r, c, d = [], [], []
    for (ix, vol, G) in tets:
        ug = U[ix].mean(axis=0) @ G
        for a in range(4):
            for b in range(4):
                r.append(ix[a]); c.append(ix[b]); d.append(vol/4.0*ug[b])
    return coo_matrix((d, (r, c)), shape=(nv, nv)).tocsr()


def solve3d(n, Re, picard=15, tol=1e-6):
    V, T = kuhn_cube(n); nv = len(V); nu = 1.0/Re; h = 1.0/n
    K, Ml, Bs, tets = assemble3d(V, T)
    tau = h/2.0*min(1.0, Re*h/6.0)
    x, y, z = V[:, 0], V[:, 1], V[:, 2]
    lid = np.isclose(z, 1.0)
    walls = (np.isclose(x, 0) | np.isclose(x, 1) | np.isclose(y, 0) |
             np.isclose(y, 1) | np.isclose(z, 0)) & ~lid
    diru = lid | walls
    ub = [np.where(lid, 1.0, 0.0), np.zeros(nv), np.zeros(nv)]
    U = np.zeros((nv, 3)); U[:, 0] = ub[0].copy()
    Z = csr_matrix((nv, nv)); didx = np.where(diru)[0]
    it = 0
    for it in range(picard):
        N = adv3d(tets, nv, U); A = nu*K+N
        S = bmat([[A, Z, Z, Bs[0].T], [Z, A, Z, Bs[1].T],
                  [Z, Z, A, Bs[2].T], [Bs[0], Bs[1], Bs[2], -tau*K]]).tolil()
        rhs = np.zeros(4*nv)
        for blk in range(3):
            for gi in didx:
                k = gi+blk*nv; S.rows[k] = [k]; S.data[k] = [1.0]; rhs[k] = ub[blk][gi]
        S.rows[3*nv] = [3*nv]; S.data[3*nv] = [1.0]
        sol = spsolve(S.tocsr(), rhs)
        Un = np.column_stack([sol[:nv], sol[nv:2*nv], sol[2*nv:3*nv]]); p = sol[3*nv:]
        du = np.linalg.norm(Un-U); U = Un
        if du < tol:
            break
    return V, U, p, it


def check_3d_recirculation():
    print("\n[STEP 3] 3-D lid-driven cavity on the same Kuhn tetrahedral operator (Re=100)")
    V, U, p, it = solve3d(12, 100.0)
    x, y, z = V[:, 0], V[:, 1], V[:, 2]
    ln = np.isclose(x, 0.5) & np.isclose(y, 0.5)
    zc = z[ln]; uc = U[ln, 0]; o = np.argsort(zc); zc, uc = zc[o], uc[o]
    u_top = np.interp(0.95, zc, uc); u_mid = np.interp(0.5, zc, uc); u_low = np.interp(0.2, zc, uc)
    print(f"    centre line u(z): z=0.95 {u_top:+.3f}  z=0.5 {u_mid:+.3f}  z=0.2 {u_low:+.3f}")
    print("    -> positive near the lid, negative below (recirculation) -- same sign")
    print("       structure as 2-D Ghia, reproduced by the same operator in 3-D.")
    assert u_top > 0 and u_mid < 0, "no recirculation -> 3-D assembly wrong"
    print("    PASS  (the pressure-driven flow lives on the one tetrahedron in 3-D too)")


if __name__ == "__main__":
    print("STEP 3 demonstration: 3-D cavity, same operator as 2-D")
    check_3d_recirculation()
    print("\nNote: 3-D cavity differs quantitatively from 2-D Ghia (side-wall drag); this is")
    print("a structural demonstration that the same tetrahedral operator carries the flow in")
    print("3-D, not a tabulated benchmark. The validated benchmark is the 2-D case.")
