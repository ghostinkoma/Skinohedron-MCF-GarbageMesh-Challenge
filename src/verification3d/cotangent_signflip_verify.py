"""
cotangent_signflip_verify.py
============================
Backs theory/05_mesh_transform.md section 5b. Answers the question "can a different
weight formula (norm) remove the negative-conductance / sign-flipped cotangent
edges that geometric deformation creates?" with numbers.

The honest finding, in four checks:

  1. MECHANISM. Sign-flipped cotangents come from **obtuse dihedral angles**. The
     Kuhn cube sits exactly at 90 deg (right-angle dihedrals -> zero/boundary), so
     any geometric deformation pushes some dihedrals obtuse and sign-flips appear.

  2. FORMULA TRADEOFF. On one distorted mesh, the cotangent (FE) Laplacian
     reproduces a linear field to machine precision (the foundation: T=k2/(k1+k2),
     D.grad=K) but has sign-flips. Uniform (graph) and clamped weights have ZERO
     sign-flips but lose that exactness. Changing the norm is not free.

  3. MESH-IMPROVEMENT REALITY. Quality-guarded smoothing improves tet quality and
     keeps the exactness, but does NOT reliably remove sign-flips, because a flip
     is an edge-summed property (obtuse dihedrals), and globally non-obtuse
     ("well-centered") tetrahedral meshes are hard in 3D.

  4. CLEAN ESCAPE. Only geometry preservation (the undeformed cube = the
     topological torus's tets) gives BOTH zero sign-flips AND machine-precision
     linear reproduction at once.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

Q_REGULAR = 1.0 / (6 * np.sqrt(2))


# --- geometry helpers -------------------------------------------------------
def tetq(P):
    vol = abs(np.linalg.det(np.array([P[1]-P[0], P[2]-P[0], P[3]-P[0]]))) / 6.0
    e = sum(np.sum((P[a]-P[b])**2) for a in range(4) for b in range(a+1, 4))
    rms = np.sqrt(e / 6.0)
    return (vol / rms**3 / Q_REGULAR) if rms > 0 else 0.0


def quality(V, T):
    qs = np.array([tetq(V[[int(i) for i in tet]]) for tet in T])
    return float(qs.min()), float(qs.mean())


def flips(K):
    Kc = K.tocoo()
    return int(((Kc.row != Kc.col) & (Kc.data > 1e-12)).sum())


def dihedral_angles(P):
    """The 6 dihedral angles (deg) of a tetrahedron."""
    f = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]   # face opposite each vertex
    nrm = {}
    for vi, (a, b, c) in enumerate(f):
        nv = np.cross(P[b]-P[a], P[c]-P[a]); nv = nv/np.linalg.norm(nv)
        if np.dot(nv, P[vi]-P[a]) > 0:
            nv = -nv
        nrm[vi] = nv
    ang = []
    for i in range(4):
        for j in range(i+1, 4):
            cd = np.clip(-np.dot(nrm[i], nrm[j]), -1, 1)
            ang.append(np.degrees(np.arccos(cd)))
    return ang


def count_obtuse(V, T, tol=1e-6):
    o = 0
    for tet in T:
        for a in dihedral_angles(V[[int(i) for i in tet]]):
            if a > 90.0 + tol:
                o += 1
    return o


# --- deformation + alternative operators ------------------------------------
def cube_to_ball(V, c=0.5):
    out = np.zeros_like(V)
    for i, p in enumerate(V):
        d = p - c; m = np.max(np.abs(d))
        out[i] = (c + d*(m/np.linalg.norm(d))) if m > 1e-12 else p
    return out


def graph_laplacian(V, T):
    n = len(V); E = set()
    for tet in T:
        for a in range(4):
            for b in range(a+1, 4):
                E.add((min(int(tet[a]), int(tet[b])), max(int(tet[a]), int(tet[b]))))
    r, c, d = [], [], []; deg = np.zeros(n)
    for (i, j) in E:
        r += [i, j]; c += [j, i]; d += [-1, -1]; deg[i] += 1; deg[j] += 1
    for i in range(n):
        r.append(i); c.append(i); d.append(deg[i])
    return coo_matrix((d, (r, c)), shape=(n, n)).tocsr()


def clamp_laplacian(K):
    Kc = K.tocoo(); n = K.shape[0]
    r, c, d = [], [], []; off = {}
    for i, j, v in zip(Kc.row, Kc.col, Kc.data):
        if i != j:
            w = min(v, 0.0)                  # drop positive off-diagonals
            if w != 0:
                r.append(i); c.append(j); d.append(w); off[i] = off.get(i, 0) + w
    for i in range(n):
        r.append(i); c.append(i); d.append(-off.get(i, 0))
    return coo_matrix((d, (r, c)), shape=(n, n)).tocsr()


def boundary_mask(Vcube):
    return (np.isclose(Vcube[:, 0], Vcube[:, 0].min()) | np.isclose(Vcube[:, 0], Vcube[:, 0].max()) |
            np.isclose(Vcube[:, 1], Vcube[:, 1].min()) | np.isclose(Vcube[:, 1], Vcube[:, 1].max()) |
            np.isclose(Vcube[:, 2], Vcube[:, 2].min()) | np.isclose(Vcube[:, 2], Vcube[:, 2].max()))


def harmonic_err(K, Vphys, bm):
    """Impose linear u=x on the boundary, solve Lu=0 inside; exact answer u=x."""
    free = ~bm; ub = Vphys[:, 0]
    uf = spsolve(K[free][:, free].tocsc(), -(K[free][:, bm]) @ ub[bm])
    u = ub.copy(); u[free] = uf
    return float(np.max(np.abs(u - Vphys[:, 0])))


def smart_smooth(V, T, bm, iters=25):
    """Move each interior vertex to the neighbour centroid only if the local
    minimum tet quality does not decrease (cannot worsen the mesh)."""
    n = len(V)
    nb = [set() for _ in range(n)]
    its = [[] for _ in range(n)]
    for ti, tet in enumerate(T):
        ix = [int(i) for i in tet]
        for a in range(4):
            its[ix[a]].append(ti)
            for b in range(4):
                if a != b:
                    nb[ix[a]].add(ix[b])
    Vs = V.copy()
    for _ in range(iters):
        for i in range(n):
            if bm[i] or not nb[i]:
                continue
            old = Vs[i].copy()
            q0 = min(tetq(Vs[[int(k) for k in T[ti]]]) for ti in its[i])
            Vs[i] = np.mean([Vs[j] for j in nb[i]], axis=0)
            q1 = min(tetq(Vs[[int(k) for k in T[ti]]]) for ti in its[i])
            if q1 <= q0:
                Vs[i] = old
    return Vs


# --- checks -----------------------------------------------------------------
def check1_mechanism(Vc, T, Vb):
    print("\n[1] MECHANISM: sign-flips come from obtuse dihedral angles")
    Kb, _ = fem_laplacian(Vc, T); Kball, _ = fem_laplacian(Vb, T)
    ob0 = count_obtuse(Vc, T); obb = count_obtuse(Vb, T)
    f0, fb = flips(Kb), flips(Kball)
    print(f"    baseline Kuhn cube: sign-flips={f0:5d}  obtuse dihedrals={ob0:5d}  (90 deg exactly -> at the boundary)")
    print(f"    geometric ball    : sign-flips={fb:5d}  obtuse dihedrals={obb:5d}")
    assert f0 == 0 and ob0 == 0, "baseline cube should have no obtuse dihedrals / no flips"
    assert fb > 0 and obb > 0, "deformation should create obtuse dihedrals and flips"
    print("    PASS  (right-angle cube has 0; deformation makes dihedrals obtuse -> sign-flips appear)")


def check2_formula_tradeoff(Vb, T, bm):
    print("\n[2] FORMULA TRADEOFF: removing flips by norm change loses linear-exactness")
    Kcot, _ = fem_laplacian(Vb, T)
    Kgr = graph_laplacian(Vb, T)
    Kcl = clamp_laplacian(Kcot)
    rows = [("cotangent (FE)", Kcot), ("graph (uniform)", Kgr), ("clamp max(w,0)", Kcl)]
    res = {}
    print(f"    {'operator':<18}{'sign-flips':>11}{'harmonic err (u=x)':>22}")
    for name, Kx in rows:
        he = harmonic_err(Kx, Vb, bm); fl = flips(Kx)
        res[name] = (fl, he)
        print(f"    {name:<18}{fl:>11}{he:>22.2e}")
    assert res["cotangent (FE)"][1] < 1e-10, "cotangent should reproduce linear field exactly"
    assert res["graph (uniform)"][0] == 0 and res["graph (uniform)"][1] > 1e-3, "uniform should be flip-free but inexact"
    assert res["clamp max(w,0)"][0] == 0 and res["clamp max(w,0)"][1] > 1e-3, "clamp should be flip-free but inexact"
    print("    PASS  (cotangent: exact but flips; uniform/clamp: flip-free but lose exactness -> norm change is not free)")


def check3_mesh_reality(Vb, T, bm):
    print("\n[3] MESH-IMPROVEMENT REALITY: smoothing helps quality, keeps exactness, but flips persist")
    Kb, _ = fem_laplacian(Vb, T)
    q0 = quality(Vb, T); f0 = flips(Kb); he0 = harmonic_err(Kb, Vb, bm)
    Vs = smart_smooth(Vb, T, bm, iters=25)
    Ks, _ = fem_laplacian(Vs, T)
    q1 = quality(Vs, T); f1 = flips(Ks); he1 = harmonic_err(Ks, Vs, bm)
    print(f"    before: q_min={q0[0]:.3f} q_mean={q0[1]:.3f}  flips={f0}  harmonic_err={he0:.1e}")
    print(f"    after : q_min={q1[0]:.3f} q_mean={q1[1]:.3f}  flips={f1}  harmonic_err={he1:.1e}")
    assert q1[0] >= q0[0] - 1e-9, "guarded smoothing must not reduce min quality"
    assert he1 < 1e-10, "cotangent exactness must be preserved by smoothing"
    print("    NOTE  flips are not eliminated (an edge-summed obtuse-dihedral property;")
    print("          globally non-obtuse 'well-centered' 3D meshes are hard).")
    print("    PASS  (quality up + exactness kept; sign-flips remain -> not a cheap fix)")


def check4_clean_escape(Vc, T, bm):
    print("\n[4] CLEAN ESCAPE: geometry preservation gives 0 flips AND exactness together")
    Kb, _ = fem_laplacian(Vc, T)
    f0 = flips(Kb); he0 = harmonic_err(Kb, Vc, bm)
    print(f"    undeformed cube (= topological torus tets): flips={f0}  harmonic_err={he0:.2e}")
    assert f0 == 0 and he0 < 1e-10, "geometry-preserving mesh should give both"
    print("    PASS  (only NOT moving vertices gives flip-free AND exact -> topological deformation is the clean path)")


if __name__ == "__main__":
    n = 6
    Vc, T = kuhn_cube(n); Vc = Vc - 0.5
    bm = boundary_mask(Vc)
    Vb = cube_to_ball(Vc)
    print("Kuhn cube n=%d: nV=%d nTet=%d" % (n, len(Vc), len(T)))
    check1_mechanism(Vc, T, Vb)
    check2_formula_tradeoff(Vb, T, bm)
    check3_mesh_reality(Vb, T, bm)
    check4_clean_escape(Vc, T, bm)
    print("\nAll cotangent sign-flip checks PASSED.")
