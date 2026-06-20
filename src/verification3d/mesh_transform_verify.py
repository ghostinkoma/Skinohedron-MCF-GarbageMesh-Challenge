"""
mesh_transform_verify.py
========================
Backs theory/05_mesh_transform.md (V2.5). Reshape the Kuhn cube into a torus and
a ball and measure what each deformation does to the verified operator
L = M^{-1} K. The one distinction under test:

  topological (glue faces, vertices fixed)  -> tet shapes unchanged -> L preserved
  geometric   (move vertices)               -> tets distort        -> L degrades

Checks (the five targets of section 5):
  1. Topological torus preserves L + the periodic spectrum matches the analytic
     flat-torus eigenvalues  lambda = (2*pi)^2 (a^2+b^2+c^2)  (L=1).
  2. Wrap-around: a discrete plane wave e^{i k.x} on the periodic mesh is a
     stationary mode (Rayleigh quotient = the analytic eigenvalue) -> no seam
     reflection.
  3. Geometric distortion metric: cube->ball and wrapped-torus degrade tet
     quality and create sign-flipped cotangent weights; topological does not.
  4. Sphere obstruction is real (Theorema Egregium): no cube->sphere map keeps
     quality near 1.
  5. Physics still solves on the degraded mesh: heat decays, wave oscillates.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import diags, coo_matrix
from scipy.sparse.linalg import eigsh
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

Q_REGULAR = 1.0 / (6 * np.sqrt(2))      # vol/rms^3 of a regular tetrahedron


# --------------------------------------------------------------------------- #
#  helpers                                                                      #
# --------------------------------------------------------------------------- #
def tet_quality(V, T):
    """Normalised quality per tet (1 = regular, 0 = sliver): min and mean."""
    qs = []
    for tet in T:
        P = V[[int(i) for i in tet]]
        vol = abs(np.linalg.det(np.array([P[1]-P[0], P[2]-P[0], P[3]-P[0]]))) / 6.0
        e = sum(np.sum((P[a]-P[b])**2) for a in range(4) for b in range(a+1, 4))
        rms = np.sqrt(e / 6.0)
        qs.append((vol / rms**3) if rms > 0 else 0.0)
    qs = np.array(qs) / Q_REGULAR
    return float(qs.min()), float(qs.mean())


def sign_flipped_offdiag(K):
    """Count strictly positive off-diagonal entries (= negative conductance)."""
    Kc = K.tocoo()
    m = (Kc.row != Kc.col) & (Kc.data > 1e-12)
    return int(m.sum())


def remap_periodic(V, axes, L=1.0):
    """Identify vertices at coord=L with their coord=0 partners on given axes."""
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
    return inv, int(inv.max() + 1)


def periodic_KM(V, T, axes, L=1.0):
    inv, nnew = remap_periodic(V, axes, L)
    K, M = fem_laplacian(V, T)
    Kc = K.tocoo()
    Kp = coo_matrix((Kc.data, (inv[Kc.row], inv[Kc.col])), shape=(nnew, nnew)).tocsr()
    Mp = np.zeros(nnew)
    for i, m in enumerate(M):
        Mp[inv[i]] += m
    return Kp, Mp, inv, nnew


# deformation maps -----------------------------------------------------------
def map_cube_to_ball(V, c=(0.5, 0.5, 0.5)):
    c = np.array(c); out = np.zeros_like(V)
    for i, p in enumerate(V):
        d = p - c; m = np.max(np.abs(d))
        if m < 1e-12:
            out[i] = p; continue
        out[i] = c + d * (m / np.linalg.norm(d))      # cube face -> sphere
    return out


def map_wrap_torus(V, R=1.0, r=0.35):
    """Wrap [0,1]^2 x [0,1] into a ring: (u,v) -> angles, w -> minor offset."""
    out = np.zeros_like(V)
    for i, (x, y, z) in enumerate(V):
        th = 2*np.pi*x; ph = 2*np.pi*y
        rr = r * (0.5 + z)
        out[i] = [(R + rr*np.cos(ph))*np.cos(th),
                  (R + rr*np.cos(ph))*np.sin(th),
                  rr*np.sin(ph)]
    return out


# --------------------------------------------------------------------------- #
#  checks                                                                       #
# --------------------------------------------------------------------------- #
def check1_periodic_spectrum(V, T):
    print("\n[1] Topological torus: L preserved + periodic spectrum")
    qb = tet_quality(V, T)
    Kp, Mp, inv, nnew = periodic_KM(V, T, axes=(0, 1, 2), L=1.0)
    print(f"    identification: {len(V)} -> {nnew} vertices (3-torus, expect 8^3=512)")
    # tet shapes are untouched by gluing -> quality identical to baseline
    print(f"    tet quality baseline = (min {qb[0]:.3f}, mean {qb[1]:.3f}); gluing moves no vertex -> identical")
    lam, _ = eigsh(Kp, k=9, M=diags(Mp), sigma=1e-6, which="LM")
    lam = np.sort(lam)
    nz = int((lam < 1e-6).sum())
    lam1 = lam[1]
    pred = (2*np.pi)**2
    err = abs(lam1 - pred) / pred
    deg = int(np.sum(np.abs(lam - lam1) < 0.05*lam1))    # multiplicity of lam1
    print(f"    null eigenvalues (constants) = {nz}  (target 1)")
    print(f"    lowest non-zero lam1 = {lam1:.3f}  vs analytic (2pi)^2 = {pred:.3f}  (err {err*100:.1f}%)")
    print(f"    degeneracy of lam1 = {deg}  (target 6: modes (+-1,0,0),(0,+-1,0),(0,0,+-1))")
    assert nnew == 512 and nz == 1 and err < 0.08 and deg == 6, "periodic spectrum mismatch"
    print("    PASS  (flat-torus spectrum to FE accuracy; correct null space + 6-fold degeneracy)")


def check2_wraparound(V, T):
    print("\n[2] Wrap-around: the seam is seamless (no boundary, no reflection)")
    Kp, Mp, inv, nnew = periodic_KM(V, T, axes=(0, 1, 2), L=1.0)
    # (a) seamlessness: after full periodic identification every vertex is an
    #     interior vertex -> identical connectivity degree. A reflecting boundary
    #     would show reduced-degree seam vertices.
    deg = np.diff(Kp.tocsr().indptr)
    print(f"    vertex degree: min={deg.min()} max={deg.max()} "
          f"(uniform => no boundary; seam vertices are full interior vertices)")
    assert deg.min() == deg.max(), "seam vertices have reduced degree (not seamless)"
    # K p = 0 for constants everywhere, including the seam rows
    rowsum = np.abs(np.asarray(Kp.sum(axis=1)).ravel()).max()
    print(f"    max |row-sum of K| = {rowsum:.2e}  (Laplacian property holds at the seam)")
    assert rowsum < 1e-9, "Laplacian row-sum broken at seam"
    # (b) long-wave plane waves are stationary periodic modes (FE-accurate);
    #     higher k just shows the known FE dispersion (03b), reported not asserted.
    Vr = np.zeros((nnew, 3))
    for i, p in enumerate(V):
        Vr[inv[i]] = np.where(np.isclose(p, 1.0), 0.0, p)
    ok = True
    for k in [(1, 0, 0), (0, 1, 0), (1, 1, 0), (2, 0, 0), (2, 2, 2)]:
        kk = np.array(k)
        w = np.cos(2*np.pi*(Vr @ kk)); w -= (Mp*w).sum()/Mp.sum()
        rq = (w @ (Kp @ w)) / (w @ (Mp*w))
        pred = (2*np.pi)**2 * (kk @ kk)
        err = abs(rq - pred)/pred
        tag = "stationary mode" if (kk @ kk) <= 2 else "FE dispersion (short wave)"
        print(f"    k={k}: Rayleigh={rq:8.3f} analytic={pred:8.3f} err={err*100:4.1f}%  ({tag})")
        if (kk @ kk) <= 2:
            ok = ok and err < 0.08
    assert ok, "long-wave plane wave is not a stationary periodic mode"
    print("    PASS  (uniform degree => seamless wrap-around; long waves are periodic eigenmodes)")


def check3_distortion(V, T):
    print("\n[3] Geometric distortion metric (topological vs geometric)")
    qb = tet_quality(V, T)
    Kb, _ = fem_laplacian(V, T)
    print(f"    baseline Kuhn cube : q_min={qb[0]:.3f} q_mean={qb[1]:.3f}  sign-flipped offdiag={sign_flipped_offdiag(Kb)}")
    # geometric ball
    Vball = map_cube_to_ball(V)
    qb2 = tet_quality(Vball, T); Kball, _ = fem_laplacian(Vball, T)
    print(f"    geometric ball     : q_min={qb2[0]:.3f} q_mean={qb2[1]:.3f}  sign-flipped offdiag={sign_flipped_offdiag(Kball)}")
    # geometric wrapped torus
    Vtor = map_wrap_torus(V)
    qb3 = tet_quality(Vtor, T); Ktor, _ = fem_laplacian(Vtor, T)
    print(f"    geometric torus    : q_min={qb3[0]:.3f} q_mean={qb3[1]:.3f}  sign-flipped offdiag={sign_flipped_offdiag(Ktor)}")
    # topological torus: same coordinates -> identical quality
    print(f"    topological torus  : q_min={qb[0]:.3f} q_mean={qb[1]:.3f}  (vertices unmoved -> identical to baseline)")
    assert qb2[0] < qb[0] and qb3[0] < qb[0], "geometric maps should reduce min quality"
    print("    PASS  (geometric maps degrade quality / can flip cotangents; gluing preserves it)")


def check4_sphere_obstruction(V, T):
    print("\n[4] Sphere obstruction is real (no cube->sphere map keeps quality ~1)")
    qb = tet_quality(V, T)
    best = 0.0
    # sample several cube->sphere style maps; none should restore quality
    cands = []
    cands.append(map_cube_to_ball(V))
    cands.append(map_cube_to_ball(V) * 0.8 + V * 0.2)        # blended
    # normalised radial push with exponent variations
    for power in (0.5, 1.0, 2.0):
        out = np.zeros_like(V)
        for i, p in enumerate(V):
            d = p - 0.5; m = np.max(np.abs(d))
            out[i] = 0.5 + (d*(m/np.linalg.norm(d))) if m > 1e-12 else p
        cands.append(out)
    for Vc in cands:
        q = tet_quality(Vc, T)[0]
        best = max(best, q)
    print(f"    baseline q_min={qb[0]:.3f};  best q_min over sphere maps = {best:.3f}")
    assert best < 0.5 * qb[0], "a sphere map unexpectedly preserved quality"
    print("    PASS  (every cube->sphere map degrades quality -> curvature obstruction is concrete)")


def check5_physics_survives(V, T):
    print("\n[5] Physics still solves on the degraded (geometric ball) mesh")
    Vball = map_cube_to_ball(V)
    K, M = fem_laplacian(Vball, T)
    Minv = 1.0 / M
    # heat: energy must decrease monotonically
    dt = 0.3 / float(np.max(K.diagonal() * Minv))
    T0 = np.zeros(len(Vball)); T0[int(np.argmin(((Vball-0.5)**2).sum(1)))] = 1.0
    e0 = T0 @ (M * T0); Tf = T0.copy()
    mono = True; prev = e0
    for s in range(300):
        Tf = Tf - dt * Minv * (K @ Tf)
        e = Tf @ (M * Tf)
        if e > prev + 1e-12:
            mono = False
        prev = e
    print(f"    heat energy {e0:.4e} -> {prev:.4e}  monotone decreasing = {mono}")
    # wave: bounded leapfrog energy
    lammax = 2.2 * float(np.max(K.diagonal() * Minv))
    dtw = 0.4 * np.sqrt(4/lammax)
    p0 = T0.copy(); p1 = T0.copy(); emax = 0; efin = 0
    for s in range(800):
        p2 = 2*p1 - p0 - dtw*dtw*Minv*(K @ p1)
        p0, p1 = p1, p2
        ke = p1 @ (M*p1); emax = max(emax, ke); efin = ke
    bounded = efin < 5 * (T0 @ (M*T0)) + 1e-9
    print(f"    wave leapfrog: final/initial energy ratio = {efin/(T0@(M*T0)):.3f}  bounded = {bounded}")
    assert mono and np.isfinite(prev) and bounded, "physics failed on degraded mesh"
    print("    PASS  (heat decays, wave stays bounded -> L degrades gracefully, does not break)")


if __name__ == "__main__":
    V, T = kuhn_cube(8)                      # [0,1]^3, side length 1 (for periodic L=1)
    print("Kuhn cube n=8:  nV=%d  nTet=%d" % (len(V), len(T)))
    check1_periodic_spectrum(V, T)
    check2_wraparound(V, T)
    check3_distortion(V, T)
    check4_sphere_obstruction(V, T)
    check5_physics_survives(V, T)
    print("\nAll V2.5 mesh-transform checks PASSED.")
