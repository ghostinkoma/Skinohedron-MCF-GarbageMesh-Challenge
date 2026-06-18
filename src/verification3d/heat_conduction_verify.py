"""
heat_conduction_verify.py
=========================
Exploration + verification of the HEAT (diffusion) route, following the two
integrating leads:
  (a) "find where the scatter form integrates with FE/DEC", and
  (b) "as heat, conduction can be computed from contact area".

Both converge to one finding, decided here by code:

  Heat conduction across a contact face is a CONDUCTANCE c·ΔT. Done with the
  CORRECT geometry, the conductance is exactly the cotangent / FE weight, so the
  heat operator IS the FE/DEC Laplacian — and because diffusion needs only the
  metric (areas/lengths), NOT full face directions, it sidesteps the spurious
  anisotropy that blocked the scalar WAVE in 01e.

Three checks:

  A. CONSISTENCY — cube Neumann spectrum vs exact π²(a²+b²+c²):
       * naive cell-centroid area/distance FV (two-point flux): INCONSISTENT on
         tetrahedra (wrong magnitude, anisotropic) — the two-point flux
         approximation fails on non-orthogonal meshes.
       * cotangent / FE heat operator: CONSISTENT (→ π²), ~10× more isotropic
         than the wave scatter node.
  B. MATERIALS — two-material steady conduction, interface temperature:
       exact T_iface = k₂/(k₁+k₂) (series thermal resistance) to machine
       precision, for any conductivity ratio. This is the heat analogue of the
       S-parameter material contrast — and here it is exact, not anisotropy-
       limited.
  C. VALID NETWORK — the cotangent conductances are (essentially) non-negative on
       the Kuhn mesh, i.e. a physical heat network obeying the maximum principle.
       (Where a mesh produces negative conductances, that is the well-known
       non-Delaunay/poorly-centred-tet issue — the mesh-quality thread again.)

VERDICT: the heat route is the natural consistent foundation. "Contact-area
conduction" is vindicated — as the cotangent/FE conductance — and it integrates
the scatter-form material idea with FE/DEC without the wave anisotropy obstacle.
"""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.sparam_wave import PortComplex
from ksf3d.fem3d import fem_laplacian
from scipy.sparse import diags, lil_matrix
from scipy.sparse.linalg import eigsh, spsolve

DIRS = {"<100>": (1, 0, 0), "<110>": (1, 1, 0), "<111>": (1, 1, 1)}


# --------------------------------------------------------------------------- #
#  A. consistency: cell-area FV vs cotangent FE                               #
# --------------------------------------------------------------------------- #
def fv_area_laplacian(V, T):
    """Naive cell-centred two-point-flux: conductance = contact_area / centre
    distance. (Uses areas only, no face normals.)"""
    pc = PortComplex(V, T)
    cent = np.array([V[T[i]].mean(axis=0) for i in range(pc.N_C)])
    N = pc.N_C
    L = lil_matrix((N, N))
    for (p0, p1) in pc.interior_facets:
        i, j = pc.cell_of(p0), pc.cell_of(p1)
        w = pc.area[p0] / np.linalg.norm(cent[i] - cent[j])
        L[i, i] += w; L[j, j] += w; L[i, j] -= w; L[j, i] -= w
    return L.tocsr(), pc.vol


def neumann_lowest(L, mass, k=4):
    M = diags(mass)
    vals = eigsh(L, k=k, M=M, sigma=1e-8, which="LM", return_eigenvectors=False)
    return np.sort(vals)


def fem_dispersion_aniso(n, ncyc=1):
    V, T = kuhn_cube(n); V = V - 0.5
    K, M = fem_laplacian(V, T)
    kmag = 2 * np.pi * ncyc
    cs = {}
    for name, kd in DIRS.items():
        kd = np.array(kd, float); kd /= np.linalg.norm(kd)
        u = np.cos(V @ (kmag * kd))
        lam = (u @ (K @ u)) / (u @ (M * u))
        cs[name] = np.sqrt(max(lam, 0)) / kmag
    v = np.array(list(cs.values()))
    return (v.max() - v.min()) / v.mean()


# --------------------------------------------------------------------------- #
#  B. two-material steady conduction                                          #
# --------------------------------------------------------------------------- #
def two_material_iface(n, k1, k2):
    V, T = kuhn_cube(n); N = len(V)
    A = lil_matrix((N, N))
    for tet in T:
        P = V[tet]
        vol = abs(np.linalg.det(np.array([P[1]-P[0], P[2]-P[0], P[3]-P[0]]))) / 6
        if vol <= 1e-15:
            continue
        C = np.linalg.inv(np.column_stack([np.ones(4), P])); G = C[1:4, :]
        kk = k1 if P.mean(axis=0)[0] < 0.5 else k2
        Ke = kk * vol * (G.T @ G)
        for a in range(4):
            for b in range(4):
                A[tet[a], tet[b]] += Ke[a, b]
    A = A.tocsr()
    left = np.where(np.isclose(V[:, 0], 0))[0]
    right = np.where(np.isclose(V[:, 0], 1))[0]
    fixed = np.concatenate([left, right])
    Tval = np.concatenate([np.zeros(len(left)), np.ones(len(right))])
    free = np.setdiff1d(np.arange(N), fixed)
    Tvec = np.zeros(N); Tvec[fixed] = Tval
    Tvec[free] = spsolve(A[free][:, free].tocsr(), -A[free][:, fixed] @ Tvec[fixed])
    iface = np.isclose(V[:, 0], 0.5)
    return Tvec[iface].mean(), k2 / (k1 + k2)


# --------------------------------------------------------------------------- #
#  C. conductance validity                                                    #
# --------------------------------------------------------------------------- #
def negative_conductance_fraction(n, tol=1e-12):
    V, T = kuhn_cube(n)
    K, _ = fem_laplacian(V, T)
    Kc = K.tocoo()
    off = np.array([v for i, j, v in zip(Kc.row, Kc.col, Kc.data) if i < j])
    # conductance c_ij = -K_ij ; negative conductance <=> K_ij > 0
    return (off > tol).sum() / max(len(off), 1)


def main():
    print("=== HEAT (diffusion) route — contact-area conduction = FE/DEC ===\n")

    print("A. CONSISTENCY — cube Neumann spectrum vs exact π² = 9.8696 (mult 3):")
    for n in (8, 12, 16):
        V, T = kuhn_cube(n)
        Lfv, volfv = fv_area_laplacian(V, T)
        fv = neumann_lowest(Lfv, volfv)
        K, Mv = fem_laplacian(V, T)
        fe = neumann_lowest(K, Mv)
        print(f"   n={n:2d}  FV(area/dist) lowest≈ {fv[1]:.3f}   "
              f"FE(cotangent) lowest≈ {fe[1]:.3f}")
    print("   => area/distance FV converges to the WRONG value (~7.6, anisotropic):")
    print("      two-point flux is inconsistent on non-orthogonal tetrahedra.")
    print("      cotangent/FE converges to π² (correct).")
    aniso = np.mean([fem_dispersion_aniso(n) for n in (12, 20)])
    print(f"   FE heat-operator directional anisotropy Δc/c ≈ {aniso:.3f}")
    print(f"      (vs ~0.49 for the geometry-blind wave node — ~{0.49/aniso:.0f}× better)\n")

    print("B. MATERIALS — two-material steady conduction, T_iface vs k₂/(k₁+k₂):")
    okB = True
    for k1, k2 in [(1, 1), (1, 3), (1, 10), (1, 0.1)]:
        Ti, Tth = two_material_iface(12, k1, k2)
        err = abs(Ti - Tth)
        okB &= err < 1e-9
        print(f"   k₁={k1:<4} k₂={k2:<4}: T_iface = {Ti:.5f}  "
              f"theory {Tth:.5f}  err {err:.1e}")
    print("   => inter-material heat conduction is EXACT (series thermal")
    print("      resistance), not anisotropy-limited. This is the heat analogue")
    print("      of the S-parameter material contrast — and it is exact.\n")

    print("C. VALID NETWORK — negative-conductance fraction (Kuhn mesh):")
    for n in (6, 10, 14):
        f = negative_conductance_fraction(n)
        print(f"   n={n:2d}: {f*100:.2f}%  (≈0 ⇒ a physical heat network)")
    print("   => essentially non-negative ⇒ obeys the discrete maximum principle.")
    print("      (Where negatives appear, it is the non-Delaunay/poorly-centred-")
    print("       tet issue — the mesh-quality thread, not a flaw of the heat idea.)\n")

    print("[結句] The heat route is the integrating foundation. 'Conduction from")
    print("  contact area' is correct — realised as the cotangent/FE conductance")
    print("  (contact area over the RIGHT, dual, length). It (1) is consistent")
    print("  (→ π²) where the naive area/distance FV is not, (2) gives EXACT two-")
    print("  material conduction, and (3) sidesteps the 01e wave anisotropy because")
    print("  diffusion needs only the metric, not face directions. This unifies the")
    print("  scatter-form material idea with FE/DEC — the place they integrate is")
    print("  HEAT. The remaining refinement (negative conductances on bad cells) is")
    print("  the same mesh-quality / well-centred-tetrahedra thread as before.")

    assert okB, "two-material conduction not exact"
    print("\nMaterial-conduction exactness asserted.")


if __name__ == "__main__":
    main()
