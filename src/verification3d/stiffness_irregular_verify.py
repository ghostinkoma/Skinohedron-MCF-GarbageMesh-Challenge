"""
stiffness_irregular_verify.py  --  convergence of the stiffness K on irregular, non-periodic meshes

Most exact results elsewhere in the notebook use fields that lie in the P1 space or
periodic sinusoids on the uniform Kuhn mesh (implementation checks, see theory/README.md,
"How to read the claims"). This script tests the stiffness K of fem_laplacian on a
problem that is NOT representable exactly, on a NON-periodic domain, on uniform AND
irregular meshes, and measures convergence ORDERS.

Problem: -Laplace(u) = f on the unit cube, Dirichlet data from the exact solution
    u(x, y, z) = sin(2x + 1) * exp(y) * cos(z),   f = 4u
(smooth, not polynomial, not periodic). Load vector b = M f (lumped). Errors use a
degree-2 four-point quadrature on every tetrahedron:
    L2 error  ||u - u_h||_L2 / ||u||_L2
    H1 error  ||grad(u - u_h)||_L2 / ||grad u||_L2
and the lowest Dirichlet eigenvalue of (K, M) against 3 pi^2.

Meshes: the Kuhn cube (uniform), and the same cube with every interior vertex moved by
a random offset of up to 20% of h in each coordinate (fixed seed; the perturbation is
the same fraction of h at every level, so the mesh does not become regular under
refinement). No tetrahedron is inverted (checked). The irregular meshes have thousands
of sign-flipped cotangent weights (positive off-diagonal K entries), i.e. the discrete
maximum principle is lost.

  CHECK A  uniform mesh:   L2 ~ h^2, H1 ~ h^1, eigenvalue ~ h^2
  CHECK B  irregular mesh: the same orders, despite the sign-flipped weights

Known theory (theory/README.md): H1 ~ h from Cea's lemma + interpolation, L2 ~ h^2 from
Aubin-Nitsche, eigenvalue ~ h^2 from Babuska-Osborn. These need shape regularity, not the
maximum principle, which is why B converges. (Contrast: the POINTWISE consistency of the
operator on a persistently irregular surface mesh does not converge, paper section 6;
that is a different, strong-form quantity.)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve, eigsh
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

NS = (8, 12, 16, 24)
JITTER = 0.2

# degree-2 four-point rule on a tetrahedron (barycentric coordinates, equal weights)
_A, _B = 0.5854101966249685, 0.1381966011250105
BARY = np.array([[_A, _B, _B, _B], [_B, _A, _B, _B], [_B, _B, _A, _B], [_B, _B, _B, _A]])


def u_exact(X):
    return np.sin(2 * X[..., 0] + 1) * np.exp(X[..., 1]) * np.cos(X[..., 2])


def grad_exact(X):
    s, c = np.sin(2 * X[..., 0] + 1), np.cos(2 * X[..., 0] + 1)
    e, cz, sz = np.exp(X[..., 1]), np.cos(X[..., 2]), np.sin(X[..., 2])
    return np.stack([2 * c * e * cz, s * e * cz, -s * e * sz], axis=-1)


def f_rhs(X):
    return 4.0 * u_exact(X)          # -Laplace(u) = (4 - 1 + 1) u


def build_mesh(n, jitter, seed=1):
    V, T = kuhn_cube(n)
    V = V.astype(float)
    sign0 = np.sign(np.linalg.det(V[T][:, 1:] - V[T][:, :1]))
    if jitter > 0:
        rng = np.random.default_rng(seed)
        interior = ~((np.isclose(V, 0.0) | np.isclose(V, 1.0)).any(axis=1))
        V[interior] += jitter / n * rng.uniform(-1, 1, size=(interior.sum(), 3))
    sign = np.sign(np.linalg.det(V[T][:, 1:] - V[T][:, :1]))
    assert np.all(sign == sign0), "a tetrahedron was inverted by the perturbation"
    return V, T


def solve(n, jitter):
    V, T = build_mesh(n, jitter)
    K, M = fem_laplacian(V, T)
    bnd = (np.isclose(V, 0.0) | np.isclose(V, 1.0)).any(axis=1)
    free = ~bnd
    u = np.zeros(len(V)); u[bnd] = u_exact(V[bnd])
    rhs = M * f_rhs(V) - K[:, bnd] @ u[bnd]
    u[free] = spsolve(K[free][:, free].tocsc(), rhs[free])

    P = V[T]
    vol = np.abs(np.linalg.det(P[:, 1:] - P[:, :1])) / 6.0
    Q = np.einsum("qk,tkd->tqd", BARY, P)                 # quadrature points
    w = vol[:, None] / 4.0
    uh_q = np.einsum("qk,tk->tq", BARY, u[T])
    l2 = np.sqrt(np.sum(w * (u_exact(Q) - uh_q) ** 2) / np.sum(w * u_exact(Q) ** 2))
    C = np.linalg.inv(np.concatenate([np.ones((len(T), 4, 1)), P], axis=2))
    grad_h = np.einsum("tdk,tk->td", C[:, 1:, :], u[T])  # constant P1 gradient per tet
    gq = grad_exact(Q)
    h1 = np.sqrt(np.sum(w * np.sum((gq - grad_h[:, None, :]) ** 2, axis=-1))
                 / np.sum(w * np.sum(gq ** 2, axis=-1)))

    lam = eigsh(K[free][:, free].tocsc(), k=1, M=diags(M[free]), sigma=0.0,
                which="LM", return_eigenvectors=False)[0]
    eig = abs(lam - 3 * np.pi ** 2) / (3 * np.pi ** 2)

    Kc = K.tocoo()
    flips = int(np.sum((Kc.row < Kc.col) & (Kc.data > 1e-14)))
    return l2, h1, eig, flips


def run_family(tag, jitter):
    rows = [solve(n, jitter) for n in NS]
    h = 1.0 / np.array(NS, float)
    order = lambda k: float(np.polyfit(np.log(h), np.log([r[k] for r in rows]), 1)[0])
    for n, (l2, h1, eig, flips) in zip(NS, rows):
        print(f"    n={n:>2}: L2 {l2:.3e}   H1 {h1:.3e}   eig {eig:.3e}   sign-flipped edges {flips}")
    o = (order(0), order(1), order(2))
    print(f"    orders: L2 ~ h^{o[0]:.2f}   H1 ~ h^{o[1]:.2f}   eigenvalue ~ h^{o[2]:.2f}")
    return o, rows


def check(o):
    l2, h1, eig = o
    assert 1.8 < l2 < 2.2, f"L2 order {l2:.2f} not ~2"
    assert 0.9 < h1 < 1.1, f"H1 order {h1:.2f} not ~1"
    assert 1.8 < eig < 2.2, f"eigenvalue order {eig:.2f} not ~2"


def main():
    print("Stiffness K (fem_laplacian): -Laplace u = f, u = sin(2x+1) e^y cos z, Dirichlet, unit cube")

    print("\n[A] uniform Kuhn mesh")
    oA, _ = run_family("uniform", 0.0)
    check(oA)
    print("    PASS  -- L2 O(h^2), H1 O(h), eigenvalue O(h^2) (Aubin-Nitsche / Cea / Babuska-Osborn)")

    print(f"\n[B] irregular mesh (interior vertices moved up to {JITTER:.0%} of h, no inverted tets)")
    oB, rowsB = run_family("irregular", JITTER)
    assert all(r[3] > 0 for r in rowsB), "expected sign-flipped cotangent weights"
    check(oB)
    print("    PASS  -- the same optimal orders, although the maximum principle is lost")

    print("\nConclusion: on a non-periodic problem whose solution is not in P1, the stiffness")
    print("K converges at the optimal Galerkin rates on uniform AND irregular meshes; the")
    print("thousands of sign-flipped weights cost the maximum principle, not convergence.")


if __name__ == "__main__":
    main()
