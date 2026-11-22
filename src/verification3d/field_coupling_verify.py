"""
field_coupling_verify.py  --  gravity as a field on the same operator, and field-driven NS

A deliberately CONSERVATIVE verification of operator structure only. Extends 06e (Coulomb):
gravity is the same Poisson operator, and the field enters NS as a body force. What is
verified:

  CHECK A  Gravity rides the SAME Poisson operator as Coulomb. Newtonian gravity
    K Phi = -4 pi G * b(rho) and electrostatics K phi = b(rho_q) / eps are the same
    operator K (K ~ -Laplacian) with different constants. The self-gravity potential of a
    uniform ball (R=1 m, water) is solved with the project's operator itself
    (fem_laplacian on the Kuhn cube, box [-2,2]^3, exact exterior potential on the box
    boundary). The hydrostatic central pressure rho*(Phi(R) - Phi(0)) converges to the
    analytic (2/3) pi G rho^2 R^2 under refinement (about first order, limited by the
    staircase approximation of the ball's surface on the cube mesh).

  CHECK B  Reference magnitudes ONLY (closed-form, NOT an NS solve). The same water
    droplet (R=1 m) has very different characteristic pressure under three fields:
    zero-g (Young-Laplace 2 sigma/R), uniform earth-g (rho g 2R), self-gravity (the FE
    value from CHECK A). This illustrates that the field sets the body-force scale; it
    does not solve Navier-Stokes with a free surface.

  CHECK C  Stitching = linear superposition, tested for real. A gravity source and a
    Coulomb source are solved separately on the same K (Kuhn cube, fem_laplacian), and
    then the COMBINED source is solved once. The combined potential equals the sum of the
    separate potentials to machine precision; since the body force -rho grad(phi) is
    linear in phi, the forces superpose as well.

  CHECK D  The operator is dimension-independent. The same P1 cotangent/gradient operator
    is exact for a linear field in 2-D (triangle) and 3-D (tetrahedron): u^T K u = 1 for
    u = x on the unit cell. The n-simplex is the minimal cell fixing an affine field by its
    n+1 vertices -- the natural carrier for one scalar field per vertex in any dimension.

Scope: this verifies operator STRUCTURE only -- that gravity shares the Poisson operator
with Coulomb (A), and that potentials from several sources superpose on it (C). B is an
illustration with closed-form magnitudes; no NS solve is performed in this script.

(History: an earlier version solved CHECK A with a 1-D radial finite-difference
scheme that did not use the project's K, presented B as an NS result, and computed
CHECK C as norm(x - x), which is identically zero. Corrected here; see issue #3, A2.)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse.linalg import spsolve
from cavity_pspg_verify import mesh_square, assemble
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian

G = 6.674e-11
RHO = 1000.0
R = 1.0
SIGMA = 0.072


def _box_mesh(n, L):
    """Kuhn cube mapped to the box [-L, L]^3 (n must be a multiple of 4 so that
    the origin and the point (R,0,0) with R = L/2 are mesh vertices)."""
    V, T = kuhn_cube(n)
    return (V - 0.5) * 2.0 * L, T


def _load_vector(V, T, dens_tet):
    """Consistent P1 load vector b_i = sum_t vol_t/4 * dens_t (piecewise-constant
    density per tetrahedron)."""
    P = V[T]
    vol = np.abs(np.linalg.det(P[:, 1:] - P[:, :1])) / 6.0
    b = np.zeros(len(V))
    np.add.at(b, T.ravel(), np.repeat(vol * dens_tet / 4.0, 4))
    return b


def _ball_potential_exact(r):
    M_tot = RHO * (4.0 / 3.0) * np.pi * R**3
    inside = -2.0 * np.pi * G * RHO * (R**2 - r**2 / 3.0)
    outside = -G * M_tot / np.maximum(r, 1e-300)
    return np.where(r < R, inside, outside)


def _solve_self_gravity(n, L=2.0 * R):
    V, T = _box_mesh(n, L)
    K, _ = fem_laplacian(V, T)
    cen = V[T].mean(axis=1)
    rho_t = np.where(np.linalg.norm(cen, axis=1) < R, RHO, 0.0)
    b = _load_vector(V, T, rho_t)
    r = np.linalg.norm(V, axis=1)
    bnd = np.isclose(np.abs(V), L).any(axis=1); free = ~bnd
    phi = np.zeros(len(V)); phi[bnd] = _ball_potential_exact(r[bnd])
    # K ~ -Laplacian, and Laplacian(Phi) = 4 pi G rho  =>  K Phi = -4 pi G b
    rhs = -4.0 * np.pi * G * b - K[:, bnd] @ phi[bnd]
    phi[free] = spsolve(K[free][:, free].tocsc(), rhs[free])
    i0 = int(np.argmin(r))
    iR = int(np.argmin(np.linalg.norm(V - np.array([R, 0.0, 0.0]), axis=1)))
    P_c = RHO * (phi[iR] - phi[i0])          # dp = -rho dPhi, p(R) = 0
    phi_err = np.max(np.abs(phi - _ball_potential_exact(r))) / np.max(np.abs(_ball_potential_exact(r)))
    return P_c, phi_err


def checkA_gravity_is_poisson():
    print("\n[A] Gravity rides the same Poisson operator K as Coulomb (project operator)")
    print("    fem_laplacian on the Kuhn cube, box [-2,2]^3, uniform water ball R=1 m")
    P_exact = (2.0 / 3.0) * np.pi * G * RHO**2 * R**2
    errs = []; P_c = None
    for n in (8, 16, 24, 32):
        P_c, phi_err = _solve_self_gravity(n)
        e = abs(P_c - P_exact) / P_exact
        errs.append(e)
        print(f"    n={n:>2}: central pressure {P_c:.4e} Pa vs analytic {P_exact:.4e}"
              f"  err {e*100:5.2f}%   |Phi-Phi_exact|_max/|Phi|_max = {phi_err:.2e}")
    print("    gravity  : K Phi = -4 pi G b(rho),  body force -rho grad(Phi)")
    print("    Coulomb  : K phi = b(rho_q)/eps  ,  body force -rho_q grad(phi)")
    assert all(e2 < e1 for e1, e2 in zip(errs, errs[1:])), "no convergence"
    assert errs[-1] < 0.01, "central pressure not within 1% at the finest mesh"
    print("    PASS  -- same Poisson operator K, different constant; converges under refinement")
    return P_c                                 # finest-mesh FE value


def checkB_field_changes_ns(p_self_fe):
    print("\n[B] Reference magnitudes only (closed form; NOT an NS solve)")
    p_st = 2 * SIGMA / R
    p_hydro = RHO * 9.8 * (2 * R)
    print(f"    zero-g (Young-Laplace 2 sigma/R):  {p_st:.4f} Pa (isotropic)")
    print(f"    earth-g (rho g 2R):                {p_hydro:.1f} Pa (top-bottom asymmetric)")
    print(f"    self-gravity (FE value from [A]):  {p_self_fe:.4e} Pa (spherically symmetric)")
    print("    -> the field sets the body-force scale; differing by ~8 orders of magnitude.")
    print("       ILLUSTRATION: no Navier-Stokes / free-surface solve is performed here.")


def checkC_stitching_superposition():
    print("\n[C] Stitching = linear superposition of potentials on one operator (Kuhn cube)")
    n = 12; V, T = kuhn_cube(n); K, M = fem_laplacian(V, T)
    x, y, z = V[:, 0], V[:, 1], V[:, 2]
    bnd = (np.isclose(V, 0.0) | np.isclose(V, 1.0)).any(axis=1); free = ~bnd
    Kf = K[free][:, free].tocsc()
    src_g = -4.0 * np.pi * (M * np.ones(len(V)))                  # gravity-like source
    rho_q = np.exp(-((x - 0.5)**2 + (y - 0.5)**2 + (z - 0.5)**2) / (2 * 0.1**2))
    rho_q -= np.sum(M * rho_q) / np.sum(M)
    src_c = M * rho_q                                             # Coulomb-like source
    phi_g = np.zeros(len(V)); phi_g[free] = spsolve(Kf, src_g[free])
    phi_c = np.zeros(len(V)); phi_c[free] = spsolve(Kf, src_c[free])
    phi_both = np.zeros(len(V)); phi_both[free] = spsolve(Kf, (src_g + src_c)[free])
    err = np.linalg.norm(phi_both - (phi_g + phi_c)) / np.linalg.norm(phi_both)
    print(f"    ||phi(src_g + src_c) - (phi(src_g) + phi(src_c))|| / ||phi|| = {err:.2e}")
    assert err < 1e-12
    print("    PASS  -- potentials from several sources superpose on one operator; the body")
    print("            force -rho grad(phi) is linear in phi, so the forces superpose too")


def checkD_dimension_independent():
    print("\n[D] The operator is dimension-independent (n-simplex = minimal affine cell)")
    n = 32; V, T = mesh_square(n); K, Ml, Bx, By, tris = assemble(V, T)
    x = V[:, 0]; KE2 = x@(K@x)
    V3, T3 = kuhn_cube(8); K3, M3 = fem_laplacian(V3, T3); x3 = V3[:, 0]
    KE3 = x3@(K3@x3)
    print(f"    2-D triangle: u^T K u for u=x = {KE2:.4f}  (analytic 1.0)")
    print(f"    3-D tetrahedron: u^T K u for u=x = {KE3:.4f}  (analytic 1.0)")
    assert abs(KE2-1) < 1e-9 and abs(KE3-1) < 1e-9
    print("    PASS  -- same P1 gradient operator in 2-D and 3-D; the n-simplex's n+1 vertices")
    print("            fix an affine field uniquely -> natural carrier per dimension")


if __name__ == "__main__":
    print("Verifiable core of 'different fields on one tetrahedron' (operator structure only)")
    p_self = checkA_gravity_is_poisson()
    checkB_field_changes_ns(p_self)
    checkC_stitching_superposition()
    checkD_dimension_independent()
    print("\nConclusion (what is rigorously established):")
    print("  - gravity = Coulomb = the same Poisson operator K (different constant),")
    print("    solved with the project's operator and converging to the analytic value;")
    print("  - potentials from several sources superpose linearly on one operator;")
    print("  - the operator is dimension-independent (the n-simplex is the minimal affine cell).")
    print("  (operator structure only; [B] is a closed-form illustration, not an NS solve.)")
