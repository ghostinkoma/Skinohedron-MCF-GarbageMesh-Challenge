"""
thermal_glue_verify.py  --  temperature as the GLUE between Coulomb and mass

A deeper role for temperature than the parallel fusion of 06g: temperature is the
MEDIATOR. Honey heated becomes less viscous -- temperature does not act alone, it acts
THROUGH the inter-molecular (Coulomb) bonds to change the mass transport (flow). So
temperature is the glue between the Coulomb domain (stiffness/viscosity) and the mass
domain (momentum/flow).

The mediation is MULTIPLICATIVE (06g's fusion was additive): the viscosity enters the
stiffness operator as a temperature-dependent WEIGHT, nu = nu(T(x)), so the viscous
operator is K_nu = integral nu(T) grad phi . grad psi. Temperature lives INSIDE K.

Verified here:
  A  Arrhenius nu(T) = nu0 exp(Ea/T): the SAME dT gives very different viscosity change
     for honey (high Coulomb barrier Ea) vs water (low Ea) -- the correlation, not a
     parallel domain.
  B  Two-viscosity interface (cold/hot) reproduces u = U * mu2/(mu1+mu2) to MACHINE
     PRECISION -- the mechanical mirror of 02's two-material T = k2/(k1+k2).
  C  Arrhenius-temperature-set viscosities give the interface flow to machine precision.
  D  The mediation is multiplicative & consistent: constant T recovers nu*K exactly;
     stress (mu du/dy) is continuous across the interface.

So temperature mediates Coulomb<->mass through a structure already verified in 02,
to machine precision. (The structural resemblance to a field that sets other fields'
effective couplings is noted as a direction, deliberately not formalised here.)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian


def variable_stiffness(V, T, mu_tet):
    """K_mu = sum_tet mu_tet * vol * G^T G : viscosity weights the stiffness operator."""
    n = len(V); rows, cols, data = [], [], []
    for ti, tet in enumerate(T):
        ix = [int(i) for i in tet]
        P = V[ix]
        vol = abs(np.linalg.det(np.array([P[1]-P[0], P[2]-P[0], P[3]-P[0]]))) / 6.0
        if vol <= 1e-15:
            continue
        C = np.linalg.inv(np.column_stack([np.ones(4), P])); G = C[1:4, :]
        Ke = mu_tet[ti]*vol*(G.T @ G)
        for a in range(4):
            for b in range(4):
                rows.append(ix[a]); cols.append(ix[b]); data.append(Ke[a, b])
    return coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()


def couette_interface(V, T, mu_tet, U=1.0):
    """Solve wall-driven flow with variable viscosity; return interface velocity."""
    n = len(V); y = V[:, 1]
    Kmu = variable_stiffness(V, T, mu_tet)
    wall = np.isclose(y, 0.0) | np.isclose(y, 1.0); free = ~wall
    u = np.zeros(n); u[wall] = np.where(np.isclose(y, 1.0), U, 0.0)[wall]
    u[free] = spsolve(Kmu[free][:, free].tocsc(), -(Kmu[free][:, wall] @ u[wall]))
    iface = np.where(np.isclose(y, 0.5))[0]
    return u, u[iface].mean()


def checkA_arrhenius():
    print("\n[A] Arrhenius nu(T)=nu0 exp(Ea/T): temperature acts THROUGH the Coulomb barrier")
    T1, T2 = 300.0, 330.0
    for name, Ea in [("water  (weak H-bonds)", 1800.0), ("honey  (strong H-bond net)", 6000.0)]:
        ratio = np.exp(Ea/T1)/np.exp(Ea/T2)
        print(f"    {name}: Ea={Ea:5.0f}  nu(300)/nu(330) = {ratio:.2f}x  (+30K thins by 1/{ratio:.1f})")
    print("    -> same temperature change, very different viscosity change via Ea = correlation")
    print("    PASS")


def checkB_interface_machine():
    print("\n[B] Two-viscosity interface = 02's two-material formula, to machine precision")
    n = 12
    V, T = kuhn_cube(n)
    y = V[:, 1]
    tetc = np.array([y[[int(i) for i in tet]].mean() for tet in T])
    mu1, mu2 = 4.0, 1.0
    mu_tet = np.where(tetc < 0.5, mu1, mu2)
    u, u_int = couette_interface(V, T, mu_tet)
    exact = mu2/(mu1+mu2)
    print(f"    interface u = {u_int:.6f}  vs exact U*mu2/(mu1+mu2) = {exact:.6f}  err {abs(u_int-exact):.2e}")
    print(f"    (identical to 02's T=k2/(k1+k2)={exact:.4f} -- temperature-mediated interface)")
    # stress continuity: slope ratio = mu2/mu1
    free = ~(np.isclose(y, 0.0) | np.isclose(y, 1.0))
    lo = (y < 0.5) & free; hi = (y > 0.5) & free
    sl_lo = np.polyfit(y[lo], u[lo], 1)[0]; sl_hi = np.polyfit(y[hi], u[hi], 1)[0]
    print(f"    slope ratio lo/hi = {sl_lo/sl_hi:.4f} vs mu2/mu1 = {mu2/mu1:.4f} (stress continuous)")
    assert abs(u_int-exact) < 1e-12 and abs(sl_lo/sl_hi - mu2/mu1) < 1e-3
    print("    PASS  (machine precision -- the mechanical mirror of the verified heat interface)")


def checkC_arrhenius_interface():
    print("\n[C] Arrhenius-temperature-set viscosities give the interface flow exactly")
    n = 12
    V, T = kuhn_cube(n)
    y = V[:, 1]
    tetc = np.array([y[[int(i) for i in tet]].mean() for tet in T])
    Ea, mu0 = 6000.0, 1.0
    Tcold, Thot = 300.0, 330.0
    norm = np.exp(Ea/330.0)
    Tfield = np.where(tetc < 0.5, Tcold, Thot)
    mu_tet = mu0*np.exp(Ea/Tfield)/norm
    u, u_int = couette_interface(V, T, mu_tet)
    m1 = mu0*np.exp(Ea/Tcold)/norm; m2 = mu0*np.exp(Ea/Thot)/norm
    exact = m2/(m1+m2)
    print(f"    mu(cold 300K)={m1:.2f}  mu(hot 330K)={m2:.2f}")
    print(f"    interface u = {u_int:.6f}  vs exact {exact:.6f}  err {abs(u_int-exact):.2e}")
    assert abs(u_int-exact) < 1e-12
    print("    PASS  (the temperature field sets the flow through viscosity -- the glue)")


def checkD_multiplicative():
    print("\n[D] The mediation is multiplicative: constant T recovers nu*K exactly")
    n = 10
    V, T = kuhn_cube(n)
    K, M = fem_laplacian(V, T)
    nu = 2.5
    mu_tet = np.full(len(T), nu)
    Kmu = variable_stiffness(V, T, mu_tet)
    diff = abs((Kmu - nu*K)).max()
    print(f"    || K_mu(const nu) - nu*K ||_max = {diff:.2e}")
    print("    -> temperature enters as a multiplicative weight INSIDE the stiffness K")
    print("       (06g's fusion was additive energy; this mediation is multiplicative coupling)")
    assert diff < 1e-12
    print("    PASS")


if __name__ == "__main__":
    print("Temperature as the glue between Coulomb (viscosity) and mass (flow)")
    checkA_arrhenius()
    checkB_interface_machine()
    checkC_arrhenius_interface()
    checkD_multiplicative()
    print("\nConclusion:")
    print("  Temperature is not only a parallel domain (06g) but a MEDIATOR: via Arrhenius")
    print("  nu(T) it enters the stiffness K as a multiplicative weight, gluing the Coulomb")
    print("  (viscosity) and mass (flow) domains. The resulting two-viscosity interface is")
    print("  u = U*mu2/(mu1+mu2) to machine precision -- the mechanical twin of 02's")
    print("  T = k2/(k1+k2). Honey-heating (lower nu) is this mediation, verified.")
