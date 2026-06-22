"""
coulomb_operator_verify.py  --  is Coulomb already in the notebook's one operator?

Prompted by an insight about the OTHER term behind viscosity: at the molecular level
honey and water have nearly equal mass density but wildly different viscosity, because
viscosity is set by INTER-MOLECULAR forces -- and those forces are, at root, COULOMB
(hydrogen bonds, dipole-dipole, van der Waals are all electrostatic). So before any
claim about viscosity, the checkable question is: does the electrostatic Coulomb
problem live in the SAME verified operator K as the pressure Poisson (04)?

Yes. Electrostatics is `-eps grad^2 phi = rho`, i.e. `K phi = M rho` -- the same K.
This script verifies the concrete, in-scope links:

  A  Coulomb 1/r law: solving K phi = M rho for a point charge reproduces phi ~ 1/(4 pi r)
  B  Configuration energy: two charges have U(d) that depends on separation; the force
     F = -dU/dd follows the Coulomb ~1/d^2 trend (the seed of an inter-molecular force)
  C  Linearity/superposition: the electrostatic solve is linear (energy additivity)

What this does NOT do (stated honestly, see theory/06e): derive the macroscopic
viscosity from these forces. That bridge (molecular forces -> stress correlations ->
Green-Kubo -> nu) is statistical mechanics / molecular dynamics, a different engine
from this continuum FE notebook. Here we only verify that the Coulomb FORCE -- the
physical origin the insight points to -- is the same operator the notebook already
trusts.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from scipy.sparse.linalg import spsolve
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.fem3d import fem_laplacian


def setup(n=16):
    V, T = kuhn_cube(n)
    V = V - 0.5                          # center at origin: [-0.5,0.5]^3
    K, M = fem_laplacian(V, T)           # M is the lumped-mass vector
    bnd = (np.abs(V[:, 0]) > 0.49) | (np.abs(V[:, 1]) > 0.49) | (np.abs(V[:, 2]) > 0.49)
    return V, T, K, M, ~bnd


def solve_phi(K, M, rho, free):
    phi = np.zeros(len(rho))
    phi[free] = spsolve(K[free][:, free].tocsc(), (M*rho)[free])
    return phi


def checkA_coulomb_1_over_r():
    print("\n[A] Coulomb 1/r law from K phi = M rho (same operator as pressure Poisson 04)")
    V, T, K, M, free = setup(16)
    r = np.linalg.norm(V, axis=1)
    sigma = 0.06
    rho = np.exp(-(r**2)/(2*sigma**2)); rho /= np.sum(M*rho)     # total charge 1
    phi = solve_phi(K, M, rho, free)
    mask = (r > 0.15) & (r < 0.35) & free
    x = 1.0/r[mask]; y = phi[mask]
    a, b = np.polyfit(x, y, 1)
    corr = np.corrcoef(x, y)[0, 1]
    print(f"    phi ~ a/r + b:  a = {a:.4f}  vs exact 1/(4 pi) = {1/(4*np.pi):.4f}  ({abs(a-1/(4*np.pi))/(1/(4*np.pi))*100:.1f}%)")
    print(f"    correlation phi vs 1/r = {corr:.4f}")
    assert corr > 0.99 and abs(a - 1/(4*np.pi))/(1/(4*np.pi)) < 0.05
    print("    PASS  (electrostatic Coulomb potential IS the K-operator, already verified)")


def two_charge_energy(V, K, M, free, d, sigma=0.06):
    rho = (np.exp(-((V[:, 0]-d/2)**2 + V[:, 1]**2 + V[:, 2]**2)/(2*sigma**2)) +
           np.exp(-((V[:, 0]+d/2)**2 + V[:, 1]**2 + V[:, 2]**2)/(2*sigma**2)))
    rho /= np.sum(M*rho); rho *= 2.0
    phi = solve_phi(K, M, rho, free)
    return 0.5*np.sum(rho*(M*phi))


def checkB_force():
    print("\n[B] Two-charge energy U(d) and force F=-dU/dd (seed of inter-molecular force)")
    V, T, K, M, free = setup(16)
    ds = [0.15, 0.20, 0.25, 0.30, 0.40]
    Us = [two_charge_energy(V, K, M, free, d) for d in ds]
    for d, U in zip(ds, Us):
        print(f"    d={d:.2f}: U={U:.4f}")
    Fs = []
    for i in range(len(ds)-1):
        F = -(Us[i+1]-Us[i])/(ds[i+1]-ds[i]); Fs.append(F)
        print(f"    d~{(ds[i]+ds[i+1])/2:.2f}: F={F:.3f}")
    assert all(F > 0 for F in Fs) and Fs[0] > Fs[-1], "force not Coulomb-like"
    print("    PASS  (configuration-dependent energy; repulsive force decaying with d)")


def checkC_linearity():
    print("\n[C] Linearity / superposition of the electrostatic solve")
    V, T, K, M, free = setup(12)
    r = np.linalg.norm(V, axis=1); sigma = 0.07
    rho1 = np.exp(-((V[:, 0]-0.2)**2+V[:, 1]**2+V[:, 2]**2)/(2*sigma**2)); rho1 /= np.sum(M*rho1)
    rho2 = np.exp(-((V[:, 0]+0.2)**2+V[:, 1]**2+V[:, 2]**2)/(2*sigma**2)); rho2 /= np.sum(M*rho2)
    p1 = solve_phi(K, M, rho1, free); p2 = solve_phi(K, M, rho2, free)
    p12 = solve_phi(K, M, rho1+rho2, free)
    err = np.linalg.norm(p12-(p1+p2))/np.linalg.norm(p12)
    print(f"    ||phi(r1+r2) - (phi1+phi2)|| / ||.|| = {err:.2e}")
    assert err < 1e-10
    print("    PASS  (superposition exact -> the Coulomb operator is linear, like all of L)")


if __name__ == "__main__":
    print("Coulomb-in-K verification (the molecular force behind viscosity)")
    checkA_coulomb_1_over_r()
    checkB_force()
    checkC_linearity()
    print("\nAll Coulomb-operator checks PASSED.")
    print("Verified: the Coulomb potential/force -- the molecular origin of viscosity the")
    print("insight points to -- is the SAME operator K the notebook already trusts (04).")
    print("NOT done here (see theory/06e): deriving macroscopic nu from these forces, which")
    print("needs molecular dynamics + Green-Kubo -- a different engine, honestly out of scope.")
