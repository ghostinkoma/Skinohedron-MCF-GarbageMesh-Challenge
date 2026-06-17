"""
sparam_interface_verify.py
==========================
Step 1.5 verification (theory/01c_planewave_interface.md): does the material
model (R,T) of 01b §4 reproduce as a measured physical observable?

This script is deliberately split into what the code CAN confirm and what it
CANNOT yet — honesty over a green checkmark. Running it on your machine
reproduces both halves.

PART A — periodic plane-wave channel machinery (EXACT, asserted)
  The quasi-1D channel (long in x, periodic in y,z, absorbing x-ends) is built
  and its lossless core is verified to machine precision:
    * partner involution survives periodic re-pairing  (Π∘Π = id)
    * with reflecting ends, energy is conserved exactly (U = C S orthogonal)
  These confirm the 01c channel construction itself is sound.

PART B — quantitative R²/T² recovery (HONEST DIAGNOSTIC, not yet a clean pass)
  We launch a wave in medium 1 and account absorbed energy at the two ends.
  Findings (reproduced by the run below):
    * QUALITATIVELY correct: reflection rises monotonically as the contrast
      grows; a metal-like wall (Z2≪Z1) reflects ~96% ≈ R²; the matched case is
      the minimum.
    * QUANTITATIVELY imperfect with a simple source: there is a spurious
      reflection floor (nonzero reflection even at Z2=Z1, where R²=0), because a
      naive port injection is not a clean +x eigenmode — it radiates some energy
      backwards, and not all injected energy reaches the ends (energy-accounting
      < 1).
  CONCLUSION (recorded, not hidden): the algebraic identity R²+T²=1 is exact
  (already asserted in sparam_wave_verify §G), and the interface behaves
  correctly in trend, but a clean numerical match of the *fraction* to R²
  requires a proper directional (total-field/scattered-field) source. That is
  the real open item of Step 1.5 — surfaced precisely because we insisted on
  code instead of taking the theory doc on faith.
"""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from ksf3d.mesh3d_uniform import kuhn_box
from ksf3d.sparam_wave import PortComplex, S_apply, Connect, energy


def build_channel(nx, nt=2, Lx=2.0):
    V, T = kuhn_box(nx, nt, nt, Lx=Lx)
    pc = PortComplex(V, T)
    h = Lx / nx
    w = nt * h
    pc.make_periodic([(1, 0.0, w), (2, 0.0, w)])   # periodic in y, z
    return pc, V, T, h


# --------------------------------------------------------------------------- #
#  PART A — channel machinery, exact                                          #
# --------------------------------------------------------------------------- #
def partA():
    pc, V, T, h = build_channel(nx=20, nt=2)
    n_bnd = int(pc.is_boundary.sum())
    # involution already asserted inside make_periodic; reaching here = pass
    # energy conservation with reflecting ends (lossless)
    Z = np.ones(pc.N_C)
    C = Connect(pc, Z, rho_bnd=1.0)
    rng = np.random.default_rng(0)
    a = rng.normal(size=pc.nP)
    E0 = energy(a)
    for _ in range(400):
        a = C.apply(S_apply(a))
    rel = abs(energy(a) - E0) / E0
    return {"n_boundary_after_periodic": n_bnd, "energy_relerr": rel}


# --------------------------------------------------------------------------- #
#  PART B — interface diagnostic, honest                                      #
# --------------------------------------------------------------------------- #
def measure_interface(nx, Z1, Z2, nt=2, Lx=2.0, nsteps=800):
    pc, V, T, h = build_channel(nx, nt, Lx)
    cent = np.array([V[T[i]].mean(axis=0) for i in range(pc.N_C)])
    Z = np.where(cent[:, 0] < 0.0, Z1, Z2)
    C = Connect(pc, Z, rho_bnd=0.0)               # absorbing ends
    # directional-ish source: excite -x-facing ports of a layer in medium 1
    x0 = -0.7
    layer = np.abs(cent[:, 0] - x0) < h * 0.6
    a = np.zeros(pc.nP)
    nx_comp = pc.normal[:, 0]
    for ci in np.where(layer)[0]:
        for k in range(4):
            p = pc.pid(ci, k)
            if nx_comp[p] < -0.5:
                a[p] = 1.0
    E_in = energy(a)
    port_cell = np.arange(pc.nP) // 4
    bnd = pc.is_boundary
    leftx = cent[port_cell, 0] < 0
    Erefl = Etrans = 0.0
    for _ in range(nsteps):
        b = S_apply(a)
        ab = (b * b) * bnd
        Erefl += ab[leftx].sum()
        Etrans += ab[~leftx].sum()
        a = C.apply(b)
    tot = Erefl + Etrans
    R = (Z2 - Z1) / (Z2 + Z1)
    Tc = 2 * np.sqrt(Z1 * Z2) / (Z1 + Z2)
    return {"refl": Erefl / tot, "trans": Etrans / tot,
            "R2": R * R, "T2": Tc * Tc, "E_acct": tot / E_in}


def partB():
    rows = []
    for Z1, Z2 in [(1, 1.0), (1, 0.5), (1, 0.3), (1, 0.01)]:
        m = measure_interface(40, Z1, Z2)
        rows.append((Z1, Z2, m))
    return rows


def main():
    print("=== STEP 1.5 · plane-wave interface (theory 01c) ===\n")

    print("PART A — periodic plane-wave channel machinery (EXACT):")
    A = partA()
    print(f"  boundary ports after periodic re-pairing (x-ends only): "
          f"{A['n_boundary_after_periodic']}")
    print(f"  involution Π∘Π=id after re-pairing ............ PASS (asserted)")
    okA = A["energy_relerr"] < 1e-12
    print(f"  energy conservation, reflecting ends .......... "
          f"{'PASS' if okA else 'FAIL'}  (rel {A['energy_relerr']:.1e})")
    print("  => the 01c channel construction is sound.\n")

    print("PART B — quantitative R²/T² recovery (HONEST DIAGNOSTIC):")
    rows = partB()
    print(f'   {"Z1|Z2":>10} {"refl":>8} {"R²":>8} {"trans":>8} {"T²":>8} {"E_acct":>8}')
    for Z1, Z2, m in rows:
        print(f'   {Z1:.2f}|{Z2:<5.2f} {m["refl"]:8.3f} {m["R2"]:8.3f} '
              f'{m["trans"]:8.3f} {m["T2"]:8.3f} {m["E_acct"]:8.3f}')

    # qualitative checks (these we DO assert)
    refls = [m["refl"] for _, _, m in rows]               # Z2 = 1,0.5,0.3,0.01
    mono = all(refls[i] <= refls[i + 1] + 1e-9 for i in range(len(refls) - 1))
    metal_ok = rows[-1][2]["refl"] > 0.9                  # metal reflects strongly
    matched_min = abs(refls[0] - min(refls)) < 1e-9       # matched = least reflection
    print()
    print(f"  qualitative: reflection monotonic in contrast .. {'PASS' if mono else 'FAIL'}")
    print(f"  qualitative: metal wall reflects >90% .......... {'PASS' if metal_ok else 'FAIL'}")
    print(f"  qualitative: matched case is the minimum ....... {'PASS' if matched_min else 'FAIL'}")

    floor = rows[0][2]["refl"]   # reflection at Z2=Z1 (should be 0)
    acct = rows[0][2]["E_acct"]
    print()
    print(f"  quantitative: NOT a clean pass. Spurious reflection floor "
          f"= {floor:.3f}")
    print(f"  at the matched interface (R²=0), and energy accounting "
          f"= {acct:.3f} (<1):")
    print("  the simple port source is not a pure +x eigenmode. A clean fraction")
    print("  match to R² needs a total-field/scattered-field directional source.")
    print("  This is the genuine open item of Step 1.5 — surfaced by running the")
    print("  code rather than trusting the theory doc.\n")

    print("[結句] Verified by reproducible code: (A) the 01c plane-wave channel")
    print("  (periodic transverse pairing + absorbing ends) is constructed")
    print("  correctly and is exactly lossless; (B) the interface is")
    print("  QUALITATIVELY correct (monotone contrast, metal≈R², matched=min) but")
    print("  the QUANTITATIVE R² fraction is not yet cleanly recovered with a")
    print("  naive source. The exact identity R²+T²=1 stands (sparam_wave §G).")
    print("  Step 1.5 is therefore PARTIALLY closed: machinery yes, clean")
    print("  R²-fraction measurement still open (directional source needed).")

    # hard gate: only the things we honestly verified
    assert okA, "channel energy conservation failed"
    assert mono and metal_ok and matched_min, "qualitative interface checks failed"
    print("\nAsserted parts passed (channel exactness + qualitative interface).")


if __name__ == "__main__":
    main()
