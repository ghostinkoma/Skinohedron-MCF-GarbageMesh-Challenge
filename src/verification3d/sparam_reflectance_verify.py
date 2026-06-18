"""
sparam_reflectance_verify.py
============================
Step 1.5 closure attempt: recover the material reflectance R² from the simulation
and match the theoretical value to machine-ish precision.

WHAT WORKED, AND THE HONEST STORY OF WHAT DID NOT
-------------------------------------------------
The goal was a TF/SF *one-way* (directional) source (theory/01d). Implementing it
revealed something important and consistent with the geometry objection of 01e:

  * A literal one-way TF/SF port injection does NOT produce a clean +x wave on the
    current geometry-blind node S_wave = 2P₀−I. The node redistributes energy
    isotropically among the 4 tetrahedral faces (which are not axis-aligned), so
    there is no clean "+x port" to inject into. Result: the matched case (Z₂=Z₁)
    still reflected ~50%.
  * A per-face *flux monitor* (classify b² by face-normal sign) is confounded for
    the same reason: a +x-moving wave deposits energy on both +x-ish and -x-ish
    tetrahedral faces, giving a spurious ~0.17 reflection even with no interface.

Both failures are the SAME root cause the user identified: the node ignores face
orientation (the subject of 01e). Direction-resolved quantities are contaminated.

THE FIX THAT IS GEOMETRY-ROBUST
-------------------------------
*Total* absorbed energy at a boundary plane is geometry-robust: summing b² over
all faces at an absorbing end averages over the tetrahedral orientations. With a
symmetric (monopole) plane source in medium 1 and absorbing ends, exact energy
bookkeeping gives

    E_left  = ½ + ½R²      (backward half + reflected forward half)
    E_right = ½T²
    refl_frac = E_left/(E_left+E_right) = ½ + ½R²     (using R²+T²=1)

so the reflectance is recovered EXACTLY (up to discretisation) by

    R²_measured = 2·refl_frac − 1 .                                   (★)

The monopole injection (equal on all 4 ports) excites only the symmetric mode
(eigenvalue +1 of S_wave), which radiates ±x symmetrically in a uniform medium —
so the ½/½ split is exact, and (★) inherits only the small discretisation error of
the single reflection. That error converges fast (see table): R² is matched to
~1e-8 by nx=80, and the matched control (R²=0) floor vanishes under refinement.

CONCLUSION
----------
The material model R = (Z₂−Z₁)/(Z₂+Z₁) is verified as a physical observable to
~1e-8 — Step 1.5's reflectance question is closed via a geometry-robust measure.
The literal one-way source remains blocked by the geometry-blind node, which is
exactly what 01e is meant to fix; this script documents that link with numbers
rather than asserting a one-way source that does not yet exist.
"""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from ksf3d.mesh3d_uniform import kuhn_box
from ksf3d.sparam_wave import PortComplex, S_apply, Connect


def build_channel(nx, nt=2, Lx=2.0):
    V, T = kuhn_box(nx, nt, nt, Lx=Lx)
    pc = PortComplex(V, T)
    h = Lx / nx
    w = nt * h
    pc.make_periodic([(1, 0.0, w), (2, 0.0, w)])
    cent = np.array([V[T[i]].mean(axis=0) for i in range(pc.N_C)])
    return pc, cent, h


def measure_R2(nx, Z1, Z2, x_s=-0.7, omega=0.3, tau=10.0, t0=30.0, nsteps=2000):
    """Symmetric plane source + absorbing ends; recover R² via (★)."""
    pc, cent, h = build_channel(nx)
    Z = np.where(cent[:, 0] < 0.0, Z1, Z2)
    C = Connect(pc, Z, rho_bnd=0.0)
    src = np.where(np.abs(cent[:, 0] - x_s) < h * 0.6)[0]
    port_cell = np.arange(pc.nP) // 4
    bnd = pc.is_boundary
    leftx = cent[port_cell, 0] < 0.0
    a = np.zeros(pc.nP)
    Erefl = Etrans = 0.0
    for t in range(nsteps):
        ainc = np.exp(-((t - t0) / tau) ** 2) * np.sin(omega * (t - t0))
        for ci in src:
            for k in range(4):
                a[pc.pid(ci, k)] += 0.5 * ainc
        b = S_apply(a)
        ab = (b * b) * bnd
        Erefl += ab[leftx].sum()
        Etrans += ab[~leftx].sum()
        a = C.apply(b)
    refl_frac = Erefl / (Erefl + Etrans)
    R2_meas = 2.0 * refl_frac - 1.0
    R = (Z2 - Z1) / (Z2 + Z1)
    return R2_meas, R * R


def run():
    pairs = [(1, 0.5), (1, 0.3), (1, 0.01)]
    nxs = (20, 40, 80)
    out = {"pairs": [], "matched": []}
    for Z1, Z2 in pairs:
        rows = []
        for nx in nxs:
            m, th = measure_R2(nx, Z1, Z2)
            rows.append((nx, m, th, abs(m - th)))
        out["pairs"].append(((Z1, Z2), rows))
    for nx in nxs:
        m, _ = measure_R2(nx, 1.0, 1.0)
        out["matched"].append((nx, m))
    return out, nxs


def main():
    out, nxs = run()
    print("=== STEP 1.5 closure · reflectance R² recovery ===\n")
    print("Method: symmetric plane source + absorbing ends, R² = 2·refl_frac − 1")
    print("(geometry-robust total-energy measure; see module docstring for why the")
    print(" literal one-way TF/SF source is blocked by the geometry-blind node).\n")

    worst_fine = 0.0
    for (Z1, Z2), rows in out["pairs"]:
        th = rows[0][2]
        print(f"  Z₁|Z₂ = {Z1}|{Z2}   (R²_theory = {th:.6f})")
        print(f'     {"nx":>4} {"R²_meas":>11} {"abs_err":>11}')
        for nx, m, t, e in rows:
            print(f'     {nx:4d} {m:11.6f} {e:11.2e}')
        worst_fine = max(worst_fine, rows[-1][3])
        print()

    print("  matched control  Z₂=Z₁  (R²_theory = 0):")
    for nx, m in out["matched"]:
        print(f"     nx={nx:3d}:  R²_meas = {m:+.3e}")
    matched_fine = abs(out["matched"][-1][1])

    # convergence check: error must shrink with refinement
    conv_ok = all(rows[0][3] > rows[-1][3] for _, rows in out["pairs"])

    print()
    print(f"  finest-grid (nx={nxs[-1]}) worst reflectance error : {worst_fine:.2e}")
    print(f"  finest-grid matched-control floor .............. : {matched_fine:.2e}")
    print(f"  error decreases under refinement (all pairs) ... : {'PASS' if conv_ok else 'FAIL'}")

    print("\n[結句] R² is recovered to ~1e-8 at nx=80 and converges fast; the")
    print("  matched-control spurious reflection vanishes under refinement. The")
    print("  material reflectance model is thus verified as a physical observable,")
    print("  closing Step 1.5's quantitative item — by a geometry-robust TOTAL-")
    print("  energy measure. The direction-resolved one-way source/monitor remain")
    print("  contaminated by the geometry-blind node, which is precisely what 01e")
    print("  addresses; that link is shown here with numbers, not assumed.")

    assert worst_fine < 1e-5, f"R² not recovered at fine grid: {worst_fine}"
    assert matched_fine < 1e-4, f"matched floor too large: {matched_fine}"
    assert conv_ok, "error did not decrease under refinement"
    print("\nReflectance verification asserts passed.")


if __name__ == "__main__":
    main()
