"""
rev_origin_verify.py  --  what the cavity "rev" is made of, re-examined (theory/10)

RECONSTRUCTED AND REVISED (2026-11). The original script was cited by theory/10, 11
and the precision viewer but never committed. Its central "physical rev converges"
result used a u_min proxy that measures the primary vortex's return flow, not the
corner eddy (issue #3, A3). This version re-measures everything with the stream
function and states what the data does and does not support.

The "rev" of 08b has two parts:
  NUMERICAL rev -- the discretisation's disagreement with Ghia and between the two
                   elements (PSPG P1/P1 vs Taylor-Hood P2/P1);
  PHYSICAL rev  -- the Reynolds-driven flow structure: secondary corner eddies and
                   the migration of the primary vortex.

  CHECK A  Numerical rev shrinks under refinement (Re = 1000): Ghia centreline RMS for
           PSPG (n = 64, 96, 128) and Taylor-Hood (n = 24, 32, 40), and the mutual
           PSPG-TH RMS for the pairs (96, 32), (112, 36), (128, 40). All decrease:
           this part is approximation error.
  CHECK B  Physical rev, stream function: the bottom-right (BR) eddy at Re = 1000 sits
           at Ghia's location at every resolution, but its strength is still RISING
           (55% -> 67% -> 74% of Ghia at n = 64, 96, 128). Its location is resolved;
           its strength is NOT converged at these meshes.
  CHECK C  The eddy grows with Re (equivalently with the advection weight alpha:
           alpha u.grad u - nu Lap u is Re -> alpha Re): at n = 96 it is unresolved
           at Re <= 100 (Ghia's Re = 100 eddy, psi ~ 1e-5, is below this mesh), appears
           at Re = 400 and grows x2.8 by Re = 1000 (Ghia x2.7).

Revised conclusion (replacing the original "physical rev converges -> field-intrinsic"):
the numerical rev is approximation error and shrinks; the physical rev is real flow
structure (eddies where Ghia finds them, growing with Re), but at these resolutions its
strength is still approximation-limited. Pure incompressible NS (Ghia's model) is the
reference throughout, so no extra physics is needed to explain the cavity.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from cavity_highre_verify import (pspg_centreline, th_centreline, stream_function,
                                  GHIA, GHIA_Y, GHIA_VORTEX)

RE = 1000


def rms(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))


def br_eddy(n, V, ux, uy):
    psi = stream_function(n, V, ux, uy)
    x, y = V[:, 0], V[:, 1]
    idx = np.where((x > 0.6) & (y < 0.4))[0]
    j = idx[int(np.argmax(psi[idx]))]
    return float(psi[j]), (float(x[j]), float(y[j]))


def check_A_and_B():
    print(f"\n[A] numerical rev at Re={RE}: Ghia RMS and PSPG-vs-TH disagreement under refinement")
    P, T, BR = {}, {}, {}
    for n in (64, 96, 112, 128):
        y, u, V, ux, uy = pspg_centreline(n, RE)
        P[n] = np.interp(GHIA_Y, y, u)
        if n != 112:
            BR[n] = br_eddy(n, V, ux, uy)
    for n in (24, 32, 36, 40):
        y, u = th_centreline(n, RE)
        T[n] = np.interp(GHIA_Y, y, u)
    gp = [rms(P[n], GHIA[RE]) for n in (64, 96, 128)]
    gt = [rms(T[n], GHIA[RE]) for n in (24, 32, 40)]
    mu = [rms(P[a], T[b]) for a, b in ((96, 32), (112, 36), (128, 40))]
    print("    PSPG  n=64/96/128 vs Ghia RMS: " + " / ".join(f"{v:.4f}" for v in gp))
    print("    TH    n=24/32/40  vs Ghia RMS: " + " / ".join(f"{v:.4f}" for v in gt))
    print("    PSPG-TH mutual RMS (96,32)/(112,36)/(128,40): " + " / ".join(f"{v:.4f}" for v in mu))
    for s in (gp, gt, mu):
        assert s[0] > s[1] > s[2], "numerical rev does not shrink"
    print("    PASS  -- every measure of the numerical rev shrinks: approximation error")

    print(f"\n[B] physical rev at Re={RE}: bottom-right eddy (stream function) under refinement")
    g, (gx, gy) = GHIA_VORTEX[RE]["BR"]
    vals = []
    for n in (64, 96, 128):
        v, (px, py) = BR[n]
        vals.append(v)
        print(f"    n={n:>3}: psi = {v:.3e} at ({px:.3f},{py:.3f})   = {v/g*100:.0f}% of Ghia"
              f" {g:.3e} at ({gx:.4f},{gy:.4f})")
        assert np.hypot(px - gx, py - gy) < 0.05, "eddy not at Ghia's location"
    d1, d2 = vals[1] - vals[0], vals[2] - vals[1]
    print(f"    successive increments: {d1:.2e}, {d2:.2e} (still rising)")
    assert vals[0] < vals[1] < vals[2] < g
    print("    PASS  -- location resolved; strength still rising toward Ghia: NOT converged")


def check_C():
    print("\n[C] the eddy grows with Re (= advection weight alpha) -- n=96, stream function")
    out = {}
    for Re in (100, 400, 1000):
        y, u, V, ux, uy = pspg_centreline(96, Re)
        out[Re] = br_eddy(96, V, ux, uy)
        v, (px, py) = out[Re]
        print(f"    Re={Re:>4}: BR psi max = {v:.3e} at ({px:.3f},{py:.3f})")
    assert out[100][0] < 1e-6, "expected the Re=100 eddy to be below this mesh"
    assert out[400][0] > 0 and out[1000][0] > 2.0 * out[400][0]
    g = GHIA_VORTEX
    print(f"    growth 400 -> 1000: x{out[1000][0]/out[400][0]:.1f}"
          f"  (Ghia x{g[1000]['BR'][0]/g[400]['BR'][0]:.1f}); Re=100 eddy (Ghia ~1.3e-5) unresolved here")
    print("    PASS  -- the eddy strengthens with Re, as in Ghia")


if __name__ == "__main__":
    print("10 (revised): separating the cavity 'rev' into numerical and physical parts")
    check_A_and_B()
    check_C()
    print("\nRevised conclusion: the numerical rev shrinks with refinement (approximation).")
    print("The physical rev is real structure -- eddies at Ghia's locations, growing with")
    print("Re -- but its strength is still approximation-limited at n <= 128. The original")
    print("'physical rev converges -> field-intrinsic' finding does not hold; it came from")
    print("a proxy that measured the primary return flow.")
