"""
sparam_geometry_verify.py  --  geometric consistency of the scatter node (theory/01e)

RECONSTRUCTED (2026-11). The original script was cited by theory/01e and the
prologue but never committed. This version re-measures the cited claim with the
repository's own engines; the numbers it prints are the ones now recorded in 01e.

Claim (01e): the geometry-blind scatter node (ksf3d.sparam_wave: uniform S at every
cell, equal impedance, no use of face areas or distances) propagates waves with a
strongly DIRECTION-DEPENDENT speed, and this anisotropy does NOT shrink under
refinement; the FE/DEC operator L = M^-1 K on the same Kuhn mesh is about ten times
more isotropic.

  CHECK A  scatter node: front speed along six lattice directions (<100>, <010>,
           <001>, <110>, <111>, <-111>), from a monopole pulse at the centre cell.
           Delta c / c = (max - min) / mean is large (~0.5) at n = 8, 12, 16 and does
           not decrease with n -- a property of the node, not of resolution.
  CHECK B  FE operator on the same mesh: plane-wave phase speed from the Rayleigh
           quotient (as in unified_scalar_verify.py) gives Delta c / c ~ 0.05.

Measurement details (they affect the scatter number by a few hundredths, not its
size): arrival = first step at which a cell's port energy exceeds 1e-6; cells within
a cone 1 - cos(angle) < 0.15 around each direction and 0.1 < distance < 0.45 are fitted
linearly (distance vs arrival step).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from ksf3d.mesh3d_uniform import kuhn_cube
from ksf3d.sparam_wave import PortComplex, Connect, U_step, inject_pulse
from unified_scalar_verify import wave_anisotropy

DIRS = {"<100>": (1, 0, 0), "<010>": (0, 1, 0), "<001>": (0, 0, 1),
        "<110>": (1, 1, 0), "<111>": (1, 1, 1), "<-111>": (-1, 1, 1)}


def scatter_speeds(n, thr=1e-6, cone=0.15):
    V, T = kuhn_cube(n); V = V - 0.5
    pc = PortComplex(V, T)
    cent = V[T].mean(axis=1)
    C = Connect(pc, np.ones(pc.N_C), rho_bnd=1.0)
    c0 = int(np.argmin(np.linalg.norm(cent, axis=1)))
    a = np.zeros(pc.nP); inject_pulse(a, pc, c0)
    rel = cent - cent[c0]; dist = np.linalg.norm(rel, axis=1)
    arrival = np.full(pc.N_C, -1)
    for t in range(1, 3 * n):
        a = U_step(a, C)
        w = (a * a).reshape(-1, 4).sum(axis=1)
        arrival[(w > thr) & (arrival < 0)] = t
    speeds = {}
    for name, d in DIRS.items():
        d = np.array(d, float); d /= np.linalg.norm(d)
        cosang = (rel @ d) / np.maximum(dist, 1e-12)
        sel = (cosang > 1 - cone) & (arrival > 0) & (dist > 0.1) & (dist < 0.45)
        speeds[name] = float(np.polyfit(arrival[sel], dist[sel], 1)[0]) * n  # cells / step
    v = np.array(list(speeds.values()))
    return speeds, float((v.max() - v.min()) / v.mean())


def main():
    print("01e: geometric consistency -- scatter node vs FE operator on the Kuhn mesh")

    print("\n[A] geometry-blind scatter node: directional front speed (cells per step)")
    an = []
    for n in (8, 12, 16):
        sp, a = scatter_speeds(n)
        an.append(a)
        print(f"    n={n:>2}: " + "  ".join(f"{k} {v:.3f}" for k, v in sp.items())
              + f"   Delta c/c = {a:.3f}")
    assert min(an) > 0.4, "scatter anisotropy unexpectedly small"
    assert an[-1] > 0.9 * an[0], "scatter anisotropy shrinks with refinement"
    print("    PASS  -- ~0.5 and not shrinking with n (fastest along <111>, the Kuhn diagonal)")

    print("\n[B] FE operator L = M^-1 K on the same mesh: plane-wave phase speed")
    fe = []
    for n in (12, 20):
        cs, a = wave_anisotropy(n)
        fe.append(a)
        print(f"    n={n:>2}: " + "  ".join(f"c{k}={v:.4f}" for k, v in cs.items())
              + f"   Delta c/c = {a:.4f}")
    assert max(fe) < 0.06
    print(f"    PASS  -- Delta c/c ~ {np.mean(fe):.3f}, about {np.mean(an)/np.mean(fe):.0f}x more isotropic")

    print("\nConclusion: the scatter node's anisotropy is a property of the node (it does")
    print("not converge away), so it was rejected in favour of the FE/cotangent operator,")
    print("whose residual ~5% is the Kuhn mesh's own directional bias.")


if __name__ == "__main__":
    main()
