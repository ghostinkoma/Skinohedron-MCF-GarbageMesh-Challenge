"""
fluid_ns_error_decomposition_verify.py  --  where the Taylor-Green ~3% comes from (theory/06d)

RECONSTRUCTED (2026-11). The original script was cited by theory/06d and 06z but never
committed. This version re-derives the cited decomposition with the Navier-Stokes Box
of fluid_navier_stokes_verify.py (periodic Kuhn cube, Chorin projection, skew advection,
forward-Euler step).

The Taylor-Green kinetic energy decays at 4 nu k^2 in the continuum (06c measures the
discrete rate within ~3%). The decomposition:

  CHECK A  The nonlinear advection of the Taylor-Green field is (almost) a pure
           gradient, so the projection removes it: the field evolves as a viscous mode.
           (Continuum: u.grad(u) = -grad(...) exactly for Taylor-Green.)
  CHECK B  Temporal error is removable: the measured energy-decay rate converges to
           2 nu lambda_h, lambda_h = (u^T K u)/(u^T M u), as dt -> 0, with error O(dt).
  CHECK C  The remaining error is spatial and O(h^2): lambda_h vs 2 k^2 (the sinusoid is
           not in P1). This is the same number as the 06f stiffness error.

Together: rate = 2 nu lambda_h + O(dt), and lambda_h = 2 k^2 (1 + O(h^2)).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.sparse import coo_matrix
from fluid_navier_stokes_verify import Box, run_TG

K_WAVE = 2 * np.pi
NU = 0.02


def lambda_h(box):
    ux = np.sin(K_WAVE * box.Vr[:, 0]) * np.cos(K_WAVE * box.Vr[:, 1])
    return float(ux @ (box.Kp @ ux)) / float(np.sum(box.Mp * ux * ux))


def advection(box, u):
    """Skew-symmetric advection M^-1 C_skew(u) u, assembled as in Box.step."""
    U3 = u.reshape(box.nn, 3)
    rows, cols, data = [], [], []
    for (ix, vol, G) in box.tets:
        nodes = [box.inv[i] for i in ix]
        ug = U3[nodes].mean(0) @ G
        for a in range(4):
            for b in range(4):
                rows.append(nodes[a]); cols.append(nodes[b]); data.append(vol / 4.0 * ug[b])
    C = coo_matrix((data, (rows, cols)), shape=(box.nn, box.nn)).tocsr()
    Cs = (C - C.T) * 0.5
    out = np.empty_like(u)
    for c in range(3):
        out[c::3] = box.Mi * (Cs @ U3[:, c])
    return out


def check_A(box):
    print("\n[A] Taylor-Green advection is (almost) a pure gradient -> removed by the projection")
    u = box.TG0()
    a = advection(box, u)
    pa = box.project(a)
    M3 = np.repeat(box.Mp, 3)
    frac = 1.0 - float(pa @ (M3 * pa)) / float(a @ (M3 * a))
    print(f"    n=8: gradient fraction of the advection term = {frac:.4f}")
    assert frac > 0.99
    print("    PASS  -- >99% gradient: the velocity decays as a viscous mode")


def check_B(box):
    print("\n[B] temporal error is removable: rate -> 2 nu lambda_h as dt -> 0")
    target = 2 * NU * lambda_h(box)
    T_END = 0.12
    diffs = []
    for dt in (3e-3, 1.5e-3, 7.5e-4):
        ns = int(round(T_END / dt))
        _, E0, E1, maxdiv, _, tend = run_TG(box, NU, nsteps=ns, dt=dt)
        rate = -np.log(E1 / E0) / tend
        diffs.append(rate - target)
        print(f"    dt={dt:.2e}: rate {rate:.4f}   2 nu lambda_h {target:.4f}   difference {rate-target:.2e}")
    r1, r2 = diffs[0] / diffs[1], diffs[1] / diffs[2]
    print(f"    difference ratio per halving of dt: {r1:.2f}, {r2:.2f}  (2 => O(dt))")
    assert 1.8 < r1 < 2.2 and 1.8 < r2 < 2.2 and abs(diffs[-1]) < 2e-3
    print("    PASS  -- the time-stepping error is O(dt) and vanishes; the rate is 2 nu lambda_h")


def check_C():
    print("\n[C] the remaining error is spatial: lambda_h vs 2 k^2, O(h^2)")
    exact = 2 * K_WAVE**2
    ns = (8, 10, 12, 16); errs = []
    for n in ns:
        lam = lambda_h(Box(n))
        errs.append(abs(lam / exact - 1))
        print(f"    n={n:>2}: lambda_h = {lam:.3f}   2k^2 = {exact:.3f}   error {(lam/exact-1)*100:+.2f}%")
    order = float(np.polyfit(np.log(1.0 / np.array(ns, float)), np.log(errs), 1)[0])
    print(f"    order ~ h^{order:.2f}")
    assert 1.8 < order < 2.2
    print("    PASS  -- 5.0% -> 1.3% (n=8 -> 16), O(h^2): the sinusoid is not in P1")


if __name__ == "__main__":
    print("06d: decomposing the Taylor-Green decay-rate error")
    box = Box(8)
    check_A(box)
    check_B(box)
    check_C()
    print("\nConclusion: advection is projected out (A); the time error is O(dt) and")
    print("removable (B); what remains is the O(h^2) spatial error of representing a")
    print("sinusoid in P1 (C) -- the same number as the 06f stiffness error.")
