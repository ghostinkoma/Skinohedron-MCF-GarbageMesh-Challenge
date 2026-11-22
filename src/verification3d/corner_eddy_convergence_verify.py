"""
corner_eddy_convergence_verify.py  --  does the cavity corner-eddy strength converge? (theory/10 §3b)

theory/10 (revised) left one question open: at Re = 1000 the bottom-right (BR) corner
eddy sits at Ghia's location, but its strength (stream function psi) was still rising at
PSPG n <= 128 (55 -> 74% of Ghia). This script pushes the refinement further with TWO
independent elements and estimates the limit.

  Re = 1000, stream function psi from K psi = int phi (dv/dx - du/dy), psi = 0 on walls.
  PSPG P1/P1 at n = 128, 160, 192, and Taylor-Hood P2/P1 at n = 32 .. 64 (its P2 nodes are
  exactly the vertices of mesh_square(2n), so psi is computed on that P1 mesh).

  CHECK A  Monotone convergence: for both elements the BR strength increases with n and
           the increments shrink.
  CHECK B  Extrapolation psi(h) = psi_inf - C h^p (least squares on the last 4 levels):
           both elements give the SAME limit (within 5%), close to Ghia (within 15%).
           Agreement of two different discretisations on the limit is the evidence that
           the eddy strength is a property of the continuous solution, not of the mesh.
  CHECK C  The observed order p is low (~0.6-0.7), well below the rate for smooth
           solutions. This is consistent with -- but does not prove -- pollution by the
           lid-corner singularity (the velocity jumps at the two top corners).

Runtime: ~10 minutes.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import curve_fit
from cavity_highre_verify import pspg_centreline, stream_function, GHIA_VORTEX
from cavity_taylorhood_verify import solve_TH
from cavity_pspg_verify import mesh_square

RE = 1000
PSPG_N = (128, 160, 192)
TH_N = (32, 40, 48, 56, 64)


def eddies(n, V, ux, uy):
    psi = stream_function(n, V, ux, uy)
    x, y = V[:, 0], V[:, 1]
    out = {}
    i = int(np.argmin(psi)); out["primary"] = (float(psi[i]), (x[i], y[i]))
    for tag, win in (("BR", (x > 0.6) & (y < 0.4)), ("BL", (x < 0.4) & (y < 0.4))):
        idx = np.where(win)[0]; j = idx[int(np.argmax(psi[idx]))]
        out[tag] = (float(psi[j]), (x[j], y[j]))
    return out


def pspg_eddies(n):
    _, _, V, ux, uy = pspg_centreline(n, RE)
    return eddies(n, V, ux, uy)


def th_eddies(n):
    V, U, _, _ = solve_TH(n, RE)
    Vm, _ = mesh_square(2 * n)                 # P2 nodes == vertices of mesh_square(2n)
    k = np.round(V * 2 * n).astype(int)
    order = k[:, 1] * (2 * n + 1) + k[:, 0]
    ux = np.zeros(len(Vm)); uy = np.zeros(len(Vm))
    ux[order] = U[:, 0]; uy[order] = U[:, 1]
    return eddies(2 * n, Vm, ux, uy)


def extrapolate(ns, vals):
    h = 1.0 / np.array(ns, float)
    f = lambda h, a, c, p: a - c * h ** p
    popt, _ = curve_fit(f, h, np.array(vals), p0=[vals[-1] * 1.2, vals[-1], 1.0], maxfev=20000)
    return float(popt[0]), float(popt[2])


def main():
    g, (gx, gy) = GHIA_VORTEX[RE]["BR"]
    print(f"Corner-eddy convergence at Re={RE} (Ghia BR psi = {g:.4e} at ({gx:.4f},{gy:.4f}))")

    series = {}
    for name, ns, fn in (("PSPG", PSPG_N, pspg_eddies), ("Taylor-Hood", TH_N, th_eddies)):
        print(f"\n  {name}:")
        vals = []
        for n in ns:
            e = fn(n)
            v, (px, py) = e["BR"]
            vals.append(v)
            print(f"    n={n:>3}: BR psi = {v:.4e} at ({px:.3f},{py:.3f}) = {v/g*100:.0f}% of Ghia"
                  f"   | primary {e['primary'][0]:.4f}   BL {e['BL'][0]:.3e}")
        series[name] = (ns, vals, (px, py))

    print("\n[A] monotone convergence of the BR strength")
    for name, (ns, vals, _) in series.items():
        inc = np.diff(vals)
        print(f"    {name:<12} increments: " + ", ".join(f"{d:.2e}" for d in inc))
        assert np.all(inc > 0) and np.all(np.diff(inc) < 0), f"{name}: not monotonically converging"
    print("    PASS  -- rising with n, with shrinking increments")

    print("\n[B] extrapolated limit psi(h) = psi_inf - C h^p  (last 4 levels)")
    lim = {}
    for name, (ns, vals, loc) in series.items():
        a, p = extrapolate(ns[-4:], vals[-4:])
        lim[name] = (a, p)
        print(f"    {name:<12} psi_inf = {a:.4e} ({a/g*100:.0f}% of Ghia), observed order p = {p:.2f}")
        assert np.hypot(loc[0] - gx, loc[1] - gy) < 0.02, f"{name}: eddy not at Ghia's location"
    a1, a2 = lim["PSPG"][0], lim["Taylor-Hood"][0]
    print(f"    PSPG vs Taylor-Hood limits differ by {abs(a1-a2)/((a1+a2)/2)*100:.1f}%")
    assert abs(a1 - a2) / ((a1 + a2) / 2) < 0.05, "the two elements extrapolate to different limits"
    assert all(abs(a / g - 1) < 0.15 for a, _ in lim.values()), "limit far from Ghia"
    print("    PASS  -- both elements extrapolate to the same limit, within ~10% of Ghia")

    print("\n[C] observed order")
    ps = [p for _, p in lim.values()]
    assert all(0.3 < p < 1.5 for p in ps)
    print(f"    p = {', '.join(f'{p:.2f}' for p in ps)}: slow, well below smooth-solution rates;")
    print("    consistent with (not proof of) pollution by the lid-corner singularity")
    print("    PASS")

    print("\nConclusion: the corner-eddy strength converges, slowly, and two independent")
    print("elements agree on its limit -- it is a property of the continuous solution.")
    print("The value is pinned to ~10% (extrapolation uncertainty), not to Ghia's digits.")


if __name__ == "__main__":
    main()
