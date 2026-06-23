# 08 · P0 — A Trustworthy Pressure, Validated on the Lid-Driven Cavity

**Status:** **verified**, `src/verification3d/cavity_pspg_verify.py` (2-D, Ghia benchmark)
and `cavity3d_demo.py` (3-D). This is the load-bearing test named by the audit (`06z`):
the `06z` finding was that the equal-order P1 pressure is inf-sup-unstable / checkerboard
and **not trustworthy**, and that the lid-driven cavity against published data would
either validate or puncture the foundation. It validates — once the pressure is
stabilised. Done in three steps.

---

## 1. Step 1 — pressure as a stabilised domain, vs an external standard

The equal-order P1 Navier–Stokes is **PSPG-stabilised** (the pressure gets a
`−τK` Laplacian term, `τ` the standard cell-Reynolds parameter) and the steady
lid-driven cavity is solved by Picard iteration at **Re = 100**, then compared to the
canonical **Ghia, Ghia & Shin (1982)** centreline data.

| `n` | u-centreline RMS vs Ghia | v-centreline RMS vs Ghia | pressure 2nd-diff |
|----|------|------|------|
| 48 | 0.0095 | 0.0098 | 0.0011 |
| 64 | 0.0069 | 0.0064 | 0.0008 |

Both centrelines **converge to the published benchmark** (`~7×10⁻³` RMS at `n=64`, still
decreasing), and the pressure's second difference — a checkerboard indicator — is **small
and decreasing**, i.e. the pressure is **smooth, not checkerboard**. The `06z`
"pressure-not-trustworthy" gap is closed by stabilisation, and — crucially — the closure
is checked against an **external standard**, not an internal invariant.

This is the first time the project's fluid output is validated against data it did not
itself generate. It is the brick that *could* have failed, and did not.

---

## 2. Step 2 — pressure as a correlation on the one operator

The stabilisation is not a new machine. The PSPG block is exactly **`−τK`**, with the
**same stiffness `K`** that carries heat (`02`/`03`), the scalar wave (`03`), Coulomb
(`06e`), and temperature (`06g`/`06h`). Verified: that `K` is a consistent Laplacian
(row sums `≤10⁻¹⁵`). So the cavity pressure is **another reading of `L = M⁻¹K`** — a
correlation domain on the one tetrahedral operator, joining the domain map of `07`
rather than sitting beside it. Pressure was always one of the operator's readings
(`04`); the cavity shows it is a *trustworthy* one when stabilised.

---

## 3. Step 3 — the same tetrahedron, in 3-D

The identical PSPG assembly on the **Kuhn-cube tetrahedral** mesh solves the **3-D**
lid-driven cavity (`Re=100`, `n=12`). On the central line the streamwise velocity is

`u(z=0.95) = +0.69`, `u(z=0.5) = −0.13`, `u(z=0.2) = −0.08`,

— **positive near the moving lid, negative below (recirculation)** — the same sign
structure as 2-D Ghia, with a smooth pressure. (The 3-D cavity differs *quantitatively*
from 2-D because of side-wall drag; this is a structural demonstration, not a tabulated
benchmark — the validated benchmark is the 2-D case.) The same tetrahedral operator
carries the pressure-driven flow in 2-D and 3-D alike.

---

## 4. What is and isn't claimed

**Verified:**
- PSPG-stabilised P1 reproduces the Ghia (1982) `Re=100` cavity on both centrelines
  (`u,v` RMS `~7×10⁻³` at `n=64`, converging) with a **smooth** pressure.
- The stabilisation operator is the project's own `K`; pressure is a correlation domain
  on `L=M⁻¹K`.
- The same assembly reproduces the 3-D cavity recirculation structure.

**Not claimed:**
- **Only `Re=100` is benchmarked.** Higher Re (400, 1000, 3200) and the secondary
  corner vortices are **not** yet matched — the next rung if this direction continues.
- **PSPG is a stabilisation, not an inf-sup-stable element.** A Taylor–Hood (P2–P1)
  pair would be the "proper" route; PSPG is the pragmatic one and is what is validated
  here. The two should eventually be cross-checked.
- The 3-D result is a **structural demonstration**, not a benchmark (no standard 1-D
  table for the 3-D cavity at this Re without side-wall-specific references).
- No turbulence, no time-dependence; steady laminar only.
- No new mathematics: PSPG, Picard linearisation, and the Ghia benchmark are classical;
  the contribution is closing the project's own pressure gap against external data and
  showing the pressure rides the same operator.

---

## 5. Files

- Theory: this note, on `06c`/`06z`/`07`.
- Verification: `src/verification3d/cavity_pspg_verify.py` (2-D, steps 1–2, Ghia),
  `src/verification3d/cavity3d_demo.py` (3-D, step 3).
- Next rung: `Re=400/1000` cavity and a Taylor–Hood cross-check (see `07a` P0/P1).
