# 08b · Higher Reynolds — Re=400 & 1000, Both Elements, and the "Rev"

> **Known theory.** the corner eddies are the Moffatt (1964) eddies; the benchmark values are Ghia, Ghia & Shin (1982), Table I and Table V. See `theory/README.md`, *Where the results sit in known theory*, for references.

**Status:** **verified**, `src/verification3d/cavity_highre_verify.py`. Pushes the
cross-validated cavity (`08`, `08a`) to `Re=400` and `Re=1000` with **both** the PSPG and
the Taylor–Hood elements, and searches for the **"rev"** — the Reynolds-driven change in
flow structure, and the point where the discretisation begins to strain. Steady laminar
throughout (true unsteady transition is `~Re 8000`, far beyond here).

---

## 1. Both elements vs Ghia (1982)

| Re | PSPG (n=64) centreline RMS | Taylor–Hood (n=32) RMS |
|----|------|------|
| 400 | 0.0311 | 0.0139 |
| 1000 | 0.0555 | 0.0265 |

Both are close to the published benchmark at these resolutions (one mesh per element —
this check does not refine; convergence under refinement is shown at `Re=100` in `08` /
`08a`). Taylor–Hood is markedly more accurate per node (its P2 velocity is the `06f`
enrichment). `Re=1000` is genuinely harder — thinner boundary layers — so the error is
larger, reported honestly rather than hidden.

---

## 2. The cross-check holds, and loosens — the numerical "rev"

PSPG (with `τ`) and Taylor–Hood (no parameter) still agree at higher Re, but the
agreement **loosens as Re rises**:

| Re | PSPG vs Taylor–Hood mutual RMS |
|----|------|
| 400 | 0.0062 |
| 1000 | 0.0161 |

Two independent discretisations agreeing is the evidence the physics is right (`08a`).
The *loosening* (`0.006 → 0.016`) is the **numerical rev**: at `Re=1000` the schemes begin
to strain against the boundary layers, and they start to part company — a signal, visible
in the data, that more resolution or higher order is needed before pushing further.

---

## 3. The physical "rev" — corner vortices and centre migration

> **Correction (issue #3, A3).** An earlier version of this section measured the
> "bottom-right corner vortex" as the minimum of `u` in the window `x>0.8, y<0.2`
> (`−0.0705 → −0.1413`). That minimum lies on the window edge (e.g. `(0.81, 0.19)`) and is
> the **primary vortex's return flow** along the floor, not the secondary eddy, whose
> near-wall `u` is positive. The section below uses the stream function instead.

The structure is now measured with the **stream function** `ψ` of the PSPG solution
(`K ψ = ∫φ(∂v/∂x − ∂u/∂y)`, `ψ = 0` on the walls, `u = ∂ψ/∂y`, `v = −∂ψ/∂x` as in Ghia):
the primary vortex is `min ψ < 0`; the bottom-right (BR) and bottom-left (BL) secondary
eddies are `max ψ > 0` (counter-rotating) in each bottom corner. Reference: Ghia, Ghia &
Shin (1982), Table V.

| Re | n | primary `ψ` (centre) | BR eddy `ψ` (location) | BL eddy `ψ` (location) |
|----|---|------|------|------|
| 400 | 128 | `−0.1058` (0.562, 0.609) | `4.73e-4` (0.891, 0.117) | `2.8e-6` (0.039, 0.039) |
| 400 | Ghia | `−0.1139` (0.5547, 0.6055) | `6.42e-4` (0.8906, 0.1250) | `1.42e-5` (0.0508, 0.0469) |
| 1000 | 64 | `−0.0891` (0.547, 0.578) | `9.59e-4` (0.875, 0.109) | `5.8e-5` (0.078, 0.062) |
| 1000 | 96 | `−0.0966` (0.542, 0.573) | `1.18e-3` (0.875, 0.115) | `1.04e-4` (0.073, 0.073) |
| 1000 | 128 | `−0.1011` (0.531, 0.570) | `1.30e-3` (0.867, 0.109) | `1.31e-4` (0.078, 0.070) |
| 1000 | Ghia | `−0.1179` (0.5313, 0.5625) | `1.75e-3` (0.8594, 0.1094) | `2.31e-4` (0.0859, 0.0781) |

- **The secondary eddies are where Ghia puts them, and they counter-rotate.** At `n=128`
  the BR eddy lies within `0.01` of Ghia's location at both Re; the BL eddy is close at
  `Re=1000` and barely resolved at `Re=400`.
- **They grow with Re:** BR `×2.8` from `Re=400` to `1000` (Ghia `×2.7`); BL grows by more
  than an order of magnitude (Ghia `×16`; our `Re=400` BL value is too under-resolved for a
  reliable ratio).
- **The primary vortex centre migrates toward the geometric centre:** `(0.562,0.609) →
  (0.531,0.570)`, matching Ghia's `(0.5547,0.6055)` and `(0.5313,0.5625)` to about
  `0.01`. This migration is the canonical cavity "rev".
- **The strengths are NOT converged.** Under refinement at `Re=1000` the BR eddy rises
  `55% → 67% → 74%` of Ghia's value (`n = 64, 96, 128`), the BL eddy `25% → 45% → 57%`, and
  the primary `76% → 82% → 86%`. P1/PSPG under-resolves the eddies at these meshes; the
  trend is toward Ghia, but a finer mesh or higher order is needed to pin the values.

So the "rev" the search was for appears on two levels: the **flow** reorganises (corner
vortices grow, centre moves) and the **discretisation** strains (method agreement
loosens). Both are measured, both are consistent with the literature.

---

## 4. What is and isn't claimed

**Verified:**
- Both PSPG and Taylor–Hood are close to Ghia at `Re=400` and `Re=1000` at one resolution
  each (TH more accurate per node); refinement at these Re is not tested here.
- The two methods still agree (mutual RMS `0.006`/`0.016`), with agreement loosening as Re
  rises — the discretisation's strain made visible.
- The physical rev is captured **qualitatively**: with the stream function, the secondary
  eddies sit at Ghia's locations, counter-rotate and grow with Re, and the primary centre
  migrates to the Ghia-published location.

**Not claimed:**
- **Still steady laminar.** No unsteady/turbulent transition (that is `~Re 8000`, and
  would need time integration and far finer meshes). The "rev" here is *within* the steady
  regime plus the numerical strain — not a Hopf bifurcation to time-dependence.
- `Re=1000` is **not** fully mesh-converged at these resolutions (RMS `~0.03` for TH);
  the trend is convergent but a finer mesh / higher order would tighten it.
- The eddy **strengths** are not mesh-converged (BR `74%`, BL `57%`, primary `86%` of
  Ghia at `n=128`, `Re=1000`, all still rising with `n`); only PSPG is used for the
  stream-function analysis.
- 2-D only; no new mathematics — Ghia, PSPG, Taylor–Hood are classical. The contribution
  is the dual-element validation and the explicit, measured account of the rev.

---

## 5. Files

- Theory: this note, on `08`/`08a`.
- Verification: `src/verification3d/cavity_highre_verify.py` (checks A–C; imports both
  element solvers).
- Next rung: a true time-dependent run toward unsteadiness (much higher Re, or a
  deliberately destabilised setup), and/or mesh-convergence at `Re=1000` to tighten the
  numbers — both larger efforts (see `07a`).
