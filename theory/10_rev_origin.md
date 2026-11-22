# 10 · Why Does the "Rev" Exist? — First Separation (Revised)

> **Reconstructed and revised (2026-11).** The original `10` was cited by the theory
> index, `11`, `12`, `13` and the precision viewer but never committed. It was rewritten
> from those citations, and its central result is **revised**: the original "physical rev
> converges (Cauchy difference `0.042 → 0.0002`) → field-intrinsic" came from a `u_min`
> proxy that measures the primary vortex's return flow, not the corner eddy (issue #3,
> A3; see `08b`). Every number below is produced by `src/verification3d/rev_origin_verify.py`
> (also reconstructed).
>
> **Reading guide.** "Rev" and "field-intrinsic" are labels [N]; the measured content is
> in §2–§4. See `theory/README.md`, *How to read the claims*.
>
> **Known theory.** Corner eddies: Moffatt (1964); benchmark values: Ghia, Ghia & Shin
> (1982), Tables I and V; element stability: Babuška–Brezzi, PSPG, Taylor–Hood.

**Status:** **verified** (reconstructed, conclusion revised), `src/verification3d/rev_origin_verify.py`.
Takes up the open question of `09` §5.

---

## 1. The question and the four candidates

`08b` found a Reynolds-driven change ("rev") in the lid-driven cavity: corner eddies
appear and grow, the primary vortex migrates, and the two elements (PSPG, Taylor–Hood)
agree less tightly at `Re = 1000`. `09` listed four candidate sources:

1. **approximation error** (finite `h`, element order);
2. **a missing domain** (e.g. Coulomb / EHD forces);
3. **a field-intrinsic structure** of the continuous Navier–Stokes solution;
4. **chemistry** (reaction terms).

## 2. The numerical rev shrinks — approximation error (1)

At `Re = 1000`, against Ghia's centreline:

| measure | coarse | medium | fine |
|---|---|---|---|
| PSPG vs Ghia (n = 64 / 96 / 128) | 0.0555 | 0.0404 | 0.0320 |
| Taylor–Hood vs Ghia (n = 24 / 32 / 40) | 0.0338 | 0.0265 | 0.0211 |
| PSPG vs Taylor–Hood (pairs (96,32) / (112,36) / (128,40)) | 0.0161 | 0.0142 | 0.0130 |

Every measure decreases under refinement: the discretisation's disagreement with the
benchmark and between the two elements is **approximation error**. (The precision viewer's
earlier middle value `0.0140` came from an unrecorded pair; `(112, 36)` gives `0.0142`.)

## 3. The physical rev: real structure, strength not yet converged (revised)

Measured with the stream function `ψ` (`Kψ = ∫φ(∂v/∂x − ∂u/∂y)`, `ψ = 0` on the walls),
bottom-right (BR) eddy at `Re = 1000`, PSPG:

| n | BR `ψ` | location | % of Ghia (`1.751e-3` at (0.8594, 0.1094)) |
|---|---|---|---|
| 64  | `9.59e-4` | (0.875, 0.109) | 55% |
| 96  | `1.18e-3` | (0.875, 0.115) | 67% |
| 128 | `1.30e-3` | (0.867, 0.109) | 74% |

- **Location: resolved.** The eddy sits within `0.02` of Ghia's location at every
  resolution.
- **Strength: not converged.** It is still rising (increments `2.2e-4`, `1.2e-4`), from
  below, toward Ghia. P1/PSPG under-resolves the eddy at `n ≤ 128`.

**Revision.** The original conclusion — "the physical rev *converges* to a fixed value,
so it is field-intrinsic (3)" — does **not** hold: its series (`0.0989 → 0.1413 → 0.1414`)
was the primary return flow at the edge of a measurement window. What the data supports:
the eddies are real flow structure (present where Ghia finds them, counter-rotating,
growing with Re), which is consistent with (3); but their **strength** is still dominated
by approximation error (1) at these meshes, so the two are **not yet separated** for the
eddy strength.

## 3b. Follow-up: the eddy strength does converge — slowly (added 2026-11)

`src/verification3d/corner_eddy_convergence_verify.py` refines further with **two**
elements: PSPG to `n = 192`, and Taylor–Hood (P2/P1) to `n = 64`, with `ψ` from the
Taylor–Hood velocity on its P2 node grid. BR eddy at `Re = 1000`:

| element | n | BR `ψ` | % of Ghia |
|---|---|---|---|
| PSPG | 128 / 160 / 192 | `1.304e-3` / `1.383e-3` / `1.442e-3` | 74 / 79 / 82 |
| Taylor–Hood | 32 / 40 / 48 / 56 / 64 | `1.267e-3` / `1.351e-3` / `1.410e-3` / `1.455e-3` / `1.490e-3` | 72 / 77 / 80 / 83 / 85 |

- **Monotone convergence.** Both sequences rise with shrinking increments.
- **Same limit from two elements.** Fitting `ψ(h) = ψ∞ − C hᵖ` to the last four levels
  gives `ψ∞ ≈ 1.92e-3` for both PSPG and Taylor–Hood (they differ by well under 5%), within
  ~10% of Ghia's `1.75e-3`. Two different discretisations agreeing on the limit is the
  evidence that the eddy strength is a property of the continuous solution (candidate 3),
  not of the mesh.
- **Slow.** The observed order is `p ≈ 0.6`, well below the rate for smooth solutions.
  This is consistent with — but does not prove — pollution by the lid-corner singularity
  (the velocity is discontinuous at the two top corners). It also means the extrapolated
  value is pinned only to about 10%; it does not reproduce Ghia's digits (Ghia's own value
  comes from a 129 × 129 grid).

**Updated separation.** Numerical rev: approximation error, shrinking (§2). Physical rev:
real flow structure whose location, Re-dependence **and** (now) limiting strength are
element-independent — field-intrinsic (3), resolved to ~10% — with a slowly decaying
approximation error on top.

## 4. The eddy grows with Re (the advection weight)

Scaling the advection by `α` at fixed `ν` is the same as replacing `Re` by `αRe`, so the
`α`-dependence is the `Re`-dependence. At `n = 96`:

| Re | BR `ψ` | note |
|---|---|---|
| 100  | not resolved | Ghia: `1.25e-5`, below this mesh |
| 400  | `4.16e-4` at (0.885, 0.115) | Ghia `6.42e-4` |
| 1000 | `1.18e-3` at (0.875, 0.115) | Ghia `1.75e-3` |

Growth `400 → 1000`: `×2.8` (Ghia `×2.7`). The eddy strengthens as advection strengthens,
as Moffatt's corner analysis and Ghia's data describe — within the **existing**
Navier–Stokes operator.

## 5. Candidates (2) and (4)

Ghia's benchmark is pure incompressible Navier–Stokes, and every gap to it shrinks under
refinement (§2) or moves toward it (§3). No additional domain (2) or chemistry (4) is
needed to account for the cavity. (Both remain relevant elsewhere: `12`.)

## 6. What is and isn't claimed

**Verified (reconstructed script):** the numerical rev shrinks (all three measures); the
corner eddy is at Ghia's location and grows with Re (`×2.8` vs `×2.7`); its strength rises
`55 → 67 → 74%` of Ghia at n = 64 → 128.

**Verified (follow-up, §3b):** the eddy strength converges monotonically for PSPG
(to n = 192) and Taylor–Hood (to n = 64); both extrapolate to the same limit
`≈ 1.9e-3` (within ~10% of Ghia); observed order `p ≈ 0.6`.

**Not claimed:** that the original Cauchy-difference series measured the corner eddy (it
did not — the convergence now shown is a different, stream-function result); the eddy
strength to better than ~10%; that the slow order is caused by the corner singularity
(consistent, not proven).

**Next (optional):** a corner-singularity treatment (e.g. a regularised lid) or local
refinement at the top corners, to test whether the order recovers.

## 7. Files

- Verification: `src/verification3d/rev_origin_verify.py` (checks A–C; ~5 min);
  `src/verification3d/corner_eddy_convergence_verify.py` (§3b; ~10 min).
- Uses `src/verification3d/cavity_highre_verify.py` (stream function, Ghia Table V
  constants) and the PSPG / Taylor–Hood solvers of `08`, `08a`.
