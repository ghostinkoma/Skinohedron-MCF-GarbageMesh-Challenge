# 06d · Inertial Mass and the `e^{−λt}` Structure — Decomposing the Taylor–Green Error

> **Reconstructed (2026-11).** The original `06d` was cited by the theory index, `06f`,
> `06z`, `07` and `09` but never committed. This note was rewritten from those citations;
> every number is produced by `src/verification3d/fluid_ns_error_decomposition_verify.py`
> (also reconstructed). One cited phrase could not be reproduced as stated and is
> narrowed below ("all integrators agree").
>
> **Reading guide.** "Mass and viscosity as one coin" (the original framing) is an
> interpretation [N]; the measured content is the decomposition in §2–§4. See
> `theory/README.md`, *How to read the claims*.
>
> **Known theory.** Time-stepping error `O(dt)` that vanishes as `dt → 0` is
> Lax–Richtmyer convergence of a stable, consistent scheme; the remaining `O(h²)`
> eigenvalue error is Babuška–Osborn.

**Status:** **verified** (reconstructed), `src/verification3d/fluid_ns_error_decomposition_verify.py`.
Builds on `06c` (Navier–Stokes on the periodic Kuhn cube).

---

## 1. The question

`06c` measures the Taylor–Green kinetic-energy decay rate within `~3%` of the continuum
`4νk²` (`k = 2π`, `ν = 0.02`). Where does the `~3%` come from — the nonlinear advection,
the time stepping, or the space discretisation?

## 2. Advection is projected out

For Taylor–Green the continuum advection `u·∇u` is an exact gradient, so incompressibility
removes it. Discretely (skew-symmetric advection on the n=8 Box), the fraction of the
advection term that the projection removes is **0.9999**. The velocity therefore decays as
a viscous mode, at the rate set by the discrete operator:

$$\text{rate} = 2\nu\,\lambda_h, \qquad \lambda_h = \frac{u^\top K u}{u^\top M u}.$$

## 3. The time error is `O(dt)` and removable

Same Box (n=8), same end time, forward-Euler step with projection:

| `dt` | measured rate | `2νλ_h` | difference |
|---|---|---|---|
| `3.0e-3` | 3.0060 | 2.9992 | `6.77e-3` |
| `1.5e-3` | 3.0026 | 2.9992 | `3.38e-3` |
| `7.5e-4` | 3.0009 | 2.9992 | `1.70e-3` |

The difference halves with `dt` (ratio 2.00): first-order, and it vanishes. The rate is
`2νλ_h`. *(The original phrase "all integrators agree" is narrowed: only the forward-Euler
scheme of `06c` at three step sizes is tested here.)*

## 4. What remains is spatial, `O(h²)`

| n | `λ_h` | `2k²` | error |
|---|---|---|---|
| 8  | 74.981 | 78.957 | `−5.04%` |
| 10 | 76.393 | 78.957 | `−3.25%` |
| 12 | 77.169 | 78.957 | `−2.26%` |
| 16 | 77.947 | 78.957 | `−1.28%` |

Order `h^1.98`. A sinusoid is not in P1, so its Rayleigh quotient carries the classical
`O(h²)` error — the mirror of why the linear (Couette) and nodally superconvergent
(Poiseuille) profiles of `06a` were machine-precise. These are the same numbers as the
stiffness error of `06f` (`uᵀKu` vs `2k²E`).

## 5. What is and isn't claimed

**Verified (reconstructed script):** advection 99.99% gradient for Taylor–Green; time
error `O(dt)` → 0 with rate → `2νλ_h`; spatial error `λ_h` vs `2k²` `5.0 → 1.3%`, `O(h²)`.

**Not claimed:** agreement across different time integrators (only one scheme tested);
any result beyond the Taylor–Green mode; the "one coin" reading as physics.

## 6. Files

- Verification: `src/verification3d/fluid_ns_error_decomposition_verify.py` (checks A–C; ~30 s).
- Uses the Box / `run_TG` of `src/verification3d/fluid_navier_stokes_verify.py`.
