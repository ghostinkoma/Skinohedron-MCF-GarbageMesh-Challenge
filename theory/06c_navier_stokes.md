# 06c · V3 Stage C — Full Incompressible Navier–Stokes (the nonlinear coupling)

**Status:** theory + **verification (PASS)**. Backed by
`src/verification3d/fluid_navier_stokes_verify.py` (checks A–E). Stage C couples
everything: viscosity (`03`), the pressure projection (`04`), and the convective
operator (`06b`) — now with the velocity advecting **itself**, `C(u)u`. This is the
first genuinely **nonlinear** stage, so the yardstick is the nonlinear-CFD suite
(exact flows, manufactured/conservation/convergence) rather than a single `10⁻¹⁵`.

**Headline finding (carefully stated).** `‖Bu‖ ~ 10⁻¹⁶` by itself is an *algebraic
triviality*: any projector `P = I − M⁻¹Bᵀ(BM⁻¹Bᵀ)⁻¹B` satisfies `BP = 0` by
construction, so machine-precision discrete divergence proves nothing on its own —
**this is not the claim**, and the earlier framing ("cancels the classical
constraint") over-stated it. The real, tested claims are two: (i) with a *consistent*
`(B,Bᵀ)` pair from the one per-tet `G`, the velocity stays divergence-free without an
inconsistency residual (the earlier `3.8e-4` was a non-transpose bug, not inf–sup);
and (ii) **the equal-order P1 pressure is genuinely checkerboard-unstable, but that
instability is confined to the pressure and does *not* leak into the velocity** —
shown directly (not asserted) by PSPG-invariance: stabilising the pressure
(correlation `0.0→0.94`) moves the velocity by only `~10⁻⁴`. Built on `06a`/`06b`/`04`.

---

## 0. The equations and the scheme

Incompressible Navier–Stokes:

$$
\frac{\partial u}{\partial t} + (u\cdot\nabla)u = -\frac{1}{\rho}\nabla p + \nu\nabla^2 u,
\qquad \nabla\cdot u = 0.
$$

**Chorin projection step** `uⁿ → uⁿ⁺¹`, component-wise (`u=(u_x,u_y,u_z)`):

1. **Predictor (nonlinear advection + viscosity):**
   $$
   M\,u^\ast = M\,u^{n} + \Delta t\big(-\,C(u^{n})\,u^{n} \;-\; \nu\,K\,u^{n}\big),
   $$
   `C(uⁿ)` rebuilt each step from the current velocity (the nonlinearity, quadratic
   in `u`); `C` skew-symmetric (`06b`), `−K` the verified vector Laplacian.

2. **Projection (incompressibility) — the consistent form:**
   $$
   \big(B\,M^{-1}B^{\mathsf T}\big)\,p = B\,u^\ast,
   \qquad
   u^{n+1} = u^\ast - M^{-1}B^{\mathsf T} p .
   $$

So Stage C = Stage B's predictor (now self-advecting) + an **exact** projection.
The only new physics is `u` in `C(u)` being the unknown.

---

## 1. The consistent projection — why it is machine-exact

The velocity is nodal P1 (natural for advection/viscosity). Define one discrete
divergence operator `B` mapping nodal velocity → nodal scalar, assembled from the
per-tet gradient `G`:
$$
(B\,u)_i = \sum_{\text{tet}\ni i} \tfrac{\text{vol}}{4}\,(\nabla\!\cdot u)|_{\text{tet}},
\qquad (\nabla\!\cdot u)|_{\text{tet}} = \sum_{a}\sum_{c} G_{c,a}\,u_{\,\text{node}_a,c}.
$$
Use the **exact transpose** `Bᵀ` for the pressure gradient in the correction. Then
the projection is an algebraic projector onto `ker B`:
$$
B\,u^{n+1} = B\,u^\ast - (B M^{-1}B^{\mathsf T})(B M^{-1}B^{\mathsf T})^{-1} B u^\ast = 0
$$
**exactly**, to solver precision — independent of inf–sup. *Verified:* random field
`‖Bu‖` `2.5e-1 → 4.8e-16`, idempotent `P²=P` `1.7e-13` (check A).

**The earlier `3.8e-4` was an artefact, not a wall.** A first attempt paired a
per-tet divergence with a *vertex-averaged* gradient — the two were **not
transposes**, so the projector was inexact. Diagnostic: a per-tet-closed projection
already gave `‖Du‖=9.9e-15` while the inconsistent round-trip gave `1.8e-1`. The
fix is consistency, not stabilisation: one `G`, one `B`, its exact `Bᵀ`.

---

## 2. The nonlinear ingredient

`C(u)u` is quadratic. In skew form `C_skew = ½(C−Cᵀ)`, advection conserves discrete
kinetic energy exactly in the inviscid, divergence-free limit (`uᵀC_skew u = 0`), so
energy leaves only through viscosity — the physically correct behaviour and the
basis of stability. `C(uⁿ)` is reassembled each step.

---

## 3. Verification results (`fluid_navier_stokes_verify.py`, all PASS)

- **A. Consistent projection.** `‖Bu‖ → 4.8e-16`, `P²=P` `1.7e-13`. Incompressibility
  is **machine-precision on equal-order P1**.

- **B. Taylor–Green vortex (nonlinear).** Exact 2D solution `u=(\sin kx\cos ky,
  -\cos kx\sin ky)e^{-2\nu k^2 t}`, energy `E=E_0 e^{-4\nu k^2 t}`. Measured decay
  `3.059` vs analytic `3.158` (**3.1%**, FE accuracy); **max `‖Bu‖ = 4.3e-17` over
  the entire nonlinear run** — machine-precision incompressibility *maintained*.

- **C. Inviscid energy conservation (`ν=0`).** KE drift **0.000%** (`0.5→0.5`),
  `‖Bu‖ = 3.3e-17`. Skew advection + exact projection conserve energy.

- **D. Convergence.** Decay-rate error `n=8: 4.9% → n=10: 3.1% → n=12: 2.2%`,
  falling with refinement; `‖Bu‖` machine (`≤4e-16`) at every `n`.

- **E. Pressure honesty + velocity-quality (the inf–sup investigation).** Backed by
  `fluid_ns_velocity_quality_verify.py`. Three measured facts:
  1. *Pressure is genuinely bad.* The unstabilised projection pressure correlates
     only `0.0–0.2` with the exact Taylor–Green pressure — classic checkerboard.
  2. *Velocity is genuinely clean.* High-frequency velocity energy (FFT, above ⅓
     Nyquist) is `~10⁻⁹` for Taylor–Green and `~10⁻⁴` for an *asymmetric* multi-mode
     flow, **decreasing under refinement** (`9e-4→4e-4→1e-5` for n=8,10,12) — a
     convergent discretisation effect, not a blow-up.
  3. *The decisive test — PSPG-invariance.* Stabilising the pressure raises its
     correlation `0.0→0.94` while changing the velocity by only `~10⁻⁴`. So the
     checkerboard lives in the pressure and does **not** contaminate the velocity:
     the velocity is pressure-decoupled. (The `~10⁻⁴` asymmetric high-freq content is
     therefore physical cascade, since it is PSPG-invariant.)

---

## 4. What is and isn't claimed

**Verified (Stage C):**
- Full incompressible Navier–Stokes is the self-advecting Chorin projection on the
  verified operators; only `C(u)u` is new.
- The consistent `(B,Bᵀ)` projection keeps the velocity divergence-free without the
  earlier inconsistency residual; Taylor–Green decay `4νk²` is reproduced to FE
  accuracy (3.1%, convergent); inviscid energy is conserved (0.000% drift).
- **The velocity is pressure-decoupled**: PSPG stabilisation fixes the pressure
  (correlation `0.0→0.94`) while moving the velocity by only `~10⁻⁴`, so the
  checkerboard does not leak into the velocity (shown, not asserted).

**Not claimed / conceded:**
- **`‖Bu‖~10⁻¹⁶` is not, by itself, a result.** It is an algebraic property of any
  projection; the earlier "cancels the classical constraint" framing over-stated it
  and is retracted. The substantive content is the *velocity quality* (above), not
  the divergence norm.
- **The unstabilised pressure is not trustworthy.** Equal-order P1/P1 is inf–sup
  (LBB) unstable; the pressure is checkerboard (correlation `~0.0–0.2`). A usable
  pressure needs PSPG/compatible elements — available, but the baseline pressure
  should be treated as unreliable. So "we solve correct incompressible NS" holds for
  the **velocity**, not yet for the pressure.
- **Battery is incomplete.** Tested on periodic Taylor–Green and an asymmetric
  multi-mode flow only. **Untested** (and required before any broad claim):
  wall-bounded lid-driven cavity, flow past a cylinder, transition to turbulence,
  and long-time enstrophy behaviour. Periodic smooth flows can flatter a scheme.
- **No new mathematics.** Chorin projection, the algebraic `ker B` projector, PSPG,
  Taylor–Green are classical; the contribution is the self-consistent realisation on
  this project's one `G`, with the pressure/velocity split measured honestly.
- **Bounded reach:** incompressible, Newtonian, *moderate* Reynolds number on the
  `n=8`–`12` structured cube; high-Re turbulence, compressibility, free surfaces, and
  a faithful pressure are out of scope.

---

## 5. Files

- Theory: this document, on `06`/`06a`/`06b`/`04`.
- **Verification:** `src/verification3d/fluid_navier_stokes_verify.py` (checks A–E)
  and `src/verification3d/fluid_ns_velocity_quality_verify.py` (high-freq + PSPG-invariance,
  answering the inf–sup critique).
- Next (optional): a Navier–Stokes viewer (vortex decay) reusing the fluid-viewer
  rendering; a pressure-recovery study; the V3 synthesis update.
