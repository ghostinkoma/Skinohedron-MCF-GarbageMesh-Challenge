# 06c · V3 Stage C — Full Incompressible Navier–Stokes (the nonlinear coupling)

**Status:** theory + **verification (PASS)**. Backed by
`src/verification3d/fluid_navier_stokes_verify.py` (checks A–E). Stage C couples
everything: viscosity (`03`), the pressure projection (`04`), and the convective
operator (`06b`) — now with the velocity advecting **itself**, `C(u)u`. This is the
first genuinely **nonlinear** stage, so the yardstick is the nonlinear-CFD suite
(exact flows, manufactured/conservation/convergence) rather than a single `10⁻¹⁵`.

**Headline finding.** The apparent classical limitation — *"equal-order P1/P1 is not
inf–sup stable, so the projection cannot make the velocity divergence-free to machine
precision"* — is **cancelled** here by a careful, consistent accumulation: building
the discrete divergence `B` and its **exact transpose `Bᵀ`** from the *same* per-tet
gradient `G`. The velocity is then divergence-free to machine precision (`‖Bu‖ ~
10⁻¹⁶`), **maintained through the nonlinear evolution and at every resolution**.
What inf–sup genuinely costs is confined to the **pressure** field, stated openly
below. Built on `06a`/`06b`/`04`.

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

- **E. Pressure honesty (the inf–sup residual).** The projection pressure correlates
  only `0.24` with the exact Taylor–Green pressure: equal-order P1 admits spurious
  pressure modes (in the near-null space of `Bᵀ`) that **do not pollute the
  velocity**. So the velocity is divergence-free and energy-correct, but the
  **pressure** is not faithful. This is where inf–sup actually lives.

---

## 4. What is and isn't claimed

**Verified (Stage C):**
- Full incompressible Navier–Stokes is the self-advecting Chorin projection on the
  verified operators; only `C(u)u` is new.
- With the consistent `(B,Bᵀ)` projection, the **velocity is divergence-free to
  machine precision** through the nonlinear flow and at every resolution; Taylor–Green
  decay `4νk²` is reproduced to FE accuracy (3.1%, convergent); inviscid energy is
  conserved (0.000% drift).
- The apparent inf–sup obstruction to machine-precision incompressibility is an
  artefact of inconsistent operators; a consistent transpose pair removes it.

**Not claimed:**
- **The pressure is faithful.** Equal-order P1/P1 admits spurious pressure modes
  (correlation `~0.24` with exact); inf–sup's genuine cost is here, in the pressure,
  not in the velocity divergence. A faithful pressure needs a compatible/stabilised
  pair or a pressure-recovery post-process — a separate, later option. The
  **velocity** is the trustworthy output.
- **No new mathematics.** Chorin projection, skew convection, the algebraic
  `ker B` projector, Taylor–Green are classical; the contribution is the verified,
  self-consistent realisation on this project's one `G` — and the explicit
  demonstration that consistency, not stabilisation, buys machine-precision
  incompressibility here.
- **Bounded reach:** incompressible, Newtonian, *moderate* Reynolds number on the
  `n=8`–`12` structured cube. High-Re turbulence (needs upwind/SUPG), compressibility,
  free surfaces, and faithful pressure are out of scope.

---

## 5. Files

- Theory: this document, on `06`/`06a`/`06b`/`04`.
- **Verification:** `src/verification3d/fluid_navier_stokes_verify.py` (checks A–E).
- Next (optional): a Navier–Stokes viewer (vortex decay) reusing the fluid-viewer
  rendering; a pressure-recovery study; the V3 synthesis update.
