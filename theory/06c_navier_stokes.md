# 06c · V3 Stage C — Full Incompressible Navier–Stokes (the nonlinear coupling)

**Status:** theory document (code is the next stage). Stage C couples everything:
viscosity (`03`), the pressure projection (`04`), and the convective operator
(`06b`) — now with the velocity advecting **itself**, `C(u)u`. This is the first
genuinely **nonlinear** stage, so — as flagged in `06` §4 — the yardstick changes
from "machine precision vs exact" to the standard nonlinear-CFD suite (exact
flows, manufactured solutions, conservation, convergence order). Grounding probes
(reported inline) confirm the direction; full verification is the next stage.

Built on `06a` (Stokes), `06b` (advection), `04` (projection).

---

## 0. The equations and the scheme

Incompressible Navier–Stokes:

$$
\frac{\partial u}{\partial t} + (u\cdot\nabla)u = -\frac{1}{\rho}\nabla p + \nu\nabla^2 u,
\qquad \nabla\cdot u = 0.
$$

**Chorin projection step** `uⁿ → uⁿ⁺¹`, reusing the verified operators
component-wise (`u=(u_x,u_y,u_z)`):

1. **Predictor (nonlinear advection + viscosity):**
   $$
   M\,u^\ast = M\,u^{n} + \Delta t\big(-\,C(u^{n})\,u^{n} \;-\; \nu\,K\,u^{n}\big).
   $$
   The convective matrix `C(uⁿ)` is rebuilt each step from the *current* velocity —
   this is the nonlinearity (`C(u)u`, quadratic in `u`). `C` is the skew-symmetric
   operator of `06b`; `−K` is the verified vector Laplacian.

2. **Projection (incompressibility):** solve the same Poisson as `04`,
   $$
   K\,p = \frac{\rho}{\Delta t}\,D\,u^\ast,
   \qquad
   u^{n+1} = u^\ast - \frac{\Delta t}{\rho}\,G\,p .
   $$

So Stage C = **Stage B's predictor (now self-advecting) + Stage A/`04`'s
projection.** The only new thing is that `u` in `C(u)` is the unknown, not a
prescribed field. Everything else is already verified.

---

## 1. Why the yardstick changes (and stays honest)

`(u·∇)u` is nonlinear, so there is no general exact solution and the
machine-precision paradigm cannot apply to a generic flow. Two honest consequences:

- **Verification moves to the nonlinear-CFD suite** (§3): special exact flows
  (Taylor–Green), the Method of Manufactured Solutions, conservation invariants,
  and convergence order. The discipline of 空論を重ねない is unchanged — every claim
  still meets a defined target — but the target is "right answer at the right rate,"
  not "10⁻¹⁵."

- **Equal-order P1/P1 is not inf–sup stable.** With velocity and pressure both
  piecewise-linear (collocated), the discrete projection does **not** drive
  divergence to machine precision; it reduces it to a small residual that shrinks
  with refinement. This is a classical limitation, stated plainly here rather than
  hidden. Exact discrete incompressibility would need a compatible space pair or a
  PSPG-type stabilisation — out of scope for the baseline. *Probe:* `‖D u‖ ≈ 4e-4`
  on `n=10` (small, convergent — not machine).

---

## 2. The nonlinear ingredient, concretely

`C(u)u` is quadratic: the operator `C(u)` depends on `u`, then acts on `u`. Two
properties carry over from `06b` and matter for stability:

- **Skew-symmetry → energy conservation.** With the skew form `C_skew = ½(C−Cᵀ)`,
  advection conserves discrete kinetic energy `½uᵀMu` exactly in the inviscid,
  divergence-free limit (`uᵀ C_skew u = 0`). So energy can only leave through
  viscosity — the physically correct behaviour and the basis of stability.
- **Rebuilt each step.** `C(uⁿ)` is reassembled from the current velocity every
  step (explicit treatment of advection); viscosity may be implicit for large `ν`.

---

## 3. Verification design — the nonlinear suite

Each becomes PASS/FAIL when code lands (`fluid_navier_stokes_verify.py`):

1. **Taylor–Green vortex (the primary exact flow).** In 2D the decaying vortex
   $$
   u = \big(\sin kx\,\cos ky,\;-\cos kx\,\sin ky\big)\,e^{-2\nu k^2 t}
   $$
   is an **exact solution of the full nonlinear equations**, with kinetic energy
   decaying as `E(t) = E_0\,e^{-4\nu k^2 t}`. *Target:* the solver reproduces the
   decay rate `4νk²`. *Probe (n=10, ν=0.02):* measured `3.074` vs analytic `3.158`,
   **2.7%** — the nonlinear advection+viscosity+projection coupling gives the right
   physics to FE accuracy.

2. **Method of Manufactured Solutions (MMS).** Pick a smooth divergence-free
   `u_exact`, insert it into Navier–Stokes to get the forcing `f` that makes it
   exact, drive the solver with `f`, and verify it recovers `u_exact`. This
   manufactures an exact target for an otherwise unsolvable nonlinear problem — the
   cleanest way to keep 空論を重ねない past linearity. *Target:* error → 0 with
   refinement at the design order.

3. **Energy conservation (inviscid limit).** With `ν=0` and the skew-symmetric `C`,
   discrete kinetic energy must be conserved (drift → 0 as `Δt → 0`). *Target:*
   bounded, vanishing energy drift.

4. **Incompressibility (honest).** The projection reduces `‖D u‖` to a small
   residual (probe `~4e-4`), decreasing under refinement — *not* machine precision
   on equal-order P1 (the inf–sup caveat of §1). *Target:* residual small and
   convergent; reported, not overclaimed.

5. **Convergence order.** Under mesh/time refinement the error against (1)/(2) falls
   at the scheme's design rate. *Target:* observed order matches design.

Together these are what "verified" means for a nonlinear solver: an exact flow
reproduced, a manufactured exact target met, energy conserved in the inviscid
limit, incompressibility convergent, and the right order — not a single 10⁻¹⁵.

---

## 4. What is and isn't claimed

**Design claims (Stage C, to verify in code):**
- Full incompressible Navier–Stokes is the self-advecting Chorin projection on the
  verified operators; only `C(u)u`'s nonlinearity is new.
- The Taylor–Green decay rate `4νk²` is reproduced to FE accuracy (probe 2.7%);
  inviscid energy is conserved; MMS recovers the manufactured solution.

**Not claimed:**
- **Not machine-precision incompressibility.** Equal-order P1/P1 is not inf–sup
  stable; the projection reduces divergence to a small convergent residual, not to
  `10⁻¹⁵`. This is stated openly (it is the honest cost of the collocated baseline);
  a stabilised/compatible pair is a separate, later option.
- **No new mathematics.** Chorin projection, skew-symmetric convection, Taylor–Green
  and MMS are classical CFD; the contribution is the verified realisation on this
  project's `L`.
- **Bounded reach:** incompressible, Newtonian, *moderate* Reynolds number on the
  `n=8`–`10` structured cube. High-Re turbulence (needs upwind/SUPG or a turbulence
  model), compressibility, and free surfaces are out of scope.
- Runs on the clean periodic cube first (no `01e` anisotropy issues beyond the known
  ~5% FE level, no `5b` sign-flips); shaped-domain nonlinear flow is later.

---

## 5. Files and next steps

- Theory: this document, on `06`/`06a`/`06b`/`04`.
- **Verification (next):** `src/verification3d/fluid_navier_stokes_verify.py` —
  Taylor–Green decay, MMS, inviscid energy conservation, incompressibility residual,
  convergence order.
- Then: optionally a Navier–Stokes viewer (vortex decay / lid-driven cavity) reusing
  the fluid-viewer rendering; and the V3 synthesis update.
