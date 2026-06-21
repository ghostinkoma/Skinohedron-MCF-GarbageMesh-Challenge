# 06 · Fluid Dynamics — Incompressible Navier–Stokes on the Verified Operator (V3, theory)

**Status:** theory document (code is the next stage). It opens V3 by stating the
fluid model the project will build, assembling it almost entirely from pieces
*already verified* (the operator `L = M⁻¹K`, the gradient/divergence `D∘grad = K`,
and the idempotent pressure projection `P² = P` of `04`), isolating the **one new
ingredient** that is genuinely hard — nonlinear advection — and, in keeping with
空論を重ねない, naming the exact-solution (or manufactured-solution) check for every
stage. Nothing here is asserted as numerically established yet; this is the design
the verification must meet.

Built on `03`/`03b` (the operator and its dynamics) and especially `04` (pressure,
projection). See the prologue `00a` for why the ground is ready.

---

## 0. The model

Incompressible Navier–Stokes for a velocity field `u(x,t)` and pressure `p(x,t)`:

$$
\frac{\partial u}{\partial t} + (u\cdot\nabla)u
   = -\frac{1}{\rho}\nabla p + \nu\,\nabla^2 u,
\qquad \nabla\cdot u = 0.
$$

Read term by term against what the notebook already has:

| term | meaning | discrete operator | status |
|---|---|---|---|
| `ν ∇²u` | viscous diffusion | component-wise scalar `L` (vector Laplacian) | **verified** (`03`) |
| `(1/ρ)∇p` | pressure gradient | the per-tet gradient `G` (the one inside `K`) | **verified** (`04`, `D∘grad=K`) |
| `∇·u = 0` | incompressibility | the projection `P` (`∇²p = (ρ/Δt)∇·u*`) | **verified** (`04`, `P²=P` to `2.9e-14`) |
| `(u·∇)u` | **nonlinear advection** | a new per-tet convective operator `C(u)` | **the one new, hard piece** |

So V3 is *not* built from scratch. Three of the four terms are already the verified
operator in different clothing; only advection is new. The whole of V3 is, in
essence, **the verified Chorin projection of `04` with an advection–diffusion step
in front of it.**

---

## 1. The one new ingredient — and why it changes the rules

Every physics so far (heat, wave, pressure, even acoustics) has been **linear**.
That linearity is exactly what let the project verify against exact solutions at
`10⁻¹¹`–`10⁻¹⁶`: linear fields lie in the P1 space, so the discrete operator
reproduces them to machine precision (the `T=k₂/(k₁+k₂)` and `D∘grad=K` checks).

The advection term `(u·∇)u` is **nonlinear** (the field transports itself). Two
consequences must be faced honestly:

1. **No general exact solution.** The clean "machine-precision-vs-exact" paradigm
   does not apply to a generic nonlinear flow. The verification method itself must
   change (§4) — to manufactured solutions, special exact flows, conservation
   invariants, and convergence order. This is the standard, honest CFD toolkit.

2. **Discrete advection is not naturally stable.** Centred differencing of
   `u·∇` is energy-neutral but admits unphysical oscillations at high Péclet/Reynolds
   number; a faithful scheme needs either an upwind/streamline-stabilised form
   (SUPG-type) or a skew-symmetric ("energy-conserving") discretisation. The model
   must *name* the stabilisation, not hide it. The chosen baseline is a
   skew-symmetric convective operator (it conserves discrete kinetic energy in the
   inviscid, divergence-free limit), with an upwind option for high-Re robustness.

Honesty up front: V3 is the first chapter where "verified" will mean *converges to
the right answer at the right rate*, not *matches to machine precision* — because
that is what nonlinear physics permits. The discipline is unchanged; only the yardstick is.

---

## 2. The discrete scheme (Chorin projection on the verified `L`)

Velocity is stored per-vertex, per component (`u = (u_x,u_y,u_z)`, each a P1 field).
A single time step from `uⁿ` to `uⁿ⁺¹`, reusing the verified operators:

**(1) Advection–diffusion (predictor).** Form a tentative velocity that ignores
incompressibility,
$$
M\,u^\ast = M\,u^{n} + \Delta t\big(-\,C(u^{n})\,u^{n} \;-\; \nu\,K\,u^{n}\big),
$$
where `C(u)` is the per-tet convective operator (§3) and `−K` is the (already
verified) vector Laplacian applied component-wise. (Viscosity may be treated
implicitly, `(M+\nu\Delta t\,K)u^\ast=\dots`, for stability at large `ν`.)

**(2) Pressure projection (corrector).** `u*` is not divergence-free; correct it
exactly as in `04`:
$$
K\,p = \frac{\rho}{\Delta t}\,D\,u^\ast,
\qquad
u^{n+1} = u^\ast - \frac{\Delta t}{\rho}\,G\,p .
$$
This is the verified projection `P`: it removes the gradient part and leaves
`∇·uⁿ⁺¹ = 0` to solver precision, and it is idempotent (`P²=P`, measured `2.9e-14`
in `04`). Pressure is the Lagrange multiplier enforcing incompressibility, solved
instantaneously each step (pure-Neumann `K`, pinned as in `04`).

The structure is worth stating plainly: **steps (2) are Step 2, unchanged.** V3
adds only step (1), and within it only `C(u)` is new.

---

## 3. The convective operator `C(u)`

On each tetrahedron the P1 gradient `G` (the same matrix that builds `K`) gives a
constant field gradient, so for a scalar `φ` carried by the flow,
`(u·∇)φ |_tet = (u_tet · G φ_tet)`, with `u_tet` the element-averaged velocity.
Assembled against the test functions and lumped through `M`, this yields a
convective matrix `C(u)`; for the vector momentum each component is advected the
same way.

Two properties define the baseline choice:

- **Skew-symmetry / energy conservation.** Written in the skew-symmetric form
  `½[(u·∇)φ + ∇·(uφ)]` (equal to `u·∇φ` when `∇·u=0`), the discrete `C(u)` is
  antisymmetric, so `φᵀ M C(u) φ = 0`: advection alone neither creates nor destroys
  discrete kinetic energy. This is the discrete analogue of the continuous identity
  and the basis of a stable inviscid scheme.

- **Consistency.** `C(u)` is built from the *same* per-tet `G` as `K` and `D`, so it
  inherits the project's one geometric source of truth; it is the same operator
  family, now bilinear in `(u, ·)`.

An upwind/SUPG variant adds a streamline-aligned stabilisation for high Reynolds
number; it is a controlled modification of `C(u)`, reported (not hidden) when used.

---

## 4. Verification design — the honest yardstick for nonlinear flow

Because nonlinear advection has no general exact solution, V3 verifies against the
four standard, rigorous tools. Each becomes a PASS/FAIL when code lands.

1. **Incompressibility (machine precision, still).** The constraint `∇·uⁿ⁺¹ = 0`
   is *linear* and enforced by the verified projection, so it must hold to solver
   precision every step. *Target:* `‖D uⁿ⁺¹‖ → ` machine precision (like `04`'s
   `P²=P`). This anchor is unchanged from the linear era.

2. **Special exact solutions.**
   - **Taylor–Green vortex** (2D, unsteady): `u = (\cos x \sin y,\,-\sin x \cos y)\,e^{-2\nu t}`
     is an exact decaying solution; the solver must reproduce the `e^{-2\nu t}`
     energy decay and the velocity field.
   - **Poiseuille / Couette** (steady): exact linear/parabolic profiles between
     plates under a pressure gradient or wall motion; the steady state must match.
   - **Beltrami flow** (3D): `∇×u ∥ u` gives closed-form unsteady solutions for a
     3D check.

3. **Method of Manufactured Solutions (MMS).** Choose a smooth divergence-free
   `u_exact`, insert it into Navier–Stokes to compute the forcing `f` that makes it
   exact, then verify the solver driven by `f` recovers `u_exact`. This manufactures
   an exact target for an otherwise unsolvable nonlinear problem — the cleanest way
   to keep 空論を重ねない alive past linearity.

4. **Conservation & convergence order.** In the inviscid, periodic limit the
   skew-symmetric `C(u)` must conserve discrete kinetic energy (drift → 0 with
   `Δt`); with viscosity the energy must decay at the analytic rate. Under mesh and
   time-step refinement the error against (2)/(3) must fall at the scheme's design
   order. *Target:* energy drift bounded; observed order matches design.

These four, together, are what "verified" means for a nonlinear solver: not one
machine-precision number, but a constraint held exactly, agreement with every exact
flow available, a manufactured exact target, and the right convergence rate.

---

## 5. The staged roadmap (each stage verified before the next)

V3 advances the way the whole notebook has — one isolable difficulty at a time:

- **Stage A — Stokes flow (linear).** Drop advection entirely:
  `∂u/∂t = −(1/ρ)∇p + ν∇²u`, `∇·u = 0`. This is *fully linear* — verifiable at
  near-machine precision against Poiseuille/Couette (steady) and the Stokes-limit
  Taylor–Green decay. It exercises the vector Laplacian + projection coupling with
  **no** nonlinearity, confirming the plumbing before the hard part.

- **Stage B — scalar advection–diffusion (prescribed `u`).** Solve
  `∂φ/∂t + u·∇φ = κ∇²φ` for a *passive scalar* carried by a *known* velocity field.
  This isolates and verifies the new operator `C(u)` alone: pure advection of a
  profile by uniform `u` (translation), and advection–diffusion with an exact
  Gaussian-spreading solution. The nonlinearity of NS is absent (u is prescribed),
  so `C(u)` can be checked cleanly.

- **Stage C — full incompressible Navier–Stokes.** Couple it all: momentum with
  self-advection `C(u)u`, viscosity, and projection. Verify against Taylor–Green,
  MMS, and conservation/convergence (§4). This is the genuinely nonlinear target.

Each stage has a complete, runnable verification before the next begins; a failing
stage stops the line, exactly as the scattered node and the sphere obstruction did.

---

## 6. What is and isn't claimed

**Design claims (to verify next stage):**
- Incompressible Navier–Stokes is the verified Chorin projection of `04` with an
  advection–diffusion predictor; three of its four terms are the already-verified
  operator, and only advection `C(u)` is new.
- Incompressibility remains a machine-precision check (linear, via the projection);
  the nonlinear parts are verified by exact flows, manufactured solutions, and
  conservation/convergence.

**Not claimed:**
- No numbers are established in this document — it is the model and the verification
  plan only.
- **No new mathematics.** Chorin projection, the cotangent Laplacian as vector
  viscosity, skew-symmetric convection, and the MMS/Taylor–Green verification suite
  are classical CFD. The contribution is, as ever, a *single self-consistent,
  verified* realisation on this project's own `L`.
- **Bounded reach.** This targets *incompressible, Newtonian, moderate-Reynolds*
  flow on the `n=8`-class structured mesh. Turbulence modelling (LES/RANS),
  compressible shocks, free surfaces, and non-Newtonian rheology are out of scope.
- The residual `~5%` mesh anisotropy (`01e`) and, on deformed shapes, the
  sign-flipped cotangents (`5b`) are inherited; on the cube they are absent, so the
  baseline V3 checks run on the clean cube first.

---

## 7. Files and next steps

- Theory: this document, on `03`/`03b`/`04` and the operator `L`.
- Reused, already verified: `fem_laplacian` → `K, M` (`src/ksf3d/fem3d.py`), the
  per-tet gradient `G`, the divergence `D` and projection from
  `pressure_field_verify.py`.
- **Next stage (code):** `src/verification3d/fluid_*_verify.py` per stage —
  `fluid_stokes_verify.py` (Stage A), `fluid_advection_verify.py` (Stage B, the
  `C(u)` operator), `fluid_navier_stokes_verify.py` (Stage C: Taylor–Green, MMS,
  conservation); then a fluid viewer reusing the unified/solid rendering to show the
  velocity field (streamlines / speed heatmap) on cube, sphere, and torus.
