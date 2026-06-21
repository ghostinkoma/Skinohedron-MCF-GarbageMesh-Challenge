# 06a · V3 Stage A — Stokes Flow (linear), the plumbing before the nonlinearity

**Status:** theory + **verification (PASS)** for Stage A of V3. Backed by
`src/verification3d/fluid_stokes_verify.py` (checks A–E, measured numbers inline). Stokes flow drops the
nonlinear advection term, leaving a **fully linear** problem — so it is verifiable
against exact profiles at near-machine precision, exactly as the linear era was.
Its job is to fix the *plumbing* — the vector momentum / pressure-projection /
boundary coupling — before advection (Stage B) and the full nonlinear system
(Stage C). Built on `04` (projection) and `03` (vector Laplacian); see `06` for the
V3 model and roadmap.

---

## 1. The Stokes equations

Drop `(u·∇)u` from incompressible Navier–Stokes:

$$
\frac{\partial u}{\partial t} = -\frac{1}{\rho}\nabla p + \nu\,\nabla^2 u,
\qquad \nabla\cdot u = 0.
$$

Steady Stokes is `0 = −(1/ρ)∇p + ν∇²u`, `∇·u = 0`. Everything here is linear; the
only operators are the verified vector Laplacian `−K` (component-wise) and the
verified projection of `04`.

---

## 2. Channel geometry and mixed boundaries (the new element)

Stage A introduces a boundary type the linear era mostly avoided: a **wall**
(Dirichlet velocity). Use the unit cube as a plane channel:

- `y ∈ [0,1]` is the **wall-normal** direction. Walls at `y=0` and `y=1` carry
  Dirichlet velocity: **no-slip** `u=0`, or a **driven** wall `u=U`.
- `x` (streamwise) and `z` (spanwise) are **periodic** (the flow is homogeneous
  along them).

So the boundary is *mixed*: periodic in `x,z`, Dirichlet in `y`. The discrete
solve fixes the wall DOFs and condenses them out of `K` (a standard Dirichlet
lift), the same reduction used implicitly throughout the FE assembly.

For Couette and Poiseuille the flow is **unidirectional**, `u = (u_x(y),0,0)`, so
`∇·u = ∂u_x/∂x = 0` holds automatically and the pressure decouples — which is why
Stage A is clean. (The projection is still exercised separately, §5 check E.)

---

## 3. Exact solution 1 — Couette flow (wall-driven)

Top wall moves at `U`, bottom fixed, no pressure gradient. Steady Stokes reduces to
`ν u_x'' = 0` with `u_x(0)=0`, `u_x(1)=U`:

$$
\boxed{\,u_x(y) = U\,y\,}\qquad\text{(linear profile).}
$$

A linear field lies exactly in the P1 space, so the discrete solver must recover it
to **machine precision** — the velocity analogue of the `T=k₂/(k₁+k₂)` material
test. *Probe:* `max|u_x − Uy| = 9.1\times10^{-15}`.

---

## 4. Exact solution 2 — Poiseuille flow (pressure-driven)

A constant streamwise pressure gradient `G = (1/ρ)\,dp/dx` drives flow between fixed
walls. Steady Stokes reduces to `ν u_x'' = G` with `u_x(0)=u_x(1)=0`:

$$
\boxed{\,u_x(y) = \frac{-G}{2\nu}\,y(1-y)\,}\qquad\text{(parabolic profile),}
$$

with centreline speed `u_max = -G/(8\nu)` at `y=½`. For a constant source the P1 FE
solution of the 1-D Poisson problem is **nodally exact**, so this too is recovered
to machine precision. *Probe:* `max|u_x − \text{parabola}| = 1.9\times10^{-15}`,
centreline `0.12500` (exact `0.12500`).

---

## 5. Verification design — both steady and transient

Stage A is verified two independent ways, each a PASS/FAIL. **All five pass**:
Couette steady `9.1e-15`, Poiseuille steady `1.9e-15` (centreline `0.12500`),
Couette transient decay rate `9.762` vs `νπ²=9.870` (`1.1%`), Poiseuille transient
converged `4.7e-9`, incompressibility `‖Du‖=3.0e-14` / `P²=P` `2.1e-13`.

**A. Couette steady.** Solve the wall-driven steady problem; assert
`max|u_x − Uy|` is machine precision. *Target:* `< 10⁻¹²`.

**B. Poiseuille steady.** Solve the pressure-driven steady problem; assert
`max|u_x − (-G/2ν)y(1-y)|` is machine precision and centreline `= -G/(8ν)`.
*Target:* `< 10⁻¹²`.

**C. Couette transient → steady.** Start from rest and march
`∂u/∂t = ν∇²u` with the driven-wall BC; the deviation from steady decays through
the diffusion modes, the slowest being `sin(\pi y)` with rate
$$
\sigma_1 = \nu\,\pi^2/L^2 \quad (L=1).
$$
*Target:* the transient converges to the steady profile, and the late-time decay
rate of `‖u − u_steady‖` matches `νπ²` to FE accuracy.

**D. Poiseuille transient → steady.** Same march with a constant body force `−G`
from rest; assert convergence to the parabolic steady state and the same slowest-mode
decay approach. *Target:* converges to the §4 parabola.

**E. Incompressibility under projection (the coupling anchor).** Because Couette /
Poiseuille are unidirectional, they do not exercise the projection. So take a
*non-trivial* tentative velocity `u*` (a constructed 3-D field with `∇·u* ≠ 0`),
apply the verified projection `K p = (ρ/Δt) D u*`, `u = u* − (Δt/ρ) G p`, and assert
`‖D u‖ →` machine precision and idempotency `P²=P` — confirming the vector solver's
pressure coupling is the verified projection of `04`, now wired to velocity.
*Target:* `‖D u‖ < 10⁻¹⁰`, `P²=P` to `~10⁻¹⁴`.

Together: two exact steady profiles at machine precision, two transient approaches at
the analytic decay rate, and the incompressibility/projection anchor — the full
plumbing checked before any nonlinearity.

---

## 6. What is and isn't claimed

**Verified (Stage A, machine precision / FE accuracy):**
- The wall-bounded linear Stokes solver reproduces Couette (linear) and Poiseuille
  (parabolic) to machine precision, both as steady solves and as the limit of a
  transient march at the analytic decay rate.
- The vector momentum / pressure-projection coupling is the verified projection of
  `04`, keeping `∇·u = 0` to machine precision.

**Not claimed:**
- No nonlinearity is present (that is Stage B's `C(u)` and Stage C's coupling); the
  near-machine precision here is *because* Stokes is linear and will not survive
  into the nonlinear stages, where the yardstick becomes convergence/MMS (`06` §4).
- No new mathematics: plane Couette/Poiseuille and the Dirichlet lift are classical;
  the contribution is the verified realisation on this project's `L`.
- Runs first on the clean cube (no `01e` anisotropy in the wall-normal 1-D profile,
  no `5b` sign-flips); deformed-shape channels are later.

---

## 7. Files

- Theory: this document, on `06`/`04`/`03`.
- **Verification (this stage):** `src/verification3d/fluid_stokes_verify.py` —
  checks A–E above. Reuses `fem_laplacian` (`K, M`), the per-tet gradient `G`, and
  the divergence/projection from `pressure_field_verify.py`.
- Next: Stage B (`06b`, scalar advection–diffusion, the operator `C(u)`), then the
  fluid viewer.
