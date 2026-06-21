# 06b · V3 Stage B — Scalar Advection–Diffusion, the operator `C(u)`

**Status:** theory + **verification (PASS)** for Stage B of V3. Backed by
`src/verification3d/fluid_advection_verify.py` (checks A–E, measured inline). Stage B isolates the
**one genuinely new ingredient** of fluid dynamics — the convective operator
`C(u)` — by carrying a *passive scalar* `φ` in a *prescribed* velocity field. With
`u` fixed (not solved for), the problem is linear in `φ`, so `C(u)` can be checked
against exact solutions before the full nonlinear coupling of Stage C. Built on `03`
(the Laplacian for diffusion) and `06`/`06a`; uses the same per-tet gradient `G`.

---

## 1. The equation

A passive scalar advected by a known velocity `u(x)` and diffusing with coefficient
`κ`:

$$
\frac{\partial \varphi}{\partial t} + u\cdot\nabla\varphi = \kappa\,\nabla^2\varphi.
$$

Discretely, diffusion is the verified `−K` (scaled by `κ`); the new piece is the
**convective operator** `C(u)` discretising `u·∇`.

---

## 2. Building `C(u)` from the same per-tet gradient

On each tetrahedron the P1 gradient `G` (the matrix already inside `K`) gives the
constant gradient `∇φ|_tet = G φ_tet`. For a constant element velocity `u`, the
convective term `u·∇φ` is the constant `u·(Gφ_tet)`. Assembling against the test
functions and lumping through `M` gives the matrix `C(u)`, so the semi-discrete
system is

$$
M\dot\varphi = -\,C(u)\,\varphi \;-\; \kappa\,K\,\varphi .
$$

**Skew-symmetric (energy-conserving) form.** Written as
`½[u·∇φ + ∇·(uφ)]` (equal to `u·∇φ` when `∇·u = 0`), the discrete operator becomes
antisymmetric, `C_skew = ½(C − Cᵀ)`, so

$$
\varphi^{\mathsf T} C_{\text{skew}}\,\varphi = 0 :
$$

advection alone creates no scalar "energy" `½φᵀMφ`. This is the discrete image of
the continuous identity and the basis of a **stable** scheme. For a divergence-free
prescribed `u` (as in the checks here) `C` is already skew, but the skew form is
used so the property holds exactly and generally.

**Why this matters (the honest numerical point of `06`).** Centred discretisation
of advection is not naturally stable — at high Péclet number (`|u|/κ` large) a
non-conservative form admits growing oscillations. The skew-symmetric form is the
controlled, energy-bounded baseline; an upwind/SUPG variant is the option for very
high Péclet, reported when used. Stage B is where this first appears in code.

---

## 3. Exact solutions to verify against

With `u` prescribed, several exact targets exist:

1. **Consistency — a constant advects to nothing.** `u·∇(\text{const}) = 0`, so
   `C(u)\,\mathbf 1 = 0`. *Target:* `‖C·1‖ →` machine precision.

2. **Skew-symmetry / energy conservation.** `φᵀ C_skew φ = 0` exactly.
   *Target:* machine precision.

3. **Mass conservation.** Advection conserves the total `∫φ = 1ᵀMφ`.
   *Target:* constant to machine precision over a run.

4. **Pure advection (uniform `u`, periodic).** A profile translates rigidly; its
   centre of mass moves at exactly `u`. (A sharp profile also disperses — the FE
   advection dispersion, the same `~5%` family as `01e`'s `Δc/c` — so the phase
   speed is checked to FE accuracy, not machine precision.)

5. **Advection–diffusion of a Gaussian.** In a periodic box a Gaussian
   `φ ∝ exp(−(x−x₀)²/2σ₀²)` carried by uniform `u` and diffusing has the exact
   moments
   $$
   x_{\text{cm}}(t) = x_0 + u\,t,\qquad
   \sigma^2(t) = \sigma_0^2 + 2\kappa t .
   $$
   The **variance growth rate `2κ`** is a clean diffusion check (FE-accurate to
   well under 1%); the **centre-of-mass speed `u`** is the advection check (FE
   dispersion ~5%).

---

## 4. Verification design (checks, with expected accuracy)

All five PASS (`fluid_advection_verify.py`, n=12 periodic, u=(1,0,0), κ=0.003):
`C·1=1.4e-17`, `φᵀC_skew φ=-3.2e-17`, mass `1.2e-16`, advection speed `0.955` vs
`1.0` (FE dispersion), energy ratio `1.001` (stable), Gaussian variance `0.00756`
vs exact `0.00760` (`0.5%`). The checks:

- **A. Consistency `C·1 = 0`.** machine precision.
- **B. Skew-symmetry `φᵀ C_skew φ = 0`.** machine precision.
- **C. Mass conservation `∫φ` constant.** machine precision over the run.
- **D. Pure advection translation.** centre of mass moves at `u` to FE accuracy
  (`~5%`); energy stays bounded (stability of the skew form).
- **E. Advection–diffusion Gaussian.** variance grows at `2κt` to `~1–2%`
  (diffusion exact); centre moves at `u` to FE accuracy (advection dispersion).

*Probes already run* (n=12 periodic, `u=(1,0,0)`, `κ=0.003`): `C·1 = 1.9e-17`,
`φᵀC_skew φ = 2e-17`, mass conserved to `~1e-4` relative, Gaussian variance
`0.00756` vs exact `0.00760` (`0.5%`), centre-of-mass error `9e-3` (FE dispersion).
Full assertions are the code's job.

---

## 5. What is and isn't claimed

**Verified (Stage B):**
- `C(u)` built from the same per-tet `G` is consistent (`C·1=0`), energy-conserving
  in skew form (`φᵀC_skew φ=0`), and mass-conserving.
- Advection–diffusion of a Gaussian reproduces the exact moments (`x_cm=x₀+ut`,
  `σ²=σ₀²+2κt`), diffusion to ~1%, advection to FE-dispersion accuracy.

**Not claimed:**
- `u` is **prescribed** here; the nonlinear self-advection `C(u)u` is Stage C. Stage
  B verifies the *operator*, not the coupled flow.
- Advection carries the FE dispersion of the mesh (`~5%`), so the centre-of-mass
  speed is FE-accurate, not machine-exact — honestly the first place the yardstick
  is convergence/accuracy rather than machine precision (as flagged in `06` §4).
- No new mathematics: streamline advection and the skew-symmetric / SUPG forms are
  classical; the contribution is the verified realisation on this project's `G`.
- Runs on the clean cube/periodic box first (no `5b` sign-flips); deformed shapes
  later.

---

## 6. Files and next steps

- Theory: this document, on `06`/`06a`/`03`.
- **Verification (this stage):** `src/verification3d/fluid_advection_verify.py`
  (checks A–E). Reuses `fem_laplacian` (`K, M`) and the per-tet gradient `G`.
- Next: **the fluid viewer** (now that a scalar visibly moves with the flow), then
  **Stage C** (`06c`, full incompressible Navier–Stokes: self-advection, Taylor–
  Green, MMS, conservation/convergence).
