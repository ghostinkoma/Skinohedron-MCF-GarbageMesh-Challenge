# 04 · The Pressure Field: One Poisson Solve for Liquids and Gases (Step 2, theory)

**Status:** theory + **verification (PASS)**. The model below is now backed by
`src/verification3d/pressure_field_verify.py`, which asserts each of the five
exact-solution targets of §4 to machine precision. In keeping with 空論を重ねない,
the measured residuals are recorded inline.

It builds directly on the verified operator `L = M⁻¹K` of `03`/`03b`.

---

## 0. The one-line model

Pressure is the scalar field whose Laplacian balances whatever pushes the fluid:

```
continuum:   ∇²p = ∇·f
discrete:    K p = D f            (same K as heat/wave; D = discrete divergence source)
```

`f` is "what pushes the fluid" (force or momentum per unit volume). **Three
physically different situations are this same solve**, differing only in what `f`
is and which boundary condition applies:

| situation | what `f` is | boundary | gives |
|-----------|-------------|----------|-------|
| **hydrostatic** | `f = ρ g` (gravity) | Neumann `∂p/∂n = ρ g·n` (walls) | `p = ρ g·x` (= ρgh) |
| **projection (incompressible)** | `f = (ρ/Δt) u*` (tentative momentum) | Neumann `∂p/∂n = 0` (closed) | the `p` that makes `∇·u = 0` |
| **acoustic (gas)** | linearised EOS coupling | Neumann/Dirichlet | sound waves, `c² = dp/dρ` |

The third is exactly Step 1's scalar wave on `L` — so **Step 2 contains Step 1**
as the compressible (gas) limit. That is the high-abstraction unification: one
operator, one Poisson/​wave structure, liquids and gases selected by parameters.

---

## 1. Why pressure is a Poisson solve on `L`

Momentum (incompressible, inviscid) is `ρ (∂u/∂t + u·∇u) = −∇p + f`. Two ways the
Laplacian of pressure appears, both landing on the same `K`:

**(a) Equilibrium / hydrostatic.** At rest (`u = 0`, steady), momentum reduces to
`∇p = f`. Taking the divergence: `∇²p = ∇·f`. With `f = ρg` and constant `ρ,g`,
the right side is `0`, so `p` is **harmonic** (`∇²p = 0`) with a Neumann wall
condition `∂p/∂n = ρ g·n`. The solution is the linear field `p = ρ g·x + const`.

**(b) Projection / incompressibility.** Given a tentative velocity `u*` that is
*not* divergence-free, the Chorin projection corrects it,
`u = u* − (Δt/ρ) ∇p`, and demanding `∇·u = 0` gives
`∇²p = (ρ/Δt) ∇·u*`. Pressure is the **Lagrange multiplier** that enforces
incompressibility; it is whatever harmonic-plus-source field removes the
divergence of the tentative flow.

In FE/DEC form the stiffness factors as `K = Gᵀ W G` (G = per-tetrahedron P1
gradient, W = volume weights), so the *same* `K` that ran heat and waves is the
pressure Laplacian, and `G`/`Gᵀ` supply the gradient that turns `p` into a force
and the divergence that turns a velocity field into a source. No new operator —
the gradient/divergence are already inside the assembled `K`.

---

## 2. Liquids vs gases — selected by parameters

The same solve, two equations of state:

**Liquid (incompressible).** `∇·u = 0` is an exact constraint; pressure carries no
equation of state and is solved instantaneously each step (the projection above).
Hydrostatics is **linear**: `p = ρ g h` (depth `h`).

**Gas (compressible, barotropic).** Pressure follows an EOS `p = p(ρ)` (e.g.
`p = c²ρ`, isothermal). Two consequences:

- Hydrostatics becomes **exponential**, not linear: balancing `dp/dh = −ρ g`
  with `p = c²ρ` gives the isothermal atmosphere `p(h) = p₀ exp(−g h / c²)`.
- Small disturbances obey the **wave equation** with `c² = dp/dρ` — i.e. Step 1's
  `M p̈ = −K p` with the sound speed set by the EOS. Acoustics is the gas limit of
  this same pressure field.

So the **liquid ↔ gas** switch is a parameter (compressibility / EOS), and it
continuously connects the static Poisson solve (liquid, instantaneous) to the
hyperbolic wave (gas, acoustic) already verified in Step 1.

---

## 3. Parameters (the "what to simulate" knobs)

Designed from the theory above, to be exposed (as in the unified viewer) so one
instrument covers the whole family:

| parameter | meaning | liquid | gas |
|-----------|---------|--------|-----|
| `ρ` density | inertia / weight | high, ~constant | low, varies with `p` |
| `g` gravity (vector) | body force; `g=0` removes hydrostatic gradient | on | on |
| `κ` compressibility / `c²` | EOS stiffness; `κ→0` = incompressible | `κ≈0` | finite `c²` |
| boundary | closed tank (Neumann) vs open surface (Dirichlet `p=p_atm`) | either | either |
| source `u*` | tentative velocity whose divergence drives projection | for flow | for flow |

Limiting cases the knobs must reproduce:
`g≠0, u*=0, κ=0` → hydrostatic `p=ρgh`;
`g=0, u*≠0, κ=0` → pure incompressible projection;
`g≠0, κ>0` (gas) → exponential atmosphere;
small `κ>0`, oscillating source → acoustic waves (= Step 1).

---

## 4. Verification (all five PASS to machine precision)

Implemented in `src/verification3d/pressure_field_verify.py`. First, a consistency
check confirms the gradient/divergence are the same operator as `K`:
`‖D(grad p) − K p‖ / ‖K p‖ = 3.4e-16`. Then:

1. **Hydrostatic, closed tank (liquid).** Solve `K p = D(ρg)` with Neumann walls.
   *Target:* `grad p = ρg` everywhere, `p` affine. *Measured:*
   `max|grad p − ρg| = 3.6e-11`, non-affine residual `6.5e-12`. The linear field
   `p = ρ g·x` lies in P1, so it is recovered to machine precision (like the
   `T = k₂/(k₁+k₂)` interface test of `02`). **PASS.**

2. **Incompressible projection (Helmholtz).** Project per-tet fields with
   `∇²p = ∇·w`, `u = w − ∇p`. *Measured:* a pure-gradient field projects to
   `‖Pw‖/‖w‖ = 3.8e-15` (fully removed); a random field keeps a non-trivial
   div-free part (`0.96`) with divergence `‖D(Pw)‖ = 1.4e-14`; the projector is
   idempotent, `‖Pu − u‖/‖u‖ = 2.9e-14` (`P² = P`). **PASS.**

3. **Gas hydrostatic (exponential atmosphere).** With EOS `p = c²ρ`. *Measured:*
   the gas profile `p = p₀ exp(−g h / c²)` departs from the linear liquid law for
   small `c²` (gap `8.81` at `c²=1`) and recovers it as `c²→∞` (gap `4.8e-11` at
   `c²=1e6`); column ratio `p(H)/p₀` exact. **PASS.**

4. **Acoustic limit = Step 1.** Leapfrog `M p̈ = −c² K p`. *Target:* oscillation
   `ω = c·√λ₁`. *Measured:* `ω` matches `c·√λ₁` to relative `0.000` at `c = 1, 2`
   — the gas acoustic limit reproduces Step 1's scalar wave, frequency ∝ c.
   **PASS.**

5. **Solvability / null space.** Pure-Neumann `K` is singular with a constant null
   space. *Measured:* `‖K·1‖ = 6.2e-15` (constant is the null vector) and the
   source is compatible, `|1·b| = 7.2e-15` (orthogonal to constants ⇒ solvable);
   pressure unique after pinning one reference dof. **PASS.**

---

## 5. What is and isn't claimed (Step 2 honesty)

**Verified (§4, machine precision):**
- Hydrostatic, incompressible projection, and gas acoustics are three regimes of
  one Poisson/​wave solve on the already-verified `L`.
- Liquid/gas is a parameter (compressibility / EOS); the gas acoustic limit is
  exactly Step 1's scalar wave (`ω = c√λ`).

**Not claimed:**
- No new mathematics: the pressure Poisson equation, Helmholtz/Chorin projection,
  isothermal atmosphere and linear acoustics are classical. The contribution is a
  *verified, self-consistent* realisation on this project's own `L`, with the
  gradient/divergence shown to be the very operator inside `K` (`D∘grad = K`).
- Viscosity, advection (`u·∇u`), and full Navier–Stokes are **out of scope** for
  Step 2; this is the pressure/incompressibility layer only. Momentum transport
  is later (Step 3).
- The same residual ~5% mesh anisotropy from `01e`/`03` is inherited; the linear
  hydrostatic test is near-exact regardless (affine fields lie in P1).

---

## 6. Files and next steps

- Theory: this document, on `03`/`03b` and the operator `L`.
- Operator already available: `src/ksf3d/fem3d.py` (`fem_laplacian` → `K`, `M`);
  the P1 gradient/divergence `G`/`Gᵀ` implicit in `K = GᵀWG` will be exposed for
  the source term `D f`.
- **Next stage (code):** `src/verification3d/pressure_field_verify.py` asserting
  the five targets of §4; then a dedicated pressure viewer
  (`viewer/viewer_pressure.html`) with the §3 parameters (ρ, g, compressibility,
  boundary, source) so liquid/gas/atmosphere/acoustic are all selectable in one
  instrument.
