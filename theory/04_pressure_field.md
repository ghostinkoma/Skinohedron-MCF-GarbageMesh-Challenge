# 04 · The Pressure Field: One Poisson Solve for Liquids and Gases (Step 2, theory)

**Status:** theory document only (code is the next stage). It opens Step 2 by
stating the single abstract model the project will build pressure on, and — in
keeping with 空論を重ねない — it names, for every claim, the **exact solution
that will verify it** once code exists. Nothing here is asserted as numerically
established yet; this is the design the verification must meet.

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

## 4. Verification design (what each claim will be checked against)

Code is the next stage; these are the **exact solutions** each regime will be
asserted against, so the model cannot drift into 空論:

1. **Hydrostatic, closed tank (liquid).** Solve `K p = D(ρg)` with Neumann walls.
   *Exact target:* `p = ρ g·x + const`, i.e. pressure linear in depth, gradient
   exactly `ρg`. Check recovered `p` is affine to ~machine precision (linear
   fields are in the P1 space, so this should be near-exact, like the
   `T = k₂/(k₁+k₂)` interface test of `02`).

2. **Incompressible projection (Helmholtz).** Build a known field
   `w = ∇φ + (∇×A)` on the mesh; project with `∇²p = ∇·w`, set `u = w − ∇p`.
   *Exact target:* `u` is divergence-free and equals the curl part; the removed
   part equals `∇φ`. Check `‖∇·u‖` drops to solver tolerance and the projector is
   idempotent (`P² = P`).

3. **Gas hydrostatic (exponential atmosphere).** With EOS `p=c²ρ` and gravity.
   *Exact target:* `p(h) = p₀ exp(−g h / c²)`; check the profile and that
   `c²→∞` (incompressible limit) recovers the linear liquid case.

4. **Acoustic limit = Step 1.** Linearise the gas about rest; the disturbance
   must satisfy `M p̈ = −K p` with `c² = dp/dρ`.
   *Exact target:* the dispersion/decay already verified in `03`/`03b`
   (oscillation `cos(√λ t)`), now reached as the gas limit — a consistency check
   tying Step 2 back to Step 1.

5. **Solvability / null space.** Pure-Neumann (closed tank) `K` is singular with a
   constant null space (pressure defined up to a constant). *Target:* the source
   `D f` must be orthogonal to the constant (compatibility, total flux balances),
   and pressure is pinned by one reference value. Check the compatibility
   condition and a unique solution after pinning.

---

## 5. What is and isn't claimed (Step 2 honesty)

**Design claims (to be verified next stage):**
- Hydrostatic, incompressible projection, and gas acoustics are three regimes of
  one Poisson/​wave solve on the already-verified `L`.
- Liquid/gas is a parameter (compressibility / EOS); the gas acoustic limit is
  exactly Step 1's scalar wave.

**Not claimed yet:**
- No numbers are verified in this document — it is the theory and the verification
  plan only. Each target in §4 becomes a PASS/FAIL check when code lands.
- Viscosity, advection (`u·∇u`), and full Navier–Stokes are **out of scope** for
  Step 2; this is the pressure/incompressibility layer only. Momentum transport
  is later (Step 3).
- The same residual ~5% mesh anisotropy from `01e`/`03` is inherited; the linear
  hydrostatic test is expected near-exact regardless (affine fields lie in P1).

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
