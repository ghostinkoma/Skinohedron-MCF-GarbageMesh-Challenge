# 06h · Temperature as the Glue — Mediating Coulomb and Mass

**Status:** theory + **verification (PASS)**, `src/verification3d/thermal_glue_verify.py`.
A stronger role for temperature than the parallel fusion of `06g`. There, temperature
was an additive domain (its energy summed with the kinetic). Here it is the
**mediator**: heating honey lowers its viscosity, so temperature does not act alone —
it acts **through** the inter-molecular (Coulomb) bonds to change the mass transport
(flow). Temperature is the **glue** between the **Coulomb** domain (stiffness /
viscosity, `06e`) and the **mass** domain (momentum / flow). Builds on `02` (the
two-material interface), `06d`–`06g`.

---

## 1. The mediation is multiplicative — temperature lives *inside* `K`

`06g`'s fusion was **additive** (total energy = kinetic + potential). The glue is
**multiplicative**: the viscosity enters the stiffness operator as a
**temperature-dependent weight**,
$$
K_\nu \;=\; \int \nu\big(T(x)\big)\,\nabla\varphi\cdot\nabla\psi \,dV,
$$
so temperature sits *inside* the Coulomb/stiffness operator, scaling it pointwise.
Verified: with constant `T`, `K_ν = ν·K` exactly (`‖K_ν−νK‖ = 1.3e-15`, check D). So
this is genuinely a modulation of the one operator, not a new one.

The temperature → viscosity law is **Arrhenius**, `ν(T) = ν₀ exp(E_a/T)`, where the
activation energy `E_a` is the inter-molecular (Coulomb, `06e`) bond barrier. The
*same* `ΔT` produces very different viscosity changes by material (check A):

| liquid | `E_a` | `ν(300K)/ν(330K)` |
|---|---|---|
| water (weak H-bonds) | 1800 | **1.7×** |
| honey (strong H-bond network) | 6000 | **6.2×** |

Temperature acting *through* `E_a` (Coulomb) is the correlation the insight names —
not a parallel domain, a mediating one.

---

## 2. The mediated interface IS the verified `02` structure — to machine precision

The decisive check: a flow with a cold (viscous) layer and a hot (thin) layer. Steady
variable-viscosity shear gives stress continuity `μ ∂u/∂y = const`, so the interface
velocity is the **harmonic-mean** form
$$
\boxed{\,u_{\text{interface}} = U\,\frac{\mu_2}{\mu_1+\mu_2}\,}
$$
— **identical in structure to `02`'s two-material temperature interface
`T = k₂/(k₁+k₂)`.** Verified to machine precision:

- two viscosities `μ₁=4, μ₂=1`: interface `u = 0.200000` vs exact `0.200000`,
  error `7.5e-16`; slope ratio `0.2500 = μ₂/μ₁` (stress continuous) — check B.
- **Arrhenius-set** viscosities (`μ(300K)=6.16`, `μ(330K)=1.00`): interface
  `u = 0.139652` vs exact, error `1.1e-15` — check C.

So the temperature field, through `ν(T)`, sets the flow exactly, and the governing
formula is the **mechanical twin of the already-verified heat interface**. Honey
thinning when heated is precisely this mediation — and it lands on a structure the
notebook proved long ago, to machine precision.

---

## 3. Why this matters — the glue completes the domain map

With this, the three domains relate as the insight predicted:

- **mass** `M` — inertia; exact in P1 (`06f`).
- **Coulomb** `K` — stiffness/viscosity; `O(h²)` for smooth modes, → machine precision
  by enrichment (`06f`); reproduces the electrostatic `1/r` (`06e`).
- **temperature** `T` — both a **parallel** domain carrying its own conserved energy
  (`06g`) **and** the **glue** that multiplicatively couples Coulomb to mass:
  `ν(T)` inside `K_ν`, with the decay rate `2ν(T)·λ_h`, `λ_h = K/M`. Temperature scales
  the Coulomb/mass ratio.

So temperature is not a fourth, separate axis — it is the **mediator** that makes the
mass and Coulomb axes interact, exactly as a glue should. The verified harmonic-mean
interface is the concrete fingerprint of that mediation.

*(A structural resemblance to a field whose local value sets the effective coupling of
other fields is noted only as a direction; it is deliberately **not** formalised here,
per the project's rule against writing ahead of verification.)*

---

## 4. What is and isn't claimed

**Verified:**
- Temperature enters multiplicatively inside the stiffness (`K_ν`, constant-`T` →
  `ν·K` exact); Arrhenius `ν(T)` makes the same `ΔT` act differently by Coulomb
  barrier `E_a` (honey vs water).
- The temperature-mediated two-viscosity interface is `u = U μ₂/(μ₁+μ₂)` to machine
  precision — the mechanical twin of `02`'s `T = k₂/(k₁+k₂)`, including Arrhenius-set
  viscosities.

**Not claimed:**
- No first-principles `ν(T)` from molecular Coulomb dynamics — Arrhenius is used as the
  (classical, empirical) bridge; deriving `E_a` itself needs the MD engine flagged in
  `06e`, out of scope.
- The interface result is steady, piecewise-linear (in P1, hence machine-precise);
  time-dependent thermo-viscous flow with feedback (temperature advected by the flow it
  controls) is a larger nonlinear study, named not done.
- The resemblance to a coupling-setting field is **not** developed here (deferred by
  intent); no physical claim beyond the verified mediation is made.
- No new mathematics: Arrhenius viscosity, harmonic-mean interfaces, and variable-
  coefficient FE are classical; the contribution is showing temperature mediates
  Coulomb↔mass on this project's operator, landing on the verified `02` structure.

---

## 5. Files

- Theory: this note, on `02`/`06e`/`06f`/`06g`.
- Verification: `src/verification3d/thermal_glue_verify.py` (checks A–D).
- Implied direction: time-dependent thermo-viscous flow (temperature advected by, and
  controlling, the flow) — the nonlinear feedback the glue enables.
