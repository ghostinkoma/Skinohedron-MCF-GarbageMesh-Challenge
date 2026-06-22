# 06e · Coulomb and Viscosity — the Other Term, and an Honest Multi-Scale Boundary

**Status:** part **verified** (`coulomb_operator_verify.py`), part **conceptual map
(explicitly not yet verified)**. This note takes up an insight about the *origin* of
the viscosity `ν` that `06c`/`06d` treat as a given constant: at the molecular level,
**honey and water have nearly equal mass density but vastly different viscosity**, so
viscosity cannot come from mass — it comes from **inter-molecular forces**, which are
at root **Coulomb** (hydrogen bonds, dipole–dipole, van der Waals are all
electrostatic). The note verifies the checkable link and draws the rest as an honest
map, marking exactly where this notebook's verified ground ends.

---

## 1. The verified link — Coulomb is already the operator `K`

Electrostatics is `−ε∇²φ = ρ`, i.e. `K φ = M ρ` — **the same `K`** as the pressure
Poisson of `04` and the Laplacian of `03`. So the Coulomb interaction is not a new
physics to bolt on; it is the operator the notebook already trusts, read with charge
as the source. Verified (`coulomb_operator_verify.py`):

- **A. Coulomb 1/r law.** A point charge gives `φ ~ a/r` with `a = 0.0789` vs exact
  `1/(4π) = 0.0796` (`0.8%`), correlation `0.9998`.
- **B. Configuration energy and force.** Two charges have an energy `U(d)` that
  depends on separation; the force `F = −dU/dd` is repulsive (like charges) and
  decays with `d` (`2.04 → 0.72` over `d=0.17→0.35`) — the Coulomb `~1/d²` trend.
  This is the **seed of an inter-molecular force**.
- **C. Linearity.** Superposition is exact (`3.9e-16`): the Coulomb operator is
  linear, like all of `L`.

So the *force* the insight points to — the thing that distinguishes honey from water
— is computable in the notebook's own one operator.

---

## 2. The conceptual bridge to viscosity (named, not yet verified)

Why does the Coulomb force set `ν`? The physically correct chain is:

1. **Molecular forces ← Coulomb.** Each molecule feels electrostatic forces from its
   neighbours (§1). Honey's sugars carry many hydroxyl groups → a dense
   **hydrogen-bond network** (Coulomb-derived); water has far fewer per volume.
2. **Stress ← forces.** The microscopic stress (momentum flux) is built from these
   inter-molecular forces × positions (the virial).
3. **Viscosity ← stress correlations (Green–Kubo).** The macroscopic shear viscosity
   is the time-integral of the stress autocorrelation,
   $$
   \eta = \frac{1}{V k_B T}\int_0^\infty \big\langle \sigma_{xy}(0)\,\sigma_{xy}(t)\big\rangle\,dt .
   $$
   Stronger, longer-lived inter-molecular (Coulomb) correlations → a larger integral
   → higher viscosity. **Same mass, stronger Coulomb network ⇒ higher `ν`** — exactly
   the honey/water contrast, now with a mechanism.
4. **Back to the notebook.** `ν` then enters the verified decay rate of `06d`,
   `rate = 2ν·λ_h`, where `λ_h=(uᵀKu)/(uᵀMu)` already ties `ν` to the **inertial mass
   `M`**. So the full picture is a chain: **Coulomb (`K`) → molecular forces → stress
   correlations → `ν` → (with `M`) → the flow's decay** — viscosity and inertial mass
   as two faces, and Coulomb as the root of the viscous face.

---

## 3. The honest boundary — what this notebook can and cannot verify

This is where 空論を重ねない is decisive. The notebook is a **continuum P1 finite-
element** engine on one operator `L = M⁻¹K`. Steps 1–2 of the chain (Coulomb force,
configuration energy) are *in* that engine and are verified (§1). Steps 3–4 (stress
autocorrelation, Green–Kubo, the actual number `ν` for honey vs water) require a
**molecular-dynamics / statistical-mechanics** engine — many discrete molecules,
thermal averaging, time-correlation integrals — which is a *different machine* from
this continuum solver. 

- **Verified (in scope):** Coulomb = the operator `K`; the 1/r law; configuration-
  dependent energy and force; linearity.
- **Mapped but NOT verified (out of scope here):** the derivation of macroscopic `ν`
  from molecular Coulomb forces via Green–Kubo. Asserting a *number* for honey vs
  water from this notebook would be exactly the kind of unverified leap the project
  forbids. It is named as a direction, not claimed as a result.

The value of the insight is real and is recorded honestly: it identifies the **other
term** behind viscosity (Coulomb, not mass), shows that term **already lives in the
notebook's operator**, and draws the precise line where a second engine (MD) would
have to take over to finish the derivation.

---

## 4. What is and isn't claimed

**Verified:**
- The Coulomb interaction is the same `K` operator (1/r to `0.8%`, force `~1/d²`,
  superposition exact); the molecular force that distinguishes honey from water is
  computable in the notebook's own operator.

**Not claimed:**
- **No derivation of `ν` from first principles here.** The molecular-forces→viscosity
  bridge (Green–Kubo) needs molecular dynamics + statistical mechanics, a different
  engine; this note maps it but does not verify it. No honey/water viscosity number is
  asserted.
- No new mathematics: electrostatics-as-Poisson and Green–Kubo are classical; the
  contribution is recognising that the Coulomb origin of viscosity reuses this
  project's verified `K`, and fixing the honest scope boundary.
- The continuum `ν` used in `06c`/`06d` remains an input parameter; this note explains
  *where it comes from* physically, not how to compute it within the notebook.

---

## 5. Files

- Theory: this note, on `06c`/`06d`/`04`/`03`.
- Verification: `src/verification3d/coulomb_operator_verify.py` (checks A–C).
- Possible future direction (separate engine): a small molecular-dynamics + Green–Kubo
  study to estimate `ν` from a model Coulomb network — explicitly outside the current
  continuum notebook.
