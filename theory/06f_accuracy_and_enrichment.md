# 06f · Current Accuracy, and Reaching Machine Precision by Local Enrichment

**Status:** **verified** by `fluid_highorder_accuracy_verify.py`, **which is not published in this repository** (see `theory/README.md`, *Unpublished*); the numbers below cannot be regenerated from here.
Two questions, answered together: (1) with **mass** and **Coulomb/stiffness**
considered as separate domains, how accurate is the current scheme and *where* is the
error? (2) If the two domains are genuinely different, can **extending the simplex's
local dimension and projecting the sum back to 3D** push toward machine precision?
Both have clean, measured answers. Builds on `06d` (the rate `2ν·λ_h`) and `06e`
(Coulomb = the `K` operator).

---

## 1. The two domains are different — and only one carries error

The decay rate is `2ν·λ_h`, `λ_h = (uᵀKu)/(uᵀMu)`: a ratio of the **stiffness/Coulomb**
form `uᵀKu` (the `K` operator, `06e`) over the **mass/inertia** form `uᵀMu`. Measured
against the continuous values for the Taylor–Green mode (`uᵀKu = 2k²E`, `uᵀMu = E`,
`E=¼`):

| `n` | `uᵀKu` (stiffness) | error | `uᵀM_lump u` (mass) | error |
|----|------|------|------|------|
| 8  | 18.745 | **−5.0%** | 0.2500 | **0.00%** |
| 10 | 19.098 | **−3.2%** | 0.2500 | **0.00%** |
| 12 | 19.292 | **−2.3%** | 0.2500 | **0.00%** |
| 16 | 19.487 | **−1.3%** | 0.2500 | **0.00%** |

So the domains genuinely differ:

- **Mass domain (`M`, lumped): exact *for this test*.** `uᵀM_lump u = 0.2500` equals
  the continuous integral at *every* resolution. This is because lumped mass on a
  uniform grid is the trapezoidal rule, which is exact for a whole-period trigonometric
  integrand — a property of the **test field**, not of lumped mass in general.
  `mass_quadrature_verify.py` (published) measures it: the `06f`-type sinusoid on the
  uniform Kuhn mesh is exact (`≤9e-16`), but a non-trigonometric field `eˣ·y` on the same
  mesh gives `2.1e-2 → 2.3e-3` (n=8→24, **`O(h²)`**), and the same sinusoid on a
  jittered (irregular) mesh is not exact either (`5.0e-3 → 2.1e-5`).
- **Coulomb/stiffness domain (`K`): `O(h²)`.** `uᵀKu` under-resolves the sinusoid's
  gradient energy by `5.0→1.3%`. **The entire residual lives here.**

This is the precise sense in which the current `~2–3%` is *not* a mass effect and *not*
a time-stepping effect (`06d`): it is the stiffness/Coulomb operator under-resolving a
**sinusoid's gradient** in P1.

---

## 2. Current accuracy, stated plainly

For the Taylor–Green benchmark on the structured cube, P1 lumped-mass:

- **velocity divergence:** machine precision (`‖Bu‖~10⁻¹⁶`, `06c`) — algebraic.
- **inviscid energy conservation:** machine precision (`06c`).
- **mass/inertia form `uᵀMu`:** exact for this whole-period sinusoid on the uniform
  mesh; `O(h²)` for a generic smooth field (§1).
- **decay-rate / stiffness form `uᵀKu`:** `~5%` at `n=8`, `~1.3%` at `n=16`, `O(h²)`.
- **Coulomb 1/r potential:** `~0.8%` (`06e`).

So, **for this benchmark**, the only non-machine-precision number is the
stiffness/Coulomb resolution of non-P1 (sinusoidal) fields, and it converges as `O(h²)`.

---

## 3. Reaching machine precision — the "n-dim → 3D" idea is local enrichment

The reason P1 leaves `O(h²)` is that a **sinusoid is an eigenvalue problem with no
nodal exactness**. This is the mirror of why the earlier profiles were machine-precise,
but the precise mechanism differs by case and is worth stating exactly (an audit
correction): **(i)** *linear* fields (Couette, the interfaces `T=k₂/(k₁+k₂)` and
`u=μ₂/(μ₁+μ₂)`) are **in P1**, hence exact everywhere; **(ii)** the *parabolic*
Poiseuille profile is **not** in P1, yet is **nodally exact by superconvergence** (the
P1 solution of the constant-source 1-D Poisson problem matches the exact solution at
the nodes — verified `1.9e-15`); **(iii)** the *sinusoidal* eigenvalue modes have **no
such nodal exactness** and are `O(h²)`. So "machine-precise ⟺ in P1" is too simple: the
honest rule is linear-exact, polynomial-source-nodally-exact, eigenmode-`O(h²)`. The
remedy for case (iii) is to enrich the **local** space — add degrees of freedom inside each
simplex — which is exactly the verifiable form of "extend the tetrahedron to higher
local dimension and project the sum back to 3D." This is **high-order / spectral finite
elements**: `Pn` Lagrange shape functions per element, assembled (the "projection") on
the same mesh.

Measured — **in 1-D only** (smallest eigenvalue of `−d²/dx²` for the `k=2π` mode,
`≈12` DOF; no 2-D or 3-D enrichment is implemented in this repository):

| local order | eigenvalue error |
|----|------|
| P1 | 2.30e-02 |
| P2 | 1.58e-03 |
| P3 | 1.37e-04 |
| P4 | 1.29e-05 |
| P5 | 2.64e-07 |
| P6 | 3.43e-09 |
| **P8** | **2.35e-13** |

**Spectral convergence to machine precision — in 1-D.** With a rich enough local space
the sinusoid — not in P1 — is captured to `10⁻¹³` at fixed, modest total DOF. That the
same holds after projecting (assembling) to 3-D is the classical expectation for
high-order FE on smooth fields, but it is **not demonstrated here**: only the 1-D case
is measured. The cost is more local DOF per element (denser local
coupling), the standard high-order tradeoff.

---

## 4. So how close can we get? — the honest summary

- **Mass/inertia domain:** exact for this benchmark's whole-period sinusoid on the
  uniform mesh; `O(h²)` in general (§1).
- **Coulomb/stiffness domain:** `O(h²)` in P1 (`~1–5%`), the residual of this benchmark;
  → machine precision by **local high-order enrichment** (P8 ≈ `10⁻¹³`, **1-D only**).
- **Divergence/energy:** already machine precision (`06c`).

The path to a near-machine-precision incompressible solver is therefore *not* a new
projection or a mass trick — it is **raising the local order of the stiffness/Coulomb
operator** on the same mesh. For this benchmark the mass form needs no change (it is
exact here); in general it is `O(h²)` too and would be enriched consistently with `K`.

---

## 5. What is and isn't claimed

**Verified:**
- The current stiffness error is `O(h²)` (`5.0→1.3%`, `n=8→16`); for this benchmark's
  whole-period sinusoid the lumped-mass form is exact at every `n` (a trapezoidal-rule
  property of the test; `O(h²)` for generic fields, `mass_quadrature_verify.py`).
- Local order enrichment drives the stiffness eigenvalue error spectrally to machine
  precision **in 1-D** (P8 ≈ `2e-13` at ~16 DOF).

**Not claimed:**
- A full 3-D high-order incompressible solver is **not** implemented here; the 1-D
  spectral demonstration establishes the principle, and 3-D `Pn` assembly is the
  (classical, larger) engineering step it implies — named, not built.
- That lumped mass is exact in general — it is exact only for this kind of test field
  on a uniform mesh.
- Machine precision by enrichment in 2-D or 3-D — only the 1-D case is measured.
- No new mathematics: high-order/spectral FE and trapezoidal-rule exactness are
  classical; the contribution is the explicit decomposition (mass vs Coulomb error) on
  this project's operators for this benchmark, and the 1-D enrichment demonstration.
- The `O(h²)` and spectral rates are for smooth (sinusoidal) modes; non-smooth fields
  (shocks, corners) do not enjoy spectral convergence — out of scope.

---

## 6. Files

- Theory: this note, on `06d`/`06e`.
- Verification: `fluid_highorder_accuracy_verify.py` (checks A–B) — **unpublished**, not in this repository.
- Scope of the mass-exactness: `src/verification3d/mass_quadrature_verify.py` (checks A–C, published).
- Implied future engine: a 3-D `Pn` (high-order) stiffness assembly to carry the cube
  benchmark from `~1%` to machine precision — classical, and outside the current P1
  notebook.
