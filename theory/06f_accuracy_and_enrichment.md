# 06f · Current Accuracy, and Reaching Machine Precision by Local Enrichment

**Status:** **verified**, `src/verification3d/fluid_highorder_accuracy_verify.py`.
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

- **Mass domain (`M`, lumped): exact.** `uᵀM_lump u = 0.2500` equals the continuous
  integral at *every* resolution — the inertial side carries **zero** error.
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
- **mass/inertia form `uᵀMu`:** exact (this note).
- **decay-rate / stiffness form `uᵀKu`:** `~5%` at `n=8`, `~1.3%` at `n=16`, `O(h²)`.
- **Coulomb 1/r potential:** `~0.8%` (`06e`).

So the *only* non-machine-precision number is the stiffness/Coulomb resolution of
non-P1 (sinusoidal) fields, and it converges as `O(h²)`.

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

Measured (1-D, smallest eigenvalue of `−d²/dx²` for the `k=2π` mode, `≈12` DOF):

| local order | eigenvalue error |
|----|------|
| P1 | 2.30e-02 |
| P2 | 1.58e-03 |
| P3 | 1.37e-04 |
| P4 | 1.29e-05 |
| P5 | 2.64e-07 |
| P6 | 3.43e-09 |
| **P8** | **2.35e-13** |

**Spectral convergence to machine precision.** With a rich enough local space the
sinusoid — not in P1 — is captured to `10⁻¹³` at fixed, modest total DOF. The intuition
is correct: extending the element's local dimension and projecting (assembling) to 3D
*does* reach machine precision. The cost is more local DOF per element (denser local
coupling), the standard high-order tradeoff.

---

## 4. So how close can we get? — the honest summary

- **Mass/inertia domain:** already exact (lumped P1).
- **Coulomb/stiffness domain:** `O(h²)` in P1 (`~1–5%`), the sole residual; → machine
  precision by **local high-order enrichment** (P8 ≈ `10⁻¹³`).
- **Divergence/energy:** already machine precision (`06c`).

The path to a near-machine-precision incompressible solver is therefore *not* a new
projection or a mass trick — it is **raising the local order of the stiffness/Coulomb
operator** on the same mesh, leaving the (already exact) mass domain alone. The two
domains being different is exactly what makes this clean: enrich `K`, keep `M`.

---

## 5. What is and isn't claimed

**Verified:**
- The current stiffness error is `O(h²)` (`5.0→1.3%`, `n=8→16`); the lumped-mass form
  is exact at every `n`; the domains are different.
- Local order enrichment drives the stiffness eigenvalue error spectrally to machine
  precision (P8 ≈ `2e-13` at ~16 DOF in 1-D).

**Not claimed:**
- A full 3-D high-order incompressible solver is **not** implemented here; the 1-D
  spectral demonstration establishes the principle, and 3-D `Pn` assembly is the
  (classical, larger) engineering step it implies — named, not built.
- No new mathematics: high-order/spectral FE and mass-lumping exactness are classical;
  the contribution is the explicit decomposition (mass exact / Coulomb `O(h²)`) on this
  project's operators and the verified route to machine precision by enrichment.
- The `O(h²)` and spectral rates are for smooth (sinusoidal) modes; non-smooth fields
  (shocks, corners) do not enjoy spectral convergence — out of scope.

---

## 6. Files

- Theory: this note, on `06d`/`06e`.
- Verification: `src/verification3d/fluid_highorder_accuracy_verify.py` (checks A–B).
- Implied future engine: a 3-D `Pn` (high-order) stiffness assembly to carry the cube
  benchmark from `~1%` to machine precision — classical, and outside the current P1
  notebook.
