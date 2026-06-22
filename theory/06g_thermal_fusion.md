# 06g · Fusing the Temperature Domain — the "Sum-Gain" Ticket

**Status:** **verified**, `src/verification3d/thermal_fusion_verify.py`. Takes up the
insight that decomposing physical domains and taking their **sum (相和)** is a *ticket
to fuse domains*: each domain, joined with adjoint-consistent coupling, contributes
its own exactly-conserved structure. The **temperature** domain (`02`: heat, already
in the operator) joins **mass** and **Coulomb** (`06d`/`06e`) the same way, and the
fused velocity–temperature system carries an exact invariant. Builds on
`06d`/`06e`/`06f` and `02`/`03`.

---

## 1. Temperature is already the unified operator — same error split

Heat diffusion is `∂T/∂t = κ∇²T`, the **same `K`** as Coulomb (`06e`) and pressure
(`04`). Measured, it shows the *same* split the other domains do:

- **Linear temperature `T=y`: machine precision** (`9.1e-15`, `4.1e-15` at n=8,12) —
  a P1 field, exact, exactly as Couette/Poiseuille (`06a`) and the interface
  `T=k₂/(k₁+k₂)` (`02`) were.
- **Sinusoidal thermal mode: `O(h²)`** — `λ_h` vs `π²` is `−1.3%, −0.6%, −0.3%`
  (n=8,12,16), the **stiffness/Coulomb side**, converging with refinement.

So temperature carries no new error structure: it is `M`-exact where the field is in
P1 and `K`-limited (`O(h²)`) where it is sinusoidal — *the same as mass and Coulomb.*
It fuses cleanly.

---

## 2. The fusion gain — an exactly-conserved combined invariant

Fusing two domains with **adjoint-consistent coupling** produces an invariant of the
*combined* system that is exact, even though each part oscillates. Demonstrated on the
**Boussinesq internal gravity wave** (velocity ⊕ temperature):

$$
\frac{\partial u_y}{\partial t} = \beta g\,\theta - \partial_y p,\qquad
\frac{\partial \theta}{\partial t} = -\frac{N^2}{\beta g}\,u_y,\qquad \nabla\cdot u = 0,
$$

with total energy `E = ½∫u² + ½(βg/N²)∫θ²`. The buoyancy (θ→force) and the background
advection (u_y→θ̇) are adjoint, so the two energies **exchange** while the total is
conserved. Measured (symplectic Störmer–Verlet, consistent `(B,Bᵀ)` projection):

- **total energy `0.250000 → 0.250000`, net drift `0.0001%`**, bounded swing
  `0.0002%` — no secular drift; kinetic ↔ potential exchange cleanly.
- **incompressibility `‖Bu‖ = 1.9e-17`** — machine precision, maintained throughout.

This is the **fusion gain**, precisely stated: the combined velocity–temperature
system has an **exactly conserved total energy** (to symplectic accuracy, no drift) —
the same *kind* of exact invariant as incompressibility (`06c`, from the consistent
`(B,Bᵀ)` adjoint) and skew-advection energy (`06c`). Each domain fused with an
adjoint-consistent coupling adds a machine-precision conserved structure.

---

## 3. What "accuracy improves by fusing domains" does and doesn't mean

Stated honestly, to avoid overclaiming:

- **It DOES mean** the fused system gains an **exact conserved invariant** (total
  energy here), held to machine precision by adjoint-consistent coupling + a symplectic
  integrator — structure that neither domain alone expresses. This is the real "ticket":
  domains compose because they share the one operator `L`, and each coupling done as an
  exact adjoint contributes an exact invariant.
- **It does NOT mean** the temperature domain reduces the `O(h²)` error of the
  velocity's *sinusoidal* modes (or vice versa). Each domain still has its own `O(h²)`
  for non-P1 fields. The route that lifts the smooth-mode accuracy of **all** fused
  domains together is the **local enrichment** of `06f` (P1→P8 → machine precision),
  applied to the shared stiffness operator.

So the picture is: **fusion adds exact invariants (this note); enrichment lifts smooth
accuracy (`06f`); the mass domain is already exact (`06f`).** Together they are the
program the insight predicted — every domain on the one operator, summed, each
contributing structure, all liftable to machine precision by enriching the shared `K`.

---

## 4. What is and isn't claimed

**Verified:**
- Temperature joins the unified operator with the same error split (linear exact, sine
  `O(h²)`).
- The fused velocity–temperature (Boussinesq) system conserves total energy with no
  secular drift (`0.0001%`) while incompressibility stays machine precision
  (`1.9e-17`) — an exact combined invariant from adjoint-consistent coupling.

**Not claimed:**
- Fusion does **not** reduce the per-domain `O(h²)` of sinusoidal modes; smooth-mode
  accuracy is lifted by enrichment (`06f`), not by adding domains.
- The Boussinesq run is the **linearised** internal-gravity-wave regime (clean test of
  the invariant); full nonlinear thermal convection (Rayleigh–Bénard, plumes) is a
  larger study, named not done.
- No new mathematics: Boussinesq coupling, symplectic integrators, and conserved
  energy are classical; the contribution is showing temperature fuses on this project's
  one operator and the combined invariant is exact, alongside the `06f` enrichment route.

---

## 5. Files

- Theory: this note, on `02`/`03`/`06d`/`06e`/`06f`.
- Verification: `src/verification3d/thermal_fusion_verify.py` (checks A–B).
- Implied direction: full nonlinear Rayleigh–Bénard convection on the fused operator,
  and 3-D high-order enrichment to carry all fused domains to machine precision.
