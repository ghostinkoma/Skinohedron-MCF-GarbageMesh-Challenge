# 07 · V3 Synthesis — One Operator, Four Domains, Honestly Bounded

This integrates the V3 fluid-dynamics arc (`06`–`06h`) and its audit (`06z`). It states
what was built, how the pieces fit into a single structure, and — kept deliberately in
view — exactly where the verified facts end and interpretation or future work begins.
It does **not** restate proofs; each section points to the document and script that
carry them.

> **Reading note (from the audit `06z`).** Every *kernel* below is verified and
> reproduces. The *narrative* that ties them together is broader than any single proof:
> the tests are idealised (mostly linear, single-mode, periodic, analytic). Where a
> phrase is interpretation rather than measurement, it is marked. The pressure is **not**
> trustworthy (inf-sup); the trustworthy fluid output is the **velocity**.

---

## 1. The one operator

Everything rests on a single discrete operator on the Kuhn-cube P1 mesh,
$$
L = M^{-1}K,
$$
with `K` the cotangent/FE stiffness (assembled from one per-tet gradient `G`) and `M`
the lumped mass. Verified across this notebook, the *same* `L` is read as:

| reading | equation | document |
|---|---|---|
| heat / diffusion | `Ṫ = −κ L T` | `02`, `03` |
| scalar wave | `p̈ = −c² L p` | `01`, `03` |
| pressure / incompressibility | `K p = ρ/Δt · ∇·u*` | `04`, `06c` |
| electrostatics / Coulomb | `K φ = M ρ` | `06e` |
| momentum (viscous) | `u̇ = −ν L u + …` | `06a`–`06c` |

One structure, many physics — the project's central bet, and the part that is most
solidly verified (the operator identities are exact and pressure-independent).

---

## 2. The domain map — mass, Coulomb, temperature, velocity

V3's result is not a faster solver; it is a **map of how the physical domains sit on
this one operator**, each with a *different and measured* error character:

**Mass `M` (inertia).** Exact. The lumped-mass form `uᵀMu` equals the continuous
integral at every resolution (`0.00%`, `06f`). The inertial domain carries **no**
discretisation error.

**Coulomb / stiffness `K`.** `O(h²)` for smooth modes. `uᵀKu` under-resolves a
sinusoid's gradient by `5.0→1.3%` (`n=8→16`, `06f`); the electrostatic `1/r` law is
recovered to `0.8%` (`06e`). **All** the smooth-mode error lives here, not in the mass.

**Temperature `T`.** Plays **two** roles:
- *Parallel domain (additive, `06g`).* Joined to velocity via Boussinesq, the fused
  system conserves a **total energy** (kinetic + potential) with drift `0.0001%` and no
  secular growth (symplectic), while incompressibility stays `1.9e-17`.
- *Glue (multiplicative, `06h`).* Through `ν(T)` it enters the stiffness as a weight,
  `K_ν = ∫ν(T)∇φ·∇ψ`, mediating Coulomb↔mass. The two-viscosity interface is
  `u = Uμ₂/(μ₁+μ₂)` to machine precision (`7.5e-16`) — the mechanical twin of `02`'s
  `T = k₂/(k₁+k₂)`. Honey thinning when heated *is* this mediation.

**Velocity `u`.** Divergence-free to machine precision (`‖Bu‖~10⁻¹⁶`) — but see §4: this
is algebraically trivial on its own. The *earned* statement is that the velocity is
**pressure-decoupled** (PSPG moves it only `~10⁻⁴`, `06c`).

**Pressure `p`.** Once PSPG-stabilised, the lid-driven cavity matches the Ghia (1982)
`Re=100` benchmark on both centrelines (`u,v` RMS `~7×10⁻³` at `n=64`, converging) with a
smooth pressure (`08`). The stabilisation block is `−τK` — the **same operator** — so the
pressure is a correlation domain on `L`, validated against an **external standard**, in
2-D and 3-D alike. This closes the `06z` "pressure-not-trustworthy" gap (for `Re=100`).

---

## 3. The machine-precision rule (audit-corrected)

A single principle organises every "machine precision vs `O(h²)`" result, stated
precisely (the `06z` audit corrected an earlier over-simplification):

1. **Linear fields → exact everywhere** (they live in P1): Couette, the interfaces
   `T=k₂/(k₁+k₂)` and `u=μ₂/(μ₁+μ₂)`.
2. **Polynomial-source problems → nodally exact** by superconvergence even though the
   field is *not* in P1: parabolic Poiseuille (`1.9e-15` at the nodes).
3. **Eigenvalue / sinusoidal modes → `O(h²)`** (no nodal exactness): Taylor–Green decay,
   thermal sine mode, `λ_h` vs `2k²`.

And the **route past (3) to machine precision** is *local enrichment* — raising the
element order (the verifiable form of "extend the simplex's local dimension, project the
sum back"): in 1-D, `P1 2e-2 → P8 2e-13` (`06f`). Enrich `K`; leave `M` (already exact).

This is the spine of the whole arc: the inertial `e^{-νλt}` / `e^{iωt}` time-structure is
exact (the scalar-wave structure, `06d`); the residual is purely the spatial `O(h²)` of
representing a sinusoid, removable by enrichment.

---

## 4. What is earned, and what is not (carried from the audit)

**Earned (verified, reproducing, pressure-independent):**
- the operator identities (§1) and the domain map's measured error split (§2);
- the machine-precision rule (§3) including the enrichment route to `2e-13`;
- the velocity's divergence-freeness *and* its pressure-decoupling (`~10⁻⁴`);
- the temperature fusion invariant (`0.0001%`) and the glue interface (`7.5e-16`).

**Not earned (conceded in `06z`, repeated here so the synthesis cannot drift):**
- **`‖Bu‖~10⁻¹⁶` is an algebraic triviality** of projection — not a quality result.
- **The equal-order P1 pressure is inf-sup-unstable / checkerboard** (corr `0.0–0.2`);
  PSPG is a patch, not a fix. The unstabilised pressure is **not trustworthy**; "full NS
  solved" holds for the velocity only.
- **The battery is idealised:** periodic Taylor–Green + asymmetric modes + analytic
  interfaces. No wall-bounded benchmark, no turbulence, no experiment.
- **The unifying narrative** ("two sides of one coin", "glue", any Higgs-like resonance)
  is interpretation on idealised tests; each has a verified core but the story is wider
  than the proof.
- **No first-principles `ν`** from molecular Coulomb dynamics (Green–Kubo/MD is out of
  scope, `06e`); Arrhenius `ν(T)` is used as a classical bridge.

**Versus state-of-the-art (so the position is not overstated):** as a *solver* this is
not competitive — SOTA has stable pressure, high order in 3-D, turbulence, adaptive
geometry, HPC scale, and experimental validation, none of which this has. The value here
is *epistemic*: a verification-first reconstruction on one operator, with every failure
kept. (`06z §4`.)

---

## 5. Remaining tasks

See `07a_open_problems.md` for the enumerated, prioritised list. In one line: the next
brick must be one that can **fail against an external standard** — a lid-driven cavity
with an inf-sup-stable / stabilised pressure — not another internal elaboration.

---

## 6. Map of the V3 documents

`06` overview · `06a` Stokes · `06b` advection · `06s` shapes · `06c` Navier–Stokes
(+velocity quality) · `06d` inertial mass & error · `06e` Coulomb · `06f` accuracy &
enrichment · `06g` thermal fusion · `06h` temperature as glue · `06z` audit · `07` this
synthesis · `07a` open problems. Viewer: `viewer/viewer_model.html` (the domain map).
