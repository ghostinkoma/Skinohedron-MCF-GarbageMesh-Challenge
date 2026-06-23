# 06z · V3 Verification Audit — Did We Pile Virtual-on-Virtual?

**Status:** audit record. Treats all of `06a`–`06h` as **hypothesis** and re-checks the
documents and code end-to-end, to support (or refute) the claim that V3 is *verified
fact, not virtual stacked on virtual*. Includes an honest comparison with
state-of-the-art simulation, and a Go/Stop recommendation. Re-applies the critical
("phantom-agent") lens in full.

---

## 1. Re-run of every verification script (today)

All **20** verification scripts were re-executed. Every one **reproduces**, with its
documented honest caveats intact:

| script | key result reproduced |
|---|---|
| `sparam_wave/geometry/interface/reflectance` | wave route; `interface` keeps its **documented spurious-reflection floor 0.134** (honest negative result; asserted parts pass) |
| `heat_conduction_verify` | interface `T=k₂/(k₁+k₂)` exact `2e-15`; consistency → `π²` |
| `unified_scalar_verify` | one `L=M⁻¹K` carries heat (exact materials) + wave (~5% isotropic) |
| `pressure_field_verify` | linear `p` machine-precise; `P²=P`; gas/acoustic limits |
| `dynamics_boundary_verify` | radiation cutoff band; `T⁴` shedding |
| `mesh_transform_verify` | flat-torus spectrum; curvature obstruction concrete |
| `cotangent_signflip_verify` | sign-flips ⟺ obtuse dihedrals; no cheap fix |
| `fluid_stokes / advection / shapes` | Stage A/B + cube/torus/sphere invariants |
| `fluid_navier_stokes_verify` | Stage C: incompressibility machine, energy, TG decay |
| `fluid_ns_velocity_quality_verify` | PSPG corr `0.001→0.946`, velocity moves `6.5e-5` |
| `fluid_ns_error_decomposition_verify` | rate `2ν·λ_h`; temporal removable; spatial `O(h²)` |
| `fluid_highorder_accuracy_verify` | mass exact; stiffness `O(h²)`; P8 → `2e-13` |
| `coulomb_operator_verify` | Coulomb `1/r` `0.8%`; superposition `4e-16` |
| `thermal_fusion_verify` | Boussinesq total-energy drift `0.0001%` |
| `thermal_glue_verify` | `ν(T)` interface `u=Uμ₂/(μ₁+μ₂)` machine `7e-16` |

**Finding:** the code runs and matches the documents. Nothing rotted; the negative
results are still negative and still recorded.

---

## 2. Cross-claim consistency audit

The claims triangulate rather than conflict:

- The Stage C Taylor–Green `~3%` (`06c`), the `λ_h` vs `2k²` gap (`06d`/`06f`), and the
  stiffness `uᵀKu` error (`06f`) are **the same number** seen three ways (`~3.2%` at
  `n=10`). Consistent.
- `06g` (temperature **additive** — conserved energy) and `06h` (temperature
  **multiplicative** — glue inside `K_ν`) are explicitly distinguished, not contradictory.
- The machine-precision claims were audited for mechanism and **one imprecision was
  found and fixed**: `06f` had lumped *parabolic* Poiseuille with *linear* fields as
  "in P1." A parabola is **not** in P1; Poiseuille is machine-precise by **nodal
  superconvergence** (constant-source 1-D Poisson), a different mechanism. No number was
  wrong; the *framing* was. Corrected in `06f §3`. The precise rule is now: linear →
  exact; polynomial-source → nodally exact; eigenmode → `O(h²)`.

**Finding:** internally consistent, with that single framing correction. No
contradictions; honest scoping is preserved throughout.

---

## 3. The phantom-agent critique — conceded in full, and bounded

The critique is **correct and code-faithful** (it quotes our own asserts). Stated plainly:

1. **`‖Bu‖~10⁻¹⁶` is an algebraic triviality** of any projector (`BP=0` by
   construction). **Conceded** — it is labelled exactly this way in
   `fluid_ns_velocity_quality_verify.py` and `06c`. It is *not* a claim of quality.
2. **The equal-order P1/P1 pressure is inf-sup (LBB) unstable → checkerboard.**
   **Conceded** — correlation `~0.0–0.2` with exact; PSPG is a *stabilisation patch*,
   not a resolution of inf-sup. The unstabilised pressure is **not trustworthy**.
3. **Therefore "we solve correct incompressible NS" cannot stand for the pressure.**
   **Conceded** — it holds for the **velocity** only. The velocity is shown
   pressure-decoupled (PSPG moves it `~10⁻⁴`), which is the *positive* result, but it
   does **not** rescue the pressure.
4. **Taylor–Green is too special; cavity/cylinder/turbulence untested.** **Conceded** —
   `06c` lists exactly these as untested. Asymmetric modes were added, but the standard
   wall-bounded benchmarks are absent.

What the critique does **not** overturn (and does not claim to):
- The **velocity** divergence-free constraint and the pressure-decoupling of the
  velocity are real, measured results.
- The scalar/heat/Coulomb/temperature **operator** identities and the machine-precision
  **linear/interface** results stand on their own (they are not pressure-dependent).

**Net:** the agent is right that the *pressure* and the *grand "full NS solved"* framing
are not earned. The narrower, per-document claims (velocity quality, operator identities,
domain map, enrichment route) are earned. The honest position is the narrow one.

---

## 4. Honest comparison with state-of-the-art simulation

To avoid overclaiming, where this notebook stands versus mature CFD/multiphysics:

**What state-of-the-art does that this notebook does NOT:**
- **Stable pressure:** SOTA uses inf-sup-stable pairs (Taylor–Hood P2–P1, MINI) or
  principled stabilisation (PSPG/SUPG) as standard. Here the baseline pressure is
  unreliable; PSPG appears only as a diagnostic patch.
- **High order in 3-D:** SOTA spectral/`hp`-FEM and DG reach high order routinely on
  real geometries. Here P1 only; high order is demonstrated **in 1-D**.
- **Turbulence & high Re:** SOTA has DNS/LES/RANS at high Reynolds number. Here:
  moderate-Re, laminar, smooth analytic flows.
- **Geometry & adaptivity:** SOTA uses unstructured, curved, adaptively refined meshes
  (AMR). Here: a structured Kuhn cube, `n=8–16`, plus degraded sphere/torus.
- **Scale & performance:** SOTA runs on HPC/GPU at 10⁶–10⁹ DOF with iterative/multigrid
  solvers. Here: small problems, direct solvers, one machine.
- **Boundary conditions & multiphysics maturity:** SOTA handles inflow/outflow, slip,
  immersed boundaries, FSI, combustion, MHD. Here: periodic + simple Dirichlet;
  Boussinesq only, linearised.
- **Validation:** SOTA is validated against experiments and standard benchmarks
  (lid-driven cavity across Re, cylinder Strouhal number, Ahmed body…). Here: **only
  analytic solutions** — no benchmark, no experiment.

**What this notebook has that is genuinely valuable (but is *not* a performance edge):**
- **Verification-first discipline:** every step checked against an exact solution or a
  conserved invariant; nothing advances on assertion.
- **One-operator conceptual unity:** heat, wave, pressure, Coulomb, temperature as
  readings of a single `L=M⁻¹K` — pedagogically clean and self-consistent.
- **Honest negative results preserved** (scatter-node rejection, sphere obstruction,
  sign-flips, spurious reflection floor, unreliable pressure).
- **Full reproducibility** at small scale.

**Honest positioning:** as a *solver*, this is **not competitive** with state-of-the-art
— not in accuracy order, turbulence, geometry, scale, or validation. Its value is
**epistemic**: a verification-first reconstruction showing how one operator underlies
many physics, with every claim checked and every failure kept. That is a teaching and
clarity contribution, not a computational-capability one.

---

## 5. Go / Stop — recommendation and opinion

**The honest status:** every *kernel* is verified and reproduces; the *risk* is that the
V3 narrative (`06d`–`06h`: "two sides of one coin", "glue", Higgs-like) is interpretation
layered on idealised, mostly-linear or single-mode tests. Each interpretation has a
verified numerical core, but the cumulative story is broader than the cumulative proof.

**Argument for GO (extend further):**
- The kernels are solid and compose cleanly; the domain map (mass exact / Coulomb
  `O(h²)`→enrich / temperature fuse+glue / velocity incompressible) is coherent.
- The next steps are concrete and standard (driven cavity, 3-D high-order, nonlinear
  Rayleigh–Bénard).

**Argument for STOP / consolidate (my recommendation):**
- The single most load-bearing gap — a **trustworthy pressure** and a **standard
  wall-bounded benchmark** (lid-driven cavity) — is exactly what would either validate
  or puncture the "we can simulate real flow" reading. Until that exists, extending into
  *more* physics (more `06x` notes) **risks piling interpretation on an untested base** —
  the very thing this audit is meant to guard against.
- The body is at a natural, honest summit: one operator, end-to-end verified on linear/
  interface problems to machine precision and on smooth modes to `O(h²)`, with the
  enrichment route and the domain map mapped.

**My opinion:** **consolidate now, then GO narrowly.** Concretely: (a) write the V3
synthesis that states the domain map *with* the §3/§4 caveats front-and-centre (so the
narrative cannot drift from the proof); then (b) take **one** standard benchmark —
the **lid-driven cavity** — with a **stabilised or Taylor–Hood pressure**, because that
single addition converts "velocity-only, periodic, analytic" into "a real, externally
checkable flow." If the cavity passes against published data, the foundation is earned
and broad GO is justified. If it struggles, we will have learned exactly where the
edifice meets reality — before building more on top.

This keeps faith with 空論を重ねない: the next brick is the one that can **fail against
an external standard**, not another internal elaboration.

---

## 6. Files

- This audit, over all of `06a`–`06h` and the foundation.
- Correction applied: `06f §3` (superconvergence vs P1-membership).
- Implied next brick: lid-driven cavity with inf-sup-stable / stabilised pressure.
