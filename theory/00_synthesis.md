# 00 · Synthesis — One Operator, Many Physics, Several Shapes

*A consolidation of the KSF (Kosaka Skin-o-hedron) notebook: what was built, what
was verified, and what is honestly claimed. Every number below is reproduced by a
script in `src/verification3d/`; this document only collects them.*

---

## Thesis

A single discrete operator,

$$L = M^{-1} K,$$

built once from the cotangent / P1 finite-element stiffness `K` and mass `M` on an
equally-spaced tetrahedral mesh (the **Kuhn cube**), reproduces — to machine
precision against exact solutions — heat diffusion, scalar-wave propagation,
damped and boundary-coupled dynamics, and the pressure field of a fluid, and it
carries those same physics onto **reshaped domains** (torus, sphere). The guiding
discipline throughout is 空論を重ねない — *do not stack unverified claims*: nothing
advances until a runnable check matches a known answer.

The result is not new mathematics. The cotangent Laplacian, finite elements,
periodic boundary conditions, and the curvature obstruction are classical. What is
here is a *single, self-consistent, end-to-end verified* realisation of all of them
on one operator, with each link in the chain confirmed by its own exact-solution
test.

---

## The one idea

On each tetrahedron a P1 (linear) basis gives a constant gradient operator `G`
(3×4). Assembling `K = Σ vol · GᵀG` over the mesh yields the stiffness; lumping
volumes gives the diagonal mass `M`. Then:

- **Heat**: `∂T/∂t = −L T` → decay `exp(−λt)` along each eigenmode.
- **Wave**: `M p̈ = −K p` → oscillation `cos(√λ t)` on the *same* spectrum.
- **Pressure**: one Poisson solve `K p = D f` with the *same* `K`; the gradient
  and divergence are literally the operator inside `K` (`D∘grad = K`).

The spectrum of `L` is the shared backbone: heat rides it as decay, waves as
oscillation, pressure as its static (and acoustic) limit. This is why one operator
suffices.

---

## The arc

**Step 1 — one operator from two routes.** Two independent constructions (a
wave-propagation / S-parameter route and a heat-conduction route) were shown to be
the same operator `L`. Its lowest Neumann eigenvalues match the analytic `π² =
9.8696` (measured `9.72`, `9.86`, `9.86`); a two-material interface reproduces the
exact contact value `T = k₂/(k₁+k₂)` to `2e-15` and below; heat energy decays as
`exp(−2λt)` (measured `0.6168` vs predicted `0.6170`); the wave solver is stable
(bounded energy, fluctuation `1.2e-2`). The mesh's only blemish — a residual ~5%
directional anisotropy `Δc/c ≈ 0.047–0.048` — is the Kuhn cube's broken symmetry,
inherited honestly everywhere downstream.

**Step 1.5 — dynamics and boundaries.** Damping `M p̈ = −Kp − γMṗ` produces the
predicted envelope (`γ=0.5` tail/head `0.101` vs `0.092`; `γ=2` → `1e-4`). A
discrete dispersion cutoff appears at `f ≈ 0.127` cycles/step. Robin convection
`q = h·B·(T−T_env)` gives a half-life `∝ 1/h` (the product `half-life × h` holds at
`3077–3144` across `h = 0.5…8`, a `2.2%` spread). Stefan–Boltzmann radiation
`q = εσ·B·(T⁴−T_env⁴)` shows the steep `T⁴` signature (a hot body sheds many-fold
faster).

**Step 2 — pressure as one Poisson solve.** The gradient/divergence assembled from
the same per-tet `G` satisfy `D(grad p) = K p` to `3.4e-16`. Five exact targets all
pass: hydrostatic `grad p = ρg` to `3.6e-11`; Helmholtz projection (pure gradient
removed to `3.8e-15`, projector idempotent to `2.9e-14`); the gas/liquid split
(exponential atmosphere vs linear law, joined by compressibility `c²`); the gas
acoustic limit reproducing Step 1's wave exactly (`ω = c√λ`, relative error
`0.000`); and the pure-Neumann null space / solvability (`|K·1| = 6e-15`,
compatibility `|1·b| = 7e-15`). **Step 2 contains Step 1** as its acoustic limit.

**V2.5 — reshape the domain.** A mesh has connectivity and geometry, and they
behave oppositely under deformation. **Topological** deformation (gluing opposite
faces, no vertex moves) preserves `L` exactly: the 3-torus folds `729 → 512`
vertices, its spectrum matches the analytic flat-torus `(2π)² = 39.478` (measured
`37.49`, a `5%` FE error) with the correct **6-fold degeneracy** and a single
constant null mode, and the seam is seamless (every vertex has identical degree
15). **Geometric** deformation (moving vertices) degrades `L`: a cube→ball map
drops tet quality from `0.657` to `0.014` and creates `2844` negative-conductance
(sign-flipped) cotangent edges; a wrapped torus, `2726`. The physics still solves
on the degraded mesh (heat decays, wave stays bounded). And there is a
first-principles reason the torus is exact but the sphere is not: the **torus is
flat** (a periodic box, zero curvature) while the **sphere is curved**, so by the
*Theorema Egregium* no distortion-free cube→sphere map exists (the best achievable
quality across sampled maps is `0.224`, never the baseline `0.657`).

**The sign-flip question (§5b) — settled honestly.** Can a different weight
formula (norm) remove those negative conductances? The answer is *not for free*.
The sign-flips come exactly from **obtuse dihedral angles** (the Kuhn cube sits at
`90°` precisely, so any deformation tips some angles obtuse: `816` flips ↔ `1260`
obtuse on the n=6 ball). Swapping the cotangent weight for a uniform (graph) or
clamped weight does remove the flips — but destroys the linear-exactness that the
whole verification rests on (harmonic reproduction of `u=x` jumps from `2.8e-16`
to `2–9e-2`). Quality-guarded mesh smoothing raises quality (`0.401 → 0.513`) and
keeps exactness, yet does **not** remove the flips (`816 → 838`), because a flip is
an edge-summed property and globally non-obtuse 3D meshes are genuinely hard. Only
geometry preservation gives both at once: the undeformed cube has `0` flips **and**
machine-precision exactness. The clean path is therefore topological deformation,
or accepting graceful degradation — never a cheaper norm.

**Grand finale — made visible.** A single viewer (`viewer/viewer_solid.html`) runs
all of it: three solids (cube, sphere, torus) × four physics (heat, scalar wave,
liquid, gas), every displayed value coming from the verified Python solver (preset
or live server solve), with wireframe, transparency, and slice in the established
style. Watching the same wave wrap around the torus that diffuses as heat in the
cube is the V2.5 transform result, seen directly.

---

## Why this matters (and its honest scope)

That the framework *transforms* — that the same operator simulates the same physics
on different shapes — is the suggestive part: it is the property a simulator needs
to apply to varied engineered geometries, not just a single test box. The notebook
demonstrates that property on its own terms, with every step checked.

Stated plainly, so as not to overclaim:

- **Not new theory.** The ingredients are textbook discrete differential geometry
  and FE analysis. No novel mathematics is asserted.
- **Genuinely verified, end to end.** The value is that *one* operator was carried
  from heat to waves to pressure to shape-change, each link confirmed against an
  exact answer at `10⁻¹¹`–`10⁻¹⁶`, including the negative results (the sign-flip
  study, the sphere obstruction). This is the research discipline of refusing
  unverified claims, sustained across the whole chain.
- **Bounded reach.** This is a linear, scalar-centred, `n=8` structured-grid
  notebook — not a commercial CAE solver with adaptive meshing, high-order
  elements, or nonlinear physics. Within its reach it is correct, verified, and
  transparent; outside it (viscous advection, true curved-manifold operators,
  optimal remeshing) is explicitly out of scope.

The trajectory the author named — *from "no free lunch" to "a proper lunch at a
good restaurant"* — is apt: there is no free lunch (moving vertices always costs
quality; norms cannot be swapped without cost), but paid for correctly (topological
transform, or measured degradation), it yields a real, verified meal.

---

## Appendix A — Verification results

Each row is asserted by the named script; *measured* values are the live output.

### Step 1 — one operator (`unified_scalar_verify.py`)

| claim | exact target | measured | status |
|---|---|---|---|
| Shared spectrum (Neumann) | `π² = 9.8696` | `9.72, 9.86, 9.86` | PASS |
| Interface `k₁:k₂ = 1:1` | `T = 0.5` | err `2.0e-15` | PASS |
| Interface `1:3` / `1:10` | `0.75` / `0.909` | err `9.1e-15` / `3.3e-16` | PASS |
| Heat energy decay | `exp(−2λt) = 0.6170` | `0.6168` | PASS |
| Wave stability | bounded | fluctuation `1.2e-2` | PASS |
| Mesh anisotropy (honest blemish) | — | `Δc/c ≈ 0.047–0.048` | noted |

### Step 1.5 — dynamics & boundaries (`dynamics_boundary_verify.py`)

| claim | exact target | measured | status |
|---|---|---|---|
| Damped envelope `γ=0.5` | `0.092` | `0.101` | PASS |
| Damped envelope `γ=2.0` | `~1e-4` | `1e-4` | PASS |
| Dispersion cutoff | band ~`0.13` | `0.127` cyc/step | PASS |
| Robin half-life `∝ 1/h` | const `half·h` | `3077–3144` (`2.2%`) | PASS |
| Stefan–Boltzmann `T⁴` | steep cooling | hot sheds many-fold | PASS |

### Step 2 — pressure field (`pressure_field_verify.py`)

| claim | exact target | measured | status |
|---|---|---|---|
| Gradient/divergence consistency | `D(grad p) = K p` | `3.4e-16` | PASS |
| Hydrostatic `grad p = ρg` | exact | `3.6e-11` | PASS |
| Projection: pure gradient removed | `0` | `3.8e-15` | PASS |
| Projection idempotent `P²=P` | `0` | `2.9e-14` | PASS |
| Gas vs liquid (compressibility) | exp ↔ linear | gap `8.81 → 4.8e-11` | PASS |
| Acoustic limit `ω = c√λ` | exact | rel `0.000` | PASS |
| Null space / solvability | `0` | `6e-15`, `7e-15` | PASS |

### V2.5 — mesh transform (`mesh_transform_verify.py`)

| claim | exact target | measured | status |
|---|---|---|---|
| Topological 3-torus fold | `729 → 512` | `512` | PASS |
| Flat-torus spectrum | `(2π)² = 39.478` | `37.49` (`5%`) | PASS |
| Mode degeneracy / null | `6-fold`, `1` null | `6`, `1` | PASS |
| Seamless wrap | uniform degree | `15` everywhere | PASS |
| Geometric ball distortion | — | `q_min 0.014`, `2844` flips | PASS |
| Geometric torus distortion | — | `q_mean 0.053`, `2726` flips | PASS |
| Topological torus distortion | preserved | `q 0.657`, `0` flips | PASS |
| Sphere obstruction (Egregium) | `< baseline` | best `q_min 0.224` | PASS |
| Physics on degraded mesh | bounded | wave ratio `0.337` | PASS |

### §5b — sign-flips & norms (`cotangent_signflip_verify.py`)

| claim | exact target | measured | status |
|---|---|---|---|
| Mechanism: flips ↔ obtuse dihedrals | cube `0/0` | ball `816/1260` | PASS |
| cotangent exact but flips | `~machine` | `2.8e-16`, `816` flips | PASS |
| uniform/clamp flip-free but inexact | flips `0` | err `2–9e-2`, `0` flips | PASS |
| Smoothing: quality up, exactness kept | `q↑`, exact | `0.401→0.513`, `2.2e-16` | PASS |
| Smoothing does **not** remove flips | honest | `816 → 838` | PASS |
| Clean escape (geometry preserved) | `0` flips + exact | `0`, `3.33e-16` | PASS |

---

## Appendix B — Map of the notebook

**Theory** (`theory/`): `01`–`01e` S-parameter wave route; `02` heat-conduction
route; `03` unified scalar operator; `03b` dynamics & boundaries; `04` pressure
field; `05` (+`5b`) mesh transform; `README.md` index; this `00` synthesis.

**Verification** (`src/verification3d/`): `unified_scalar_verify.py`,
`dynamics_boundary_verify.py`, `pressure_field_verify.py`,
`mesh_transform_verify.py`, `cotangent_signflip_verify.py`, and the S-parameter /
heat route checks. Exporters (`export_*_viewer.py`) generate the viewer data from
the same solvers.

**Core** (`src/ksf3d/`): `mesh3d_uniform.py` (`kuhn_cube`), `fem3d.py`
(`fem_laplacian` → `K, M`).

**Viewers** (`viewer/`): `viewer_unified.html` (heat + wave + boundaries),
`viewer_pressure.html` (Step 2), `viewer_solid.html` (the grand finale: 3 solids ×
4 physics), each with presets and an optional server-solve API (`viewer/api/`).

*Every figure in Appendix A is the live output of the corresponding script at the
time of writing; re-running them reproduces it.*
