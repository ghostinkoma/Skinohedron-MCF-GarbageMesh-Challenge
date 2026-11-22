# Theory — index and status

The discrete-geometry foundation of the KSF project, built strictly
**theory → model → code → viewer**, with every claim backed by a re-runnable
verification (空論を重ねない). Read in order; later documents depend on earlier
ones. All work below is on the equally-spaced **Kuhn cube** (`n=8`: 729 vertices,
3072 tetrahedra) unless noted.

## The arc so far (Step 1 → Step 1.5)

| # | Document | What it establishes | Backing code | Status |
|---|----------|---------------------|--------------|--------|
| 00 | `00_synthesis.md` / `.tex` / `.pdf` | **Synthesis**: the whole arc (one operator L → heat, wave, pressure, mesh transform) as a readable narrative + verification results appendix. Markdown body, LaTeX/PDF for print | (collects all verify scripts) | done |
| 00a | `00a_prologue.md` | **Prologue (bilingual EN/JA, restored)**: how the model was searched for — the rejected first instinct, conduction done right, the bet on one operator, the verified method. *Restored from commit `a174241`; the planned EN/JA split and `.tex`/`.pdf` were never committed* | (narrative; points to `01e`/`02`/`03`) | done |
| 01 | `01_sparameter_wave_propagation.md` | S-parameter / scatter view of wave propagation | `sparam_wave_verify.py` | done |
| 01b | `01b_sparameter_model.md` | the scatter model in detail | `s3d_sparameter.py` | done |
| 01c | `01c_planewave_interface.md` | plane wave across a material interface | `sparam_interface_verify.py` | done |
| 01d | `01d_directional_source.md` | directional source; reflectance closed by total-energy method (`R² = 2·refl_frac − 1`, recovered to ~1e-8) | `sparam_reflectance_verify.py` | done |
| 01e | `01e_geometric_consistency.md` | **decisive negative result (reconstructed)**: geometry-blind scatter node has `Δc/c ≈ 0.52 → 0.61` (n=8→16), not shrinking; FE/DEC operator `Δc/c ≈ 0.048` on the same mesh (~12× more isotropic) | `sparam_geometry_verify.py` | done (reconstructed 2026-11) |
| 02 | `02_heat_conduction_route.md` | heat = conduction from contact area; two-point flux is inconsistent on tets, cotangent/FE conductance is consistent (cube spectrum → π²); two-material interface `T = k₂/(k₁+k₂)` exact to 1e-15 | `heat_conduction_verify.py` | done |
| 03 | `03_unified_scalar.md` | **the consolidation**: one operator `L = M⁻¹K` drives heat (decay `exp(−λt)`) and waves (oscillation `cos(√λt)`) on the same spectrum | `unified_scalar_verify.py` | done |
| 03b | `03b_dynamics_and_boundaries.md` | **Step 1.5**: damping (`exp(−γt/2)`, γ=0 persists), dispersion cutoff (`f≈0.127` cyc/step), Robin convection (`t½ ∝ 1/h`), Stefan–Boltzmann radiation (`T⁴`); the parametrised viewer instrument | `dynamics_boundary_verify.py` | done |
| 04 | `04_pressure_field.md` | **Step 2**: pressure as one Poisson solve `K p = D f` — hydrostatic (`p=ρgh`), incompressible projection (`∇·u=0`), and gas acoustics (`c²=dp/dρ`, = Step 1 wave) as three regimes; liquid/gas by parameters | `pressure_field_verify.py` | done |
| 05 | `05_mesh_transform.md` | **V2.5**: reshape cube to torus/sphere. Topological (gluing) preserves `L` exactly (0 sign-flips); geometric (moving vertices) degrades it. Torus=flat=exact (periodic spectrum, 6-fold degeneracy); ball=forced distortion when the cube mesh is reused (boundary edges/corners; Liouville — the Egregium applies only to the boundary surface). §5b: sign-flips ⟺ obtuse dihedrals; norm change removes them but loses linear-exactness (cotangent is load-bearing) | `mesh_transform_verify.py`, `cotangent_signflip_verify.py` | done |
| 06 | `06_fluid_dynamics.md` | **V3 (theory)**: incompressible Navier–Stokes as the verified Chorin projection of `04` + an advection–diffusion predictor. Three of four terms already verified (`L`, `D∘grad=K`, `P²=P`); only nonlinear advection `C(u)` is new. Staged Stokes→advection→full NS; verified by exact flows (Taylor–Green), MMS, conservation/convergence | planned: `fluid_*_verify.py` | theory only |
| 06a | `06a_stokes_flow.md` | **V3 Stage A (verified)**: linear Stokes flow — Couette (linear) & Poiseuille (parabolic) to machine precision, steady + transient (decay rate `νπ²`), incompressibility via the verified projection. The plumbing fixed before nonlinearity | `fluid_stokes_verify.py` | done |
| 06b | `06b_advection.md` | **V3 Stage B (verified)**: scalar advection–diffusion, the new operator `C(u)` from the same per-tet `G`. Consistent (`C·1=0`), energy-conserving (skew `φᵀC_skew φ=0`), mass-conserving (all machine precision); Gaussian variance grows at `2κt` (`0.5%`), advection carries FE dispersion (`~5%`) | `fluid_advection_verify.py` | done |
| 06c | `06c_navier_stokes.md` | **V3 Stage C (verified, honestly bounded)**: full incompressible NS = self-advecting Chorin projection. `‖Bu‖~1e-16` is algebraically trivial (not the claim). Real result: **velocity is pressure-decoupled** — PSPG fixes the (checkerboard) pressure `corr 0.0→0.94` while moving velocity only `~1e-4`; Taylor–Green `4νk²` (3.1%, convergent); inviscid drift 0.000%. Pressure unreliable; cavity/cylinder/turbulence untested | `fluid_navier_stokes_verify.py`, `fluid_ns_velocity_quality_verify.py` | done |
| 06d | `06d_inertial_mass_and_error.md` | **Taylor–Green error decomposition (reconstructed)**: advection 99.99% gradient (projected out); decay rate → `2ν·λ_h`, `λ_h=(uᵀKu)/(uᵀMu)`, with time error `O(dt)` and removable (one scheme, three `dt`); residual `λ_h` vs `2k²` `5.0→1.3%` (n=8→16), `O(h²)` | `fluid_ns_error_decomposition_verify.py` | done (reconstructed 2026-11) |
| 06e | `06e_coulomb_and_viscosity.md` | **Coulomb & viscosity (part verified)**: the *other* term behind viscosity. Honey vs water — equal mass, different `ν` — points to inter-molecular **Coulomb** force. **Verified**: Coulomb = the same `K` operator (1/r to 0.8%, force ~1/d², superposition `4e-16`). **Mapped, not verified**: molecular forces → Green–Kubo → `ν` (needs MD, a different engine — honest scope boundary). Chain: Coulomb(`K`)→forces→`ν`→(with `M`)→flow | `coulomb_operator_verify.py` | done |
| 06f | `06f_accuracy_and_enrichment.md` | **Accuracy & enrichment (verified)**: for the Taylor–Green test the **lumped mass is exact** (0.00% at every n — a trapezoidal-rule property of a whole-period sinusoid on a uniform mesh; `O(h²)` for a generic field, `mass_quadrature_verify.py`), **stiffness/Coulomb O(h²)** (5.0→1.3%, n=8→16) carries this test's error. Local high-order enrichment reaches **machine precision in 1-D only**: P1 2e-2 → P8 2e-13 (2-D/3-D not demonstrated) | `fluid_highorder_accuracy_verify.py` (reconstructed 2026-11); `mass_quadrature_verify.py` | done |
| 06g | `06g_thermal_fusion.md` | **Thermal domain fusion (verified)**: temperature joins the unified operator with the same error split (linear T exact `9e-15`, sine `O(h²)` like Coulomb). **Fusion gain**: velocity⊕temperature (Boussinesq internal gravity wave) conserves **total energy** (drift `0.0001%`, no secular drift) while KE↔PE exchange, incompressibility `1.9e-17` maintained. Fusion adds exact invariants; enrichment (06f) lifts smooth accuracy of all domains | `thermal_fusion_verify.py` | done |
| 06h | `06h_temperature_as_glue.md` | **Temperature as glue (verified)**: temperature MEDIATES Coulomb↔mass (not just parallel, 06g). Multiplicative: `ν(T)` enters stiffness as weight `K_ν` (const T → `ν·K` exact). Arrhenius `ν(T)=ν₀exp(Eₐ/T)` — same ΔT, honey 6.2× vs water 1.7× via Coulomb barrier. Two-viscosity interface `u=U·μ₂/(μ₁+μ₂)` to **machine precision** (7e-16) = mechanical twin of 02's `T=k₂/(k₁+k₂)` | `thermal_glue_verify.py` | done |
| 06z | `06z_v3_audit.md` | **V3 audit (hypothesis re-check)**: all 20 verify scripts re-run, **all reproduce** (negative results intact). Cross-claim consistent; **one framing fix** (06f: parabola is superconvergent, not 'in P1'). Phantom-agent critique **conceded in full** (pressure unreliable/inf-sup, `‖Bu‖` trivial, benchmarks untested) and **bounded** (velocity & operator results stand). Honest SOTA comparison + **Go/Stop: consolidate, then GO narrowly via lid-driven cavity w/ stable pressure** | (meta: re-runs all) | done |
| 07 | `07_v3_synthesis.md` | **V3 synthesis (integration)**: one operator `L=M⁻¹K`, four domains (mass exact for the Taylor–Green test, O(h²) in general / Coulomb O(h²) / temperature fuse+glue / velocity div-free+pressure-decoupled), the machine-precision rule (linear-exact / superconvergence / eigenmode-O(h²)) + enrichment route (1-D only), with audit caveats front-and-centre (pressure NOT trustworthy; narrative wider than proof) | (synthesis) | done |
| 07a | `07a_open_problems.md` | **Open problems (remaining tasks)**: prioritised. **P0**: trustworthy pressure (inf-sup-stable), lid-driven cavity benchmark. **P1**: cylinder/Strouhal, nonlinear Rayleigh–Bénard, complex BCs, long-time enstrophy. **P2**: 3-D high-order, sign-flip cure, first-principles ν (MD), thermo-viscous feedback | (roadmap) | open |
| 08 | `08_driven_cavity.md` | **P0: trustworthy pressure (verified, external standard)**: PSPG-stabilised P1 lid-driven cavity matches **Ghia 1982** Re=100 — u,v centreline RMS `0.0095→0.0069` (n=48→64), **pressure smooth** (2nd-diff `0.0008`, not checkerboard). Step 2: stabilisation is `−τK`, the **same operator**. Step 3: 3-D recirculation on the same Kuhn tetrahedron. Closes the 06z pressure gap | `cavity_pspg_verify.py`, `cavity3d_demo.py` | done |
| 08a | `08a_taylorhood_crosscheck.md` | **Taylor–Hood cross-check (verified)**: removes the PSPG `τ`-tuning doubt. Inf-sup-stable P2-P1 (**no parameter**) independently matches Ghia Re=100 (RMS `0.0096→0.0062`, n=16→24). **Head-to-head PSPG vs TH: RMS 0.0015** — they agree more tightly with each other than with Ghia's table → match is correct physics, not τ artifact. P2 velocity = 06f enrichment | `cavity_taylorhood_verify.py` | done |
| 08b | `08b_high_reynolds.md` | **Re=400 & 1000, both elements + the 'rev' (verified)**: PSPG & Taylor–Hood both match Ghia (TH n=32 RMS `0.0139`/`0.0265`). Cross-check holds but **loosens** (mutual RMS `0.006→0.016`) = numerical rev. **Physical rev** (stream function, corrected — the earlier `u_min` proxy measured the primary return flow): corner eddies at Ghia's locations, counter-rotating, growing with Re (BR `×2.8`, Ghia `×2.7`); primary centre migrates to Ghia's `(0.55,0.61)→(0.53,0.56)`; eddy strengths under-resolved (BR `55→74%` of Ghia, n=64→128) and still rising. Steady laminar (transition ~Re 8000 untested) | `cavity_highre_verify.py` | done |
| 09 | `09_capstone.md` | **Capstone (integration)**: from the garbage-mesh question to a classical benchmark. Arc of `01`–`08b`; one operator carrying heat/wave/pressure/Coulomb/temperature; the moment the work met an external standard (Ghia cavity, 3 Re, 2 elements); honest bounds; and the open frontier — **why does the rev exist?** (4 candidate sources to separate by experiment) | (synthesis) | done |
| 10 | `10_rev_origin.md` | **Why does the rev exist? (reconstructed, conclusion REVISED)**: numerical rev (Ghia RMS, PSPG–TH disagreement) shrinks with refinement → approximation (1). Physical rev: corner eddies at Ghia's locations, growing with Re (`×2.8` vs Ghia `×2.7`), their strength rises slowly under refinement and converges: PSPG (n≤192) and Taylor–Hood (n≤64) extrapolate to the same `≈1.9e-3` (~10% of Ghia, order `~0.6`; §3b); the slow order is caused by the lid-corner singularity (a regularised lid gives order `~2.9`; §3c). The original "converges" series used a `u_min` proxy and is withdrawn; the convergence is re-established with the stream function. (2) missing domain / (4) chemistry not needed for the cavity | `rev_origin_verify.py`, `corner_eddy_convergence_verify.py`, `lid_singularity_verify.py` | done (reconstructed, revised 2026-11) |
| 11 | `11_pour_and_precision_viewers.md` | **Two forward markers (viewers)**: (a) `viewer_pour.html` — GLSL pour-vortices, round cup (no corner eddy) vs square masu (Moffatt corner eddies, grow with Re), **with NS equations** — *illustrative*, grounded in 10. (b) `viewer_precision.html` — **verified data**: enrichment P1→P8 (`2e-2→2e-13`), refinement O(h²), Ghia benchmark, rev separation. Answers 'does precision really improve?' = yes (approximation part) | (viewers) | done |
| 12 | `12_residual_domains.md` | **Residual domains (2)(4) + NS gaps (verified)**: tests hypotheses by residual fingerprint. **(2) Coulomb/EHD**: charged fluid drives flow pure NS can't (0→0.007, linear in charge) via `Kφ=Mρ_q` (06e). **(4) Reaction**: Fisher-KPP grows/propagates a front pure transport can't (total 0.060→0.125, residual 0.142). **(1)+(2) compose**: EHD flow converges under refinement. §4 enumerates 7 remaining NS gaps (free surface = biggest) | `residual_domains_verify.py` | done |
| 13 | `13_surface_tension.md` | **Surface tension = 5th domain (verified)**: tests 'heat/Coulomb can't describe it' — **correct**. Acts on interface GEOMETRY not a field: `κn=−ΔX` (MCF core `H=Δ_S x`). Young-Laplace pressure JUMP `Δp=σκ`; area-minimising (ellipse→circle); net force 0 (`6.6e-16`). Twist: **same Laplacian, on geometry**. Curvature `O(h²)` (ratio 4.00; caught+fixed a naive-operator non-convergence bug). NOT yet NS-coupled (needs free surface) | `surface_tension_verify.py` | done |
| 14 | `14_field_coupling.md` | **Gravity as a field on the same operator (verified)**: extends 06e. Gravity = Coulomb = same Poisson `K`, solved with the project's operator (self-gravity central pressure → `(2/3)πGρ²R²`, error `2.95→0.42%`, n=8→32); droplet pressures under zero-g/earth-g/self-g given as closed-form reference magnitudes (no NS solve); potentials superpose linearly (`3e-16`); operator dimension-independent (`uᵀKu=1` 2D&3D). Operator structure only | `field_coupling_verify.py` | done |
| 06s | `06s_fluid_shapes.md` | **Fluid on shapes (verified)**: the V3 fluid solves on cube/torus/sphere. Incompressibility (`‖Du‖` machine) and mass conservation are **shape-independent** — even on the degraded sphere (2844 sign-flips, projection still `1.2e-14`); advection stable (graceful). Each shape gets its natural flow (uniform / ring / solid-body rotation) | `fluid_shapes_verify.py` | done |

## Reconstructed and restored items (2026-11)

Four documents and four scripts were cited but never committed. They are now in the
repository, marked in their headers:

| item | how it was brought back | result vs the cited numbers |
|---|---|---|
| `00a_prologue.md` | **restored** from commit `a174241` (the author's text); Japanese 「が」 repaired | — |
| `01e_geometric_consistency.md` + `sparam_geometry_verify.py` | **reconstructed** from the citations | scatter node `Δc/c ≈ 0.52–0.61`, not shrinking (cited "~50%, h-independent": consistent); FE `0.048` (same) |
| `06d_inertial_mass_and_error.md` + `fluid_ns_error_decomposition_verify.py` | **reconstructed** | `5.0→1.3%`, `O(h²)`, rate `2νλ_h`: reproduced; "all integrators agree" narrowed to one scheme at three `dt` |
| `fluid_highorder_accuracy_verify.py` (for `06f`) | **reconstructed** | `06f` tables reproduced to the printed digits (P8: `2.62e-13` vs `2.35e-13`, round-off level) |
| `10_rev_origin.md` + `rev_origin_verify.py` | **reconstructed, conclusion revised** | numerical-rev and cavity RMS series reproduced; "physical rev converges" **withdrawn** (it used the `u_min` proxy) |

The reconstructions reproduce the cited numbers where the citation was correct, and say
so where it was not. They are re-implementations, not the lost originals.

## How to read the claims — three kinds

The results in this notebook are of three different kinds, and they support different
conclusions. Each note keeps its own wording; this table labels the main results so that
an implementation check is not read as evidence of physical accuracy, and an
interpretation is not read as a verified result.

- **[I] Implementation check** — holds by construction: the exact solution lies in the
  P1 space (patch test), or the quantity is an algebraic identity of the discrete
  operators. A pass shows the code is assembled correctly; it says nothing about how
  accurate the model is for general fields.
- **[A] Approximation accuracy** — the solution is *not* representable exactly, so the
  number measures discretisation error, preferably as a rate under refinement, or against
  an external benchmark. These are the results that test the model.
- **[N] Interpretation / analogy** — a reading or a name, not a measured quantity.

| kind | result (note) |
|---|---|
| [I] | linear fields and two-material interfaces exact: Couette, `T=k₂/(k₁+k₂)`, `u=Uμ₂/(μ₁+μ₂)` (`02`, `06a`, `06h`) |
| [I] | Poiseuille nodally exact — classical 1-D superconvergence (`06a`) |
| [I] | `D∘grad = K`, `P² = P`, `K·1 = 0`, solvability (`04`); `‖Bu‖ ~ 1e-16` after projection (`06c`) |
| [I] | `C·1 = 0`, skew-energy `φᵀC_skew φ = 0`, discrete mass conservation (`06b`, `06s`) |
| [I] | superposition of linear solves (`06e`, `14`); `K_ν = νK` at constant `T` (`06h`); `uᵀKu = 1` for `u = x` (`14`) |
| [I] | trace-free projection residual `~1e-19`, S-kernel unitarity / reciprocity (paper §4–5, §10–12; `01`) |
| [I] | lumped mass exact for a whole-period sinusoid on a uniform mesh — trapezoidal rule (`06f`; `mass_quadrature_verify.py` A) |
| [A] | 2-D sphere: pointwise `O(h^1.06)`, spectral (Rayleigh quotient) `O(h^2.00)`; irregular mesh diverges / erratic (paper §6) |
| [A] | Neumann spectrum → `π²` (`9.72`, `9.86`); wave anisotropy `Δc/c ≈ 0.048` (`02`, `03`) |
| [A] | flat-torus spectrum `5%` (`05`); damping rate vs `γ/2` `0.08–0.20%`; Robin half-life `∝ 1/h` (`03b`) |
| [A] | stiffness `uᵀKu` `5.0 → 1.3%` (n=8→16, `06f`); Taylor–Green decay `3.1%` (`06c`); advection dispersion `~5%` (`06b`) |
| [A] | stiffness `K` on a **non-periodic problem not in P1** (`u = sin(2x+1)·eʸ·cos z`, Dirichlet): `L²` `h^2.00`, `H¹` `h^1.00`, eigenvalue `h^2.00` on the uniform mesh **and** on irregular meshes with up to ~3×10⁴ sign-flipped weights (`stiffness_irregular_verify.py`) |
| [A] | lumped mass `O(h²)` for a generic field (`mass_quadrature_verify.py` B) |
| [A] | Coulomb `1/r` `0.8%` (`06e`); self-gravity central pressure `2.95 → 0.42%` (`14`); curvature `O(h²)` (`13`); reflectance `R²` recovery (`01d`) |
| [A] | **external benchmark**: lid-driven cavity vs Ghia (1982), Re = 100 / 400 / 1000, PSPG and Taylor–Hood; corner eddies located, strength converging slowly to `~10%` of Ghia (`08`, `08a`, `08b`, `10` §3b) |
| [A] | P1 → P8 enrichment to `2e-13` — **1-D only** (`06f`, `fluid_highorder_accuracy_verify.py`) |
| [N] | "one operator unifies heat / wave / pressure / Coulomb / gravity" — the shared ingredient is the Laplacian, a common mathematical form; it is not evidence that these phenomena are physically linked (`00`, `07`, `09`, `14`) |
| [N] | "Coulomb domain" as a name for the stiffness `K` (see below) |
| [N] | temperature as "glue", "two sides of one coin", any Higgs-like resemblance (`06d`, `06g`, `06h`, `07`) |
| [N] | the "rev" as a field-intrinsic structure (`09`, `10`) — supported by measurement since `10` §3b: two elements converge to the same eddy strength, but the label itself is an interpretation |

**On the name "Coulomb domain".** The stiffness `K` is the discrete Laplacian. It
appears in electrostatics (`Kφ = Mρ`), in heat conduction, and in the viscous term of
Navier–Stokes. The viscous Laplacian comes from momentum diffusion, not from the
electrostatic Poisson equation; intermolecular (Coulomb) forces set the *value* of the
viscosity, which this notebook does not compute (`06e`). Where the notes say "Coulomb
domain" or "Coulomb/stiffness", read **"stiffness `K` (the shared Laplacian)"**.

## Where the results sit in known theory

Most results in this notebook are instances of classical numerical-analysis theorems.
The table says which theorem each result follows from (or is consistent with), so that a
measured number can be read against what the theory already predicts. References are
listed below the table.

| result (note) | known theory it follows from | what the theory predicts |
|---|---|---|
| linear fields / two-material interfaces exact; `uᵀKu = 1` for `u = x` (`02`, `06a`, `06h`, `14`) | P1 reproduces affine functions; patch test [Irons & Razzaque 1972] | exact, any mesh |
| Poiseuille nodally exact (`06a`, `06f`) | nodal exactness of the 1-D Galerkin solution [Tong 1969; Douglas & Dupont 1974] | exact at the nodes for constant/polynomial sources in 1-D |
| stiffness / energy error `O(h²)`, `uᵀKu` `5.0 → 1.3%` (`06f`, `03`); `H¹` `O(h)` and `L²` `O(h²)` on uniform and irregular meshes (`stiffness_irregular_verify.py`) | Céa's lemma [Céa 1964] + P1 interpolation: `‖u−u_h‖_{H¹} = O(h)` ⇒ energy error `O(h²)` [Ciarlet 1978; Brenner & Scott 2008] | `O(h²)` for smooth fields |
| `L²` / pointwise-type errors; Coulomb `1/r` `0.8%`; self-gravity `2.95 → 0.42%` (`06e`, `14`) | Aubin–Nitsche duality [Aubin 1967; Nitsche 1968] | `O(h²)` in `L²` for smooth data; reduced near the non-smooth density jump (`14`) |
| eigenvalues: `π²` spectra, flat torus `5%`, sphere spectral `O(h^2.00)` (`02`, `03`, `05`, paper §6) | Babuška–Osborn spectral approximation [Babuška & Osborn 1991] | eigenvalue error = square of the `H¹` eigenfunction error: `O(h²)` |
| pointwise consistency only `O(h)`, divergence on irregular meshes; sign-flipped cotangents ↔ obtuse angles (paper §6–7, `05` §5b) | "no free lunch" for discrete Laplacians [Wardetzky et al. 2007]; discrete maximum principle needs non-obtuse meshes [Ciarlet & Raviart 1973] | no local, symmetric, linearly precise operator keeps the maximum principle on all meshes |
| lumped mass exact for a whole-period sinusoid on a uniform mesh, `O(h²)` otherwise (`06f`, `mass_quadrature_verify.py`) | trapezoidal-rule exactness for trigonometric polynomials [Trefethen & Weideman 2014]; mass-lumping error analysis [Ciarlet 1978] | exact for band-limited periodic integrands; `O(h²)` in general |
| time stepping: damping rate `0.08–0.20%`, dispersion cutoff `f ≈ 0.127`, stable leapfrog (`03b`, `06b`, `06c`) | Lax–Richtmyer equivalence (consistency + stability ⇔ convergence) [Lax & Richtmyer 1956]; CFL condition [Courant, Friedrichs & Lewy 1928] | convergence iff stable; explicit wave step stable below the CFL limit |
| equal-order P1/P1 pressure checkerboard (`06c`, `06z`) | inf-sup (LBB) condition [Babuška 1973; Brezzi 1974] | equal-order P1/P1 violates inf-sup ⇒ spurious pressure modes |
| PSPG-stabilised P1/P1 cavity (`08`) | PSPG circumvents inf-sup [Hughes, Franca & Balestra 1986] | stable, convergent pressure with a mesh-dependent `−τK` term (`τ ∝ h²` in the viscous limit) |
| Taylor–Hood P2/P1 cavity (`08a`) | Taylor–Hood is inf-sup stable [Taylor & Hood 1973; Brezzi & Fortin 1991] | stable without a parameter; velocity one order higher |
| cavity vs Ghia; corner eddies (`08`, `08a`, `08b`) | benchmark [Ghia, Ghia & Shin 1982]; corner eddies [Moffatt 1964] | infinite sequence of counter-rotating corner eddies, rapidly weakening |
| slow cavity convergence (order `~0.6`), restored to `~2.9` with a regularised lid (`10` §3b–§3c) | corner singularities limit the global FE convergence rate on quasi-uniform meshes; graded meshes recover it [Grisvard 1985] | rate below the smooth-solution rate unless the mesh is graded |
| P1 → P8 enrichment to `2e-13` (`06f`, 1-D only) | `p`-version / spectral element convergence [Babuška & Suri 1994] | exponential convergence in `p` for analytic solutions |
| cube → ball needs distortion (`05`) | Liouville's theorem on conformal maps in ℝ³; Gauss–Bonnet / Theorema Egregium for the boundary surface | no conformal cube → ball map exists |

**References**

- J.-P. Aubin. Behavior of the error of the approximate solutions of boundary value problems for linear elliptic operators by Galerkin's and finite difference methods. *Ann. Scuola Norm. Sup. Pisa* 21 (1967) 599–637.
- I. Babuška. The finite element method with Lagrangian multipliers. *Numer. Math.* 20 (1973) 179–192.
- I. Babuška and J. Osborn. Eigenvalue problems. In *Handbook of Numerical Analysis*, Vol. II, North-Holland, 1991, 641–787.
- I. Babuška and M. Suri. The p and h-p versions of the finite element method, basic principles and properties. *SIAM Review* 36 (1994) 578–632.
- S. C. Brenner and L. R. Scott. *The Mathematical Theory of Finite Element Methods*, 3rd ed. Springer, 2008.
- F. Brezzi. On the existence, uniqueness and approximation of saddle-point problems arising from Lagrangian multipliers. *RAIRO Anal. Numér.* 8 (1974) 129–151.
- F. Brezzi and M. Fortin. *Mixed and Hybrid Finite Element Methods*. Springer, 1991.
- J. Céa. Approximation variationnelle des problèmes aux limites. *Ann. Inst. Fourier* 14 (1964) 345–444.
- P. G. Ciarlet. *The Finite Element Method for Elliptic Problems*. North-Holland, 1978.
- P. G. Ciarlet and P.-A. Raviart. Maximum principle and uniform convergence for the finite element method. *Comput. Methods Appl. Mech. Engrg.* 2 (1973) 17–31.
- R. Courant, K. Friedrichs and H. Lewy. Über die partiellen Differenzengleichungen der mathematischen Physik. *Math. Ann.* 100 (1928) 32–74.
- J. Douglas Jr. and T. Dupont. Galerkin approximations for the two point boundary problem using continuous, piecewise polynomial spaces. *Numer. Math.* 22 (1974) 99–109.
- U. Ghia, K. N. Ghia and C. T. Shin. High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method. *J. Comput. Phys.* 48 (1982) 387–411.
- P. Grisvard. *Elliptic Problems in Nonsmooth Domains*. Pitman, 1985.
- T. J. R. Hughes, L. P. Franca and M. Balestra. A new finite element formulation for computational fluid dynamics: V. Circumventing the Babuška–Brezzi condition. *Comput. Methods Appl. Mech. Engrg.* 59 (1986) 85–99.
- B. M. Irons and A. Razzaque. Experience with the patch test for convergence of finite elements. In *The Mathematical Foundations of the Finite Element Method with Applications to Partial Differential Equations*, Academic Press, 1972, 557–587.
- P. D. Lax and R. D. Richtmyer. Survey of the stability of linear finite difference equations. *Comm. Pure Appl. Math.* 9 (1956) 267–293.
- H. K. Moffatt. Viscous and resistive eddies near a sharp corner. *J. Fluid Mech.* 18 (1964) 1–18.
- J. Nitsche. Ein Kriterium für die Quasi-Optimalität des Ritzschen Verfahrens. *Numer. Math.* 11 (1968) 346–348.
- C. Taylor and P. Hood. A numerical solution of the Navier–Stokes equations using the finite element technique. *Comput. Fluids* 1 (1973) 73–100.
- P. Tong. Exact solution of certain problems by finite-element method. *AIAA J.* 7 (1969) 178–180.
- L. N. Trefethen and J. A. C. Weideman. The exponentially convergent trapezoidal rule. *SIAM Review* 56 (2014) 385–458.
- M. Wardetzky, S. Mathur, F. Kälberer and E. Grinspun. Discrete Laplace operators: no free lunch. *Symposium on Geometry Processing* (2007) 33–37.

## The one-line result

A single geometric operator — the cotangent / FE Laplacian `L = M⁻¹K`, read as a
heat-conductance network — is both *consistent* (correct continuum limit) and
*nearly isotropic*, and it supports heat, scalar waves, damping, forcing, and
physical boundary conditions (convection, radiation) as structures on the same
`L`. The residual ~5% spatial anisotropy is the Kuhn mesh's own broken symmetry,
not the operator's.

## Instrument

`viewer/viewer_unified.html` (+ `data_unified.js`, emitted by
`src/verification3d/export_unified_viewer.py`) is the runnable instrument for the
whole arc: heatmap field on the mesh, wave/heat modes, frequency + damping,
environment (water / air / vacuum-radiation) + temperature + exchange rate,
hold-to-excite interaction, slice, transparency, and the complete tetrahedral
wireframe. It time-steps the exported, verified operator in the browser (Float32), so
its trajectories are interactive, not verified presets (see `STRUCTURE.md`).

## Next

**Step 2 — pressure**: theory written and **verified** (`04_pressure_field.md`,
`pressure_field_verify.py` — all five exact-solution targets PASS to machine
precision). Pressure is one Poisson solve `K p = D f` on the verified `L`,
unifying hydrostatic (`p=ρgh`), incompressible projection (`∇·u=0`), and gas
acoustics (the Step 1 wave as the compressible limit), liquid/gas by parameters.
**Next:** a dedicated pressure viewer is built (`viewer/viewer_pressure.html`
with `api/` server solve). **V2.5 — mesh transform** theory is written
(`05_mesh_transform.md`): reshape the cube to torus (flat → exact, topological)
and sphere (curved → forced distortion, geometric), measuring what each does to
`L`. Code next: `mesh_transform_verify.py` (periodic spectrum, wrap-around,
distortion metrics, sphere obstruction, graceful physics).
