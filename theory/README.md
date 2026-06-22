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
| 00a | `00a_prologue.en.md` / `.ja.md` (+ `.tex` / `.pdf`) | **Prologue (EN and JA, separate files)**: how the model was searched for — the rejected first instinct, conduction done right, the bet on one operator, the verified method; leads into V3. Points to `01e`/`02`/`03` for the settled claims | (references existing verifications) | done |
| 01 | `01_sparameter_wave_propagation.md` | S-parameter / scatter view of wave propagation | `sparam_wave_verify.py` | done |
| 01b | `01b_sparameter_model.md` | the scatter model in detail | `s3d_sparameter.py` | done |
| 01c | `01c_planewave_interface.md` | plane wave across a material interface | `sparam_interface_verify.py` | done |
| 01d | `01d_directional_source.md` | directional source; reflectance closed by total-energy method (`R² = 2·refl_frac − 1`, recovered to ~1e-8) | `sparam_reflectance_verify.py` | done |
| 01e | `01e_geometric_consistency.md` | **decisive negative result**: geometry-blind scatter node has ~50% h-independent anisotropy; full FE/DEC operator is ~10× more isotropic (Δc/c ≈ 0.048) | `sparam_geometry_verify.py` | done |
| 02 | `02_heat_conduction_route.md` | heat = conduction from contact area; two-point flux is inconsistent on tets, cotangent/FE conductance is consistent (cube spectrum → π²); two-material interface `T = k₂/(k₁+k₂)` exact to 1e-15 | `heat_conduction_verify.py` | done |
| 03 | `03_unified_scalar.md` | **the consolidation**: one operator `L = M⁻¹K` drives heat (decay `exp(−λt)`) and waves (oscillation `cos(√λt)`) on the same spectrum | `unified_scalar_verify.py` | done |
| 03b | `03b_dynamics_and_boundaries.md` | **Step 1.5**: damping (`exp(−γt/2)`, γ=0 persists), dispersion cutoff (`f≈0.127` cyc/step), Robin convection (`t½ ∝ 1/h`), Stefan–Boltzmann radiation (`T⁴`); the parametrised viewer instrument | `dynamics_boundary_verify.py` | done |
| 04 | `04_pressure_field.md` | **Step 2**: pressure as one Poisson solve `K p = D f` — hydrostatic (`p=ρgh`), incompressible projection (`∇·u=0`), and gas acoustics (`c²=dp/dρ`, = Step 1 wave) as three regimes; liquid/gas by parameters | `pressure_field_verify.py` | done |
| 05 | `05_mesh_transform.md` | **V2.5**: reshape cube to torus/sphere. Topological (gluing) preserves `L` exactly (0 sign-flips); geometric (moving vertices) degrades it. Torus=flat=exact (periodic spectrum, 6-fold degeneracy); sphere=curved=forced distortion (Theorema Egregium). §5b: sign-flips ⟺ obtuse dihedrals; norm change removes them but loses linear-exactness (cotangent is load-bearing) | `mesh_transform_verify.py`, `cotangent_signflip_verify.py` | done |
| 06 | `06_fluid_dynamics.md` | **V3 (theory)**: incompressible Navier–Stokes as the verified Chorin projection of `04` + an advection–diffusion predictor. Three of four terms already verified (`L`, `D∘grad=K`, `P²=P`); only nonlinear advection `C(u)` is new. Staged Stokes→advection→full NS; verified by exact flows (Taylor–Green), MMS, conservation/convergence | planned: `fluid_*_verify.py` | theory only |
| 06a | `06a_stokes_flow.md` | **V3 Stage A (verified)**: linear Stokes flow — Couette (linear) & Poiseuille (parabolic) to machine precision, steady + transient (decay rate `νπ²`), incompressibility via the verified projection. The plumbing fixed before nonlinearity | `fluid_stokes_verify.py` | done |
| 06b | `06b_advection.md` | **V3 Stage B (verified)**: scalar advection–diffusion, the new operator `C(u)` from the same per-tet `G`. Consistent (`C·1=0`), energy-conserving (skew `φᵀC_skew φ=0`), mass-conserving (all machine precision); Gaussian variance grows at `2κt` (`0.5%`), advection carries FE dispersion (`~5%`) | `fluid_advection_verify.py` | done |
| 06c | `06c_navier_stokes.md` | **V3 Stage C (verified, honestly bounded)**: full incompressible NS = self-advecting Chorin projection. `‖Bu‖~1e-16` is algebraically trivial (not the claim). Real result: **velocity is pressure-decoupled** — PSPG fixes the (checkerboard) pressure `corr 0.0→0.94` while moving velocity only `~1e-4`; Taylor–Green `4νk²` (3.1%, convergent); inviscid drift 0.000%. Pressure unreliable; cavity/cylinder/turbulence untested | `fluid_navier_stokes_verify.py`, `fluid_ns_velocity_quality_verify.py` | done |
| 06d | `06d_inertial_mass_and_error.md` | **Inertial mass & the `e^{-λt}` structure (verified)**: decomposes the Taylor–Green `~3%`. TG advection is ~99% gradient (projected out) → velocity is a viscous mode decay; rate = `2ν·λ_h`, `λ_h=(uᵀKu)/(uᵀMu)` (**mass & viscosity one coin**). Temporal error removable (all integrators agree); residual is `O(h²)` `λ_h` vs `2k²` (`5.0→1.3%`, n=8→16) — because `sin∉P1`, the mirror of why linear/parabolic profiles WERE machine-precise | `fluid_ns_error_decomposition_verify.py` | done |
| 06e | `06e_coulomb_and_viscosity.md` | **Coulomb & viscosity (part verified)**: the *other* term behind viscosity. Honey vs water — equal mass, different `ν` — points to inter-molecular **Coulomb** force. **Verified**: Coulomb = the same `K` operator (1/r to 0.8%, force ~1/d², superposition `4e-16`). **Mapped, not verified**: molecular forces → Green–Kubo → `ν` (needs MD, a different engine — honest scope boundary). Chain: Coulomb(`K`)→forces→`ν`→(with `M`)→flow | `coulomb_operator_verify.py` | done |
| 06f | `06f_accuracy_and_enrichment.md` | **Accuracy & enrichment (verified)**: mass and Coulomb domains differ — **lumped mass exact** (0.00% at every n), **stiffness/Coulomb O(h²)** (5.0→1.3%, n=8→16) carries all error. Local high-order enrichment (the 'n-dim simplex → 3D projection' = spectral FE) reaches **machine precision**: P1 2e-2 → P8 2e-13. Route to machine precision = enrich K, keep M | `fluid_highorder_accuracy_verify.py` | done |
| 06g | `06g_thermal_fusion.md` | **Thermal domain fusion (verified)**: temperature joins the unified operator with the same error split (linear T exact `9e-15`, sine `O(h²)` like Coulomb). **Fusion gain**: velocity⊕temperature (Boussinesq internal gravity wave) conserves **total energy** (drift `0.0001%`, no secular drift) while KE↔PE exchange, incompressibility `1.9e-17` maintained. Fusion adds exact invariants; enrichment (06f) lifts smooth accuracy of all domains | `thermal_fusion_verify.py` | done |
| 06h | `06h_temperature_as_glue.md` | **Temperature as glue (verified)**: temperature MEDIATES Coulomb↔mass (not just parallel, 06g). Multiplicative: `ν(T)` enters stiffness as weight `K_ν` (const T → `ν·K` exact). Arrhenius `ν(T)=ν₀exp(Eₐ/T)` — same ΔT, honey 6.2× vs water 1.7× via Coulomb barrier. Two-viscosity interface `u=U·μ₂/(μ₁+μ₂)` to **machine precision** (7e-16) = mechanical twin of 02's `T=k₂/(k₁+k₂)` | `thermal_glue_verify.py` | done |
| 06s | `06s_fluid_shapes.md` | **Fluid on shapes (verified)**: the V3 fluid solves on cube/torus/sphere. Incompressibility (`‖Du‖` machine) and mass conservation are **shape-independent** — even on the degraded sphere (2844 sign-flips, projection still `1.2e-14`); advection stable (graceful). Each shape gets its natural flow (uniform / ring / solid-body rotation) | `fluid_shapes_verify.py` | done |

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
wireframe.

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
