# Theory — index and status

The discrete-geometry foundation of the KSF project, built strictly
**theory → model → code → viewer**, with every claim backed by a re-runnable
verification (空論を重ねない). Read in order; later documents depend on earlier
ones. All work below is on the equally-spaced **Kuhn cube** (`n=8`: 729 vertices,
3072 tetrahedra) unless noted.

## The arc so far (Step 1 → Step 1.5)

| # | Document | What it establishes | Backing code | Status |
|---|----------|---------------------|--------------|--------|
| 01 | `01_sparameter_wave_propagation.md` | S-parameter / scatter view of wave propagation | `sparam_wave_verify.py` | done |
| 01b | `01b_sparameter_model.md` | the scatter model in detail | `s3d_sparameter.py` | done |
| 01c | `01c_planewave_interface.md` | plane wave across a material interface | `sparam_interface_verify.py` | done |
| 01d | `01d_directional_source.md` | directional source; reflectance closed by total-energy method (`R² = 2·refl_frac − 1`, recovered to ~1e-8) | `sparam_reflectance_verify.py` | done |
| 01e | `01e_geometric_consistency.md` | **decisive negative result**: geometry-blind scatter node has ~50% h-independent anisotropy; full FE/DEC operator is ~10× more isotropic (Δc/c ≈ 0.048) | `sparam_geometry_verify.py` | done |
| 02 | `02_heat_conduction_route.md` | heat = conduction from contact area; two-point flux is inconsistent on tets, cotangent/FE conductance is consistent (cube spectrum → π²); two-material interface `T = k₂/(k₁+k₂)` exact to 1e-15 | `heat_conduction_verify.py` | done |
| 03 | `03_unified_scalar.md` | **the consolidation**: one operator `L = M⁻¹K` drives heat (decay `exp(−λt)`) and waves (oscillation `cos(√λt)`) on the same spectrum | `unified_scalar_verify.py` | done |
| 03b | `03b_dynamics_and_boundaries.md` | **Step 1.5**: damping (`exp(−γt/2)`, γ=0 persists), dispersion cutoff (`f≈0.127` cyc/step), Robin convection (`t½ ∝ 1/h`), Stefan–Boltzmann radiation (`T⁴`); the parametrised viewer instrument | `dynamics_boundary_verify.py` | done |
| 04 | `04_pressure_field.md` | **Step 2 (theory)**: pressure as one Poisson solve `K p = D f` — hydrostatic (`p=ρgh`), incompressible projection (`∇·u=0`), and gas acoustics (`c²=dp/dρ`, = Step 1 wave) as three regimes; liquid/gas by parameters | planned: `pressure_field_verify.py` | theory only |

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

**Step 2 — pressure**: theory is written (`04_pressure_field.md`) — pressure as
one Poisson solve `K p = D f` on the verified `L`, unifying hydrostatic
(`p=ρgh`), incompressible projection (`∇·u=0`), and gas acoustics (the Step 1
wave as the compressible limit), with liquid/gas selected by parameters. **Code
is the next stage**: `pressure_field_verify.py` asserting the five exact-solution
targets in §4 of that document, then a dedicated pressure viewer
(`viewer/viewer_pressure.html`).
