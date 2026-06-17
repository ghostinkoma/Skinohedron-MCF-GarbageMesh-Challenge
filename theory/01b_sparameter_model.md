# 01b · Mathematical Model of the S-parameter Wave Step

**Status:** mathematical model (no code yet). Follows
[`01_sparameter_wave_propagation.md`](01_sparameter_wave_propagation.md) in the
workflow theory → **model** → code → viewer. Everything here is a precise,
implementable specification: state, exact update operators, boundary operators,
and the discrete energy functional. Each definition is tagged to the acceptance
tests of `01 §6` so the code step is a faithful transcription and the
verification is "model prediction vs. numbers."

**Design decision (primary state = port wave variables).** The first-class state
lives on **oriented ports** (`a_f, b_f` per facet side), not on cell averages.
Reason: the discrete energy `Σ|a|²` closes exactly on port variables, and the
later incompressibility constraint (Step 2) is naturally "sum of facet fluxes = 0"
per cell. Cell-average fields `(p_σ, v_σ)` are **derived** quantities for display
and for material/boundary logic. This two-layer choice is fixed for the whole
roadmap.

---

## 1. Combinatorial data (given by the mesh)

From the simplicial complex `K` (cube `kuhn_cube`, later a ball mesh) extract,
once, at build time:

- `C = {0,…,N_C−1}` — tetrahedra (cells).
- `F = {0,…,N_F−1}` — facets (triangular faces).
- **Incidence with orientation.** Each facet `f` has two sides. Define the
  oriented-port set
  ```
  P = { (σ, k) : σ ∈ C,  k ∈ {0,1,2,3} } ,   |P| = 4 N_C ,
  ```
  where `k` indexes the local face of cell `σ`. Port `(σ,k)` sits on facet
  `face(σ,k) ∈ F`.
- **Facet pairing (the connection map).** A permutation `Π : P → P`:
  - interior facet shared by `(σ,k)` and `(σ',k')`  ⇒ `Π(σ,k) = (σ',k')` and
    `Π(σ',k') = (σ,k)` (an involution on that pair);
  - boundary facet: `(σ,k)` is its own partner, `Π(σ,k) = (σ,k)`, flagged
    `boundary(σ,k)=true`.

  `Π` is read directly from the adjacency already computed in
  [`export_sparameter_viewer.py`](../src/verification3d/export_sparameter_viewer.py)
  (`neighbors[σ][k]`). `Π∘Π = id` is the discrete shadow of `∂∂=0`; **the code
  must assert it** (test 01§6.2 support).

- **Geometry per port** (constants): outward unit normal `n_{σ,k}`, facet area
  `A_{σ,k}`, cell volume `V_σ`. Used only for the continuum-limit calibration of
  the wave speed (test 01§6.3) and for display.

---

## 2. State vectors

Two real amplitudes per oriented port:

```
a ∈ ℝ^{4 N_C}   incoming amplitudes,   a_{σ,k}
b ∈ ℝ^{4 N_C}   outgoing amplitudes,   b_{σ,k}
```

(Real scalar wave; complex `S` is not needed for the real wave equation. The
phase form (5) of `01` specialises to a real orthogonal `S` below.) The full
discrete state at time level `t` is the vector `a(t) ∈ ℝ^{4 N_C}`.

Per-cell material parameter (piecewise constant, paintable):

```
Z : C → ℝ_{>0} ,    Z_σ  (air=1, metal=0.01, wood=0.3, water=0.5, …).
```

Derived cell fields (display + diagnostics), with `√Z` normalisation matching
`01 (3)`:

```
p_σ = (√Z_σ / 2)  Σ_k (a_{σ,k} + b_{σ,k}) / (n+1)        (scalar "pressure")
v_{σ,k} = (1 / (2 √Z_σ)) (a_{σ,k} − b_{σ,k})              (outward facet flux)
```

---

## 3. The scatter operator `𝒮` (local, per cell)

Within a cell, outgoing = `S_σ ·` incoming. With `T_d` symmetry and the **real
lossless reciprocal** specialisation of `01 (5)`, the only one-parameter family of
real orthogonal symmetric `4×4` scattering matrices with equal diagonal / equal
off-diagonal is

```
S(θ) = cosθ · (2 P₀ − I) + sinθ · (something)…   →   we fix the physical one:

S_wave = 2 P₀ − I ,      P₀ = (1/4) 𝟙𝟙ᵀ ,                         (S)
```

i.e. explicitly

```
        [ −½  ½  ½  ½ ]
S_wave= [  ½ −½  ½  ½ ]      S_{kk} = −½ ,  S_{k≠l} = +½ .
        [  ½  ½ −½  ½ ]
        [  ½  ½  ½ −½ ]
```

Properties (all to be asserted in code, test 01§6.2):

- **orthogonal / lossless:** `S_waveᵀ S_wave = I` (since `(2P₀−I)` is a reflection).
- **reciprocal:** `S_wave = S_waveᵀ`.
- **eigenstructure (01 (6)):** `S_wave 𝟙 = +1·𝟙` (symmetric/monopole mode,
  eigenvalue +1, multiplicity 1); on `𝟙^⊥`, `S_wave = −I` (vector/dipole modes,
  eigenvalue −1, multiplicity 3). This is the real incarnation of the `1+n`
  split and the hook for Step 2/3.

`S_wave` is the **standard TLM scattering node** (a uniform node where each port
reflects −½ and transmits +½). Stated plainly: this is not new; it is the correct
known node, now on a tetrahedral cell.

The global scatter operator is block-diagonal:

```
𝒮 : ℝ^{4N_C} → ℝ^{4N_C} ,   (𝒮 a)_{σ,·} = S_wave · a_{σ,·} .       (10b)
```

Material `Z` does **not** enter `S_wave` itself; impedance contrast enters through
the **connection** step at shared facets (next section), exactly as `01 (9)`. This
separation (uniform local node + material at the interface) is deliberate and is
what keeps energy bookkeeping clean.

---

## 4. The connect operator `𝒞` (global, per facet)

After scattering, each port's outgoing amplitude becomes the partner port's
incoming amplitude, modulated by the impedance step at that facet. For an
interior port `(σ,k)` with partner `(σ',k') = Π(σ,k)`, define the reflection /
transmission from `01 (9)–(10)`:

```
Z₁ = Z_σ ,  Z₂ = Z_{σ'} ,
R_{σ,k} = (Z₂ − Z₁)/(Z₂ + Z₁) ,   T_{σ,k} = 2√(Z₁Z₂)/(Z₁ + Z₂) ,
R² + T² = 1 .                                                     (RT)
```

The connection update:

```
a_{σ,k}(t+1) = R_{σ,k} · b_{σ,k}(t)  +  T_{σ,k} · b_{σ',k'}(t) .    (C-int)
```

(The reflected part of one's own outgoing wave + the transmitted part of the
neighbour's outgoing wave.) For a **boundary** port, with wall reflection
coefficient `ρ ∈ [−1,1]`:

```
a_{σ,k}(t+1) = ρ · b_{σ,k}(t) ,                                    (C-bnd)
   ρ = +1  hard / Neumann wall (total reflection),
   ρ = −1  soft / Dirichlet wall,
   ρ =  0  absorbing port (open boundary).
```

Define the global connect operator `𝒞 : b ↦ a` by (C-int)/(C-bnd). It is linear,
sparse (≤ 2 nonzeros per row), and its matrix is **fixed once `Z` and boundary
flags are set** (recomputed only when material is painted).

**Energy property.** With equal impedance everywhere and reflecting walls,
`𝒞` is a permutation (the pure swap `Π`), hence orthogonal. With impedance
contrast, each interior facet contributes the `2×2` block

```
B_f = [ R   T ]      B_fᵀ B_f = (R²+T²) I = I      (by (RT)),
      [ T  −R ]
```

so `𝒞` is **orthogonal** ⇒ lossless, *provided* the two sides share one `(R,T)`.
**Consistency condition (must hold in code):** the block is symmetric across the
facet, i.e. `R_{σ',k'} = −R_{σ,k}` and `T_{σ',k'} = T_{σ,k}`. The code must build
`𝒞` facet-wise (one block per facet) rather than port-wise to guarantee this; an
assertion checks `‖𝒞ᵀ𝒞 − I‖` on a random state (test 01§6.1).

---

## 5. One time step

```
a(t+1) = 𝒞 𝒮 a(t)  ≡  𝒰 a(t) ,        𝒰 = 𝒞 𝒮 .                   (U)
```

`𝒰` is the full update. Both factors are orthogonal (under §3 and the §4
consistency condition), so `𝒰` is orthogonal and the scheme is **exactly
energy-preserving** in the lossless case — not approximately. This is the central
checkable claim of Step 1.

Optional physical damping (for visualisation / lossy media): replace `𝒮` by
`γ 𝒮`, `0 < γ ≤ 1`; then `‖a(t)‖² = γ^{2t} ‖a(0)‖²` exactly, a second checkable
identity.

**Source / pulse.** A pulse at cell `σ₀` sets `a_{σ₀,k} ← a_{σ₀,k} + s/2` for all
`k` (symmetric/monopole injection), energy `Σ_k (s/2)² = s²`. Used by `fire pulse`.

---

## 6. Discrete energy functional (the meter for 01§6.1)

```
E(t) = Σ_{(σ,k) ∈ P} a_{σ,k}(t)²  =  ‖a(t)‖²_2 .                  (E)
```

Claims, each a code assertion:

- **(E1) Lossless conservation.** Equal `Z`, reflecting walls, `γ=1`:
  `E(t) = E(0)` for all `t`, to machine precision (`|E(t)−E(0)|/E(0) < 10^{-12}`).
- **(E2) Damping law.** With `γ<1`: `E(t) = γ^{2t} E(0)` to machine precision.
- **(E3) Absorbing loss.** With some `ρ=0` ports: `E(t)` monotonically
  non-increasing, decreasing only through those ports (measured outflux matches
  `ΔE`).

---

## 7. Continuum limit / calibration (the meter for 01§6.3–6.4)

On the uniform cube node, TLM theory gives the propagation speed

```
c = Δx / (√3 · Δt)        (3-D uniform TLM, link-line form)        (c)
```

with `Δx` the node spacing and `Δt` one step. We do **not** re-derive this; we
**measure** it: fire a pulse, track the wavefront radius `r(t)`, fit `r = c·t`,
and compare to `c` from (c) and to the analytic acoustic speed set by `Z`. The
material-interface test (01§6.4) fires a planar wave at a `Z₁|Z₂` slab and
measures reflected/transmitted energy fractions against `R², T²` of (RT).

Acceptance order is exactly `01 §6`: (E1)→(E2/E3)→(c)→interface→(cube vs ball).

---

## 8. Exact data the code receives (interface to the code step)

The model is fully determined by these arrays (all derivable from the existing
exporter; no physics is decided in the code):

```
N_C, N_F
face(σ,k)         : 4 N_C  → F
Pi(σ,k)           : 4 N_C  → 4 N_C        (oriented partner; involution)
boundary(σ,k)     : 4 N_C  → {0,1}
normal(σ,k)       : 4 N_C  → S²           (geometry, for c-calibration/display)
area(σ,k),V_σ     : constants             (geometry)
Z_σ               : N_C    → ℝ_{>0}       (paintable material)
rho_bnd(σ,k)      : boundary ports → [−1,1]   (wall type)
gamma             : scalar ∈ (0,1]        (global damping)
state a           : 4 N_C reals           (evolved by 𝒰 = 𝒞𝒮)
```

Operators to implement, each with an assertion:

```
S_apply(a)   : block 4×4 S_wave per cell          assert ‖SᵀS−I‖<1e−12
C_build(Z,ρ) : facet-wise (R,T) / boundary blocks  assert ‖CᵀC−I‖<1e−12 (lossless case)
U_step(a)    : a ← C_apply(S_apply(a))
energy(a)    : Σ a²                                used by E1–E3
```

---

## 9. What this model is and is not

- **Is:** an exact, orthogonal, energy-conserving discrete update on the
  tetrahedral complex; the standard TLM scattering scheme phrased on simplices,
  with impedance materials and explicit boundary operators; a clean meter (E) and
  calibration (c) for verification.
- **Is not:** novel (it is TLM/scattering, known), and not yet pressure or fluid.
  The `1+3` eigen-split of `S_wave` (§3) is the *only* structural hook carried
  forward; Step 2 will promote the multiplicity-1 mode to a pressure with a
  per-cell constraint, Step 3 will make `S` state-dependent. Those are separate
  documents.

---

### Roadmap position

```
[Step 1] S-parameter wave
   01   theory          ✓
   01b  model           ← THIS DOCUMENT
   01c  code (verification suite)   ← next, on approval
   01d  GLSL viewer                 ← after code passes 01§6
[Step 2] pressure / tank
[Step 3] fluid dynamics
```

_Next artifacts, on approval: the verification code implementing §8 operators and
the §6/§7 acceptance tests (`src/verification3d/sparam_wave_*.py`), then — only
after (E1)–(E3) and the dispersion test pass — the GLSL viewer wired to the same
`𝒰 = 𝒞𝒮` update._
