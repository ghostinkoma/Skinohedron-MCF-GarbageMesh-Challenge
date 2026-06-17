# 01 · S-parameter Wave Propagation on a Simplicial Complex

**Status:** theory draft (no code yet). This is Step 1 of the roadmap
*S-parameter → pressure (incompressible constraint) → fluid dynamics*. This
document fixes the theory **before** any model, code, or viewer, per the project
workflow (theory → model → code → viewer).

**Scope discipline.** Everything here is classical (scattering / transmission-line
theory, lattice methods). The aim is **not** novelty; it is a correct, explicit,
checkable foundation on the simplicial substrate, written so that the later steps
(pressure, fluid) can be added one difficulty at a time. Where this reduces to a
known method, that is stated plainly.

---

## 0. Notation and the object

Let `K` be a finite simplicial complex filling a domain `Ω ⊂ ℝⁿ` (here `n = 3`,
cells are tetrahedra; the cube `kuhn_cube` is the first test domain). Write:

- `C` — the set of top-dimensional cells (tetrahedra) `σ`.
- each cell `σ` has `n+1` facets (faces); for a tetrahedron, 4 triangular faces.
- a facet `f` is either **interior** (shared by exactly two cells `σ, σ'`) or
  **boundary** (belongs to one cell).

The facets are the **ports**. A tetrahedron is a 4-port element; this is the
"S-parameter tetrahedron" already built and verified in
[`s3d_sparameter.py`](../src/verification3d/s3d_sparameter.py). This document
explains *why* that object is the right local building block and how the local
ports assemble into a global propagation law.

---

## 1. Why facets are ports (the theory)

A conservation law in integral form is the starting point, not the differential
equation. For any cell `σ` and a conserved density `q` with flux `J`,

```
d/dt ∫_σ q dV  =  − ∮_∂σ J · da  =  − Σ_f  ∮_f J · da .          (1)
```

The cell talks to the rest of the world **only through its facets**: the time
change of the cell's content is the sum of facet fluxes. That is exactly the
definition of a port — a place where a quantity enters or leaves. So the port
structure is not a modelling choice bolted on; it is what the integral
conservation law (1) already says on a cell complex.

This is also the discrete exterior calculus (DEC) statement: the boundary
operator `∂` maps a cell to its facets, and `Σ_f` in (1) is the discrete `∮_∂σ`.
The relation `∂∂ = 0` (verified in
[`dec3d.py`](../src/ksf3d/dec3d.py) / [`s3d_complex.py`](../src/verification3d/s3d_complex.py))
is the discrete analogue of "the boundary of a boundary is empty," and it is what
will later guarantee that fluxes shared between cells cancel correctly.

**Wave variables.** Instead of tracking `(q, J)` directly, scattering theory uses
the **characteristic (wave) variables** at each port: an incoming amplitude `a_f`
and an outgoing amplitude `b_f`, defined so that the energy crossing the port is

```
power_f  =  |a_f|²  −  |b_f|² .                                   (2)
```

For a scalar wave with field value `p` and flux `v` (pressure/velocity, voltage/
current, …) related to a port impedance `Z`, the standard definitions are

```
a_f = ½ (p + Z v_f) / √Z ,     b_f = ½ (p − Z v_f) / √Z ,        (3)
```

with `v_f` the outward flux at facet `f`. The choice (3) is what makes (2) hold
and is the bridge between the physical pair `(p, v)` and the scattering pair
`(a, b)`. `Z` is the **material parameter**; §4 makes it the knob.

---

## 2. The scattering matrix and where conservation comes from

A single cell relates its outgoing port amplitudes to its incoming ones by a
linear map — the **scattering matrix** `S`:

```
b = S a ,        a = (a_1,…,a_{n+1})ᵀ ,   b = (b_1,…,b_{n+1})ᵀ .  (4)
```

Two structural properties of `S` encode physics:

- **Reciprocity:** `S = Sᵀ` (the medium has no preferred direction of signalling).
- **Losslessness (energy conservation):** `Sᴴ S = I`. Then, summing (2) over ports,
  `Σ_f |a_f|² = Σ_f |b_f|²` — the cell neither creates nor destroys energy.

For a cell with full tetrahedral symmetry (`T_d`), reciprocity + losslessness +
symmetry force

```
S = e^{iα} P₀ + e^{iβ} (I − P₀) ,     P₀ = (1/(n+1)) 𝟙𝟙ᵀ ,        (5)
```

i.e. all self-terms equal, all cross-terms equal. This is **exactly** the matrix
already constructed and verified (unitary to 1e-16, reciprocal, `1+3`
eigen-degeneracy = trivial ⊕ standard irrep of `T_d`) in
[`s3d_sparameter.py`](../src/verification3d/s3d_sparameter.py). The eigenvalue
split

```
symmetric mode  e^{iα}  (multiplicity 1)  — the "monopole" / mean
vector modes    e^{iβ}  (multiplicity n)  — the "dipole" directions          (6)
```

is the algebraic fact that will matter later: **the multiplicity-1 mode is the
seed of pressure; the multiplicity-n modes are the seed of a vector (momentum)
field.** Step 2 and Step 3 of the roadmap are, in essence, about giving these
modes their physical meaning. Here, in Step 1, we keep a single scalar wave and
do not yet separate them dynamically.

---

## 3. Assembly: from local ports to a global update

Each interior facet `f` is shared by two cells `σ, σ'`. The outgoing amplitude of
one cell at `f` is the incoming amplitude of the other at the next half-step:

```
a_f^{σ}(t+½) = b_f^{σ'}(t) ,     a_f^{σ'}(t+½) = b_f^{σ}(t) .      (7)
```

This **connection step** (7) is pure bookkeeping on the complex — it is where
`∂∂ = 0` guarantees that what leaves one cell through a face is exactly what
enters its neighbour, with no leftover. At boundary facets a boundary condition
sets `a_f` from `b_f` (e.g. `a_f = b_f` for a hard/Neumann wall = total
reflection, `a_f = −b_f` for a soft/Dirichlet wall, `a_f = 0` for an absorbing
port).

One full time step is therefore the alternation

```
   scatter:    b^{σ}(t)   = S_σ a^{σ}(t)            (local, per cell)
   connect:    a^{σ}(t+1) = swap of b across shared facets  (global, per facet)   (8)
```

This two-phase update (local scatter, global connect) is the **transmission-line
matrix (TLM)** scheme, here phrased on a simplicial complex instead of a regular
grid. The continuum limit of (8), as cell size `h → 0`, is the scalar wave
equation `∂²p/∂t² = c² Δp` with `c` fixed by `Z` and the cell geometry. The proof
is the standard TLM dispersion analysis and is **not reproduced** here; it is a
known result we will instead **verify numerically** when we reach the code step
(dispersion relation vs. exact `c`, and energy conservation `Σ|a|²` constant up to
boundary loss).

---

## 4. Material = impedance (the knob)

The only material input is the port impedance `Z(σ)` per cell. At an interior
facet between cells with `Z₁` and `Z₂`, the wave splits with reflection and
transmission coefficients

```
R = (Z₂ − Z₁)/(Z₂ + Z₁) ,        T = 2√(Z₁ Z₂)/(Z₁ + Z₂) ,        (9)
```

chosen so that **power is conserved**: with the `√Z` normalisation of (3),

```
R² + T² = 1 .                                                     (10)
```

(Equation (10) is the corrected, normalised form. The naïve `R = (Z₂−Z₁)/(Z₂+Z₁)`,
`T = 2Z₂/(Z₁+Z₂)` of voltage transmission does **not** square-sum to 1 — that was a
bug caught during exploration; the amplitude variables (3) with `T = 2√(Z₁Z₂)/(Z₁+Z₂)`
fix it.) Equal impedance `Z₁ = Z₂` gives `R = 0, T = 1` (free propagation); a hard
contrast `Z₂ ≪ Z₁` gives `R → −1` (a metal-like reflector). This single scalar
field `Z(σ)` is what makes regions behave as air, metal, wood, water, …

---

## 5. Map to the later steps (honest pointers, not promises)

Step 1 is **linear, scalar, unconstrained**. The roadmap adds one difficulty per
step; here is the precise map so nothing is oversold:

| ingredient (this doc) | Step 2 — pressure / tank | Step 3 — fluid (NS) |
|---|---|---|
| scalar amplitude on each port (3) | keep amplitude; **interpret the symmetric mode (6) as pressure** | pressure couples to the vector modes |
| `S` constant (linear) | still ~linear (acoustics in a closed tank) | `S` becomes **state-dependent** → nonlinear advection `(u·∇)u` |
| no constraint between ports | impose **∮_∂σ v·da = 0** per cell at steady state ⇒ discrete `∇·u = 0` (incompressibility) | same constraint, enforced every step (pressure projection) |
| one scalar field | one scalar (pressure) | scalar (pressure) **+** vector (momentum) = modes (6) get separate dynamics |
| `Z` = material | `Z` ↔ compressibility / tank stiffness | `Z` ↔ density; plus viscosity as a dissipative port term |

**What is genuinely reusable (the user's intuition, made precise):** the *port +
scattering + connect* skeleton (1),(4),(7)–(8) is **identical** across all three
steps. Only three things change, and they change one at a time:
(i) whether `S` depends on the state (linear → nonlinear),
(ii) whether a per-cell constraint is enforced (unconstrained → incompressible),
(iii) the **rank** of the port quantity (scalar → vector), which is exactly the
`1 vs n` eigen-split (6).

So "change a few parameters and it becomes fluid dynamics" is **half right**: the
*framework* carries over unchanged, but reaching Navier–Stokes requires adding
(i)–(iii), which are structural, not mere parameter values. This document exists
so that those additions are visible and testable rather than hidden.

---

## 6. What Step 1 must deliver (acceptance tests for the code step)

When we proceed to the model and code, the S-parameter wave step is "correct" iff:

1. **Energy:** with reflecting boundaries, `Σ_f |a_f|²` is constant in time to
   machine precision (lossless), and decays at the expected rate with absorbing
   ports.
2. **Reciprocity / symmetry:** the per-cell `S` satisfies `SᴴS = I`, `S = Sᵀ`
   (already verified for the `T_d` cell).
3. **Continuum limit:** the numerical dispersion relation approaches the exact
   wave speed `c` as `h → 0` on the cube; the pulse front travels at `c`.
4. **Material contrast:** at a planar `Z₁|Z₂` interface, measured reflection /
   transmission match (9) within discretisation error.
5. **Cube vs. ball:** the same tests on a ball-filling mesh agree with the cube in
   the interior (the user's "compare cube and sphere" check), with boundary
   differences attributable to geometry only.

Only after these pass do we move to Step 2 (pressure / tank).

---

### Roadmap position

```
[Step 1] S-parameter wave   ← THIS DOCUMENT (theory)
[Step 2] pressure / tank    (incompressible constraint ∇·u = 0)
[Step 3] fluid dynamics     (nonlinear advection; momentum = vector modes)
```

_Next artifact, on approval: `theory/01b_sparameter_model.md` — the concrete
mathematical model (state vector, exact update equations, boundary operators,
discrete energy functional) that the code will implement. No code until the model
is agreed._
