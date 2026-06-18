# 03 · One Operator, Two Physics: Heat and Scalar Wave (Step 1 consolidation)

**Status:** consolidation document (theory), backed by code. This closes the
"Step 1 / Step 1e / heat" arc by stating the single structure the project arrived
at, and the verification that supports it. It is the foundation on which Step 2
(pressure) will be built.

**The one-line structure.** A single geometric operator — the cotangent / FE
Laplacian `L = M⁻¹K`, which is the heat-conductance network of `02` — supports
**both** heat and scalar waves, as two time-structures of the same `L`:

```
heat (diffusion):   M ṗ  = −K p        →  modes decay   as  exp(−λ t)
scalar wave:        M p̈  = −K p        →  modes oscillate as cos(√λ t)
                                           (same eigenvalues λ of L)
```

Everything below is what makes this honest: where `L` comes from, why it is the
*correct* operator (and the naive alternatives are not), what is verified, and the
honest limits.

---

## 1. The operator `L` (recap of the route here)

`L` is the P1 finite-element / DEC cotangent Laplacian on the tetrahedral mesh:
`K` the stiffness (conductance) matrix, `M` the (lumped) mass. Two equivalent
readings, established earlier:

- **Heat reading (`02`):** `K_ij = −c_ij`, a contact conductance `c_ij` = shared
  area over the *dual* (circumcentric) length. Heat flux `Q_ij = c_ij (T_j−T_i)`.
  This is "conduction from contact area", done with the correct geometry.
- **FE/DEC reading:** `K` is the stiffness matrix; `M⁻¹K` discretises `−∇²`.

Why `L` and not the simpler nodes the project tried first:

| operator | cube spectrum | wave/heat anisotropy Δc/c | verdict |
|---|---|---|---|
| scatter node `S_wave=2P₀−I` | — | ≈ 0.49 (h-indep.) | geometry-blind, fails (01e) |
| area-weighted `S_geo=2P_A−I` | — | ≈ 0.47 | area alone insufficient (01e) |
| naive area/distance FV | → 7.59 (wrong) | anisotropic | two-point flux inconsistent (02) |
| **cotangent / FE `L`** | **→ π² (correct)** | **≈ 0.048** | **consistent, ~10× isotropic** |

So `L` is the operator that is both *consistent* (right continuum limit) and
*nearly isotropic*; the residual ~5% is the Kuhn mesh's own broken symmetry (the
"equally-spaced tetrahedra" / mesh-quality thread), not the operator's.

---

## 2. Heat on `L` (verified)

Diffusion `M ṗ = −K p`. Two exact facts
([`heat_conduction_verify.py`](../src/verification3d/heat_conduction_verify.py)):

- **Spectrum / consistency:** lowest Neumann eigenvalue → `π²` on the unit cube.
- **Materials exact:** two-material steady conduction gives interface temperature
  `T_iface = k₂/(k₁+k₂)` to machine precision (≤ 1e-14) for all conductivity
  ratios. The material contrast that motivated the whole "S-parameter with
  materials" idea is, in the heat reading, **exact** — not anisotropy-limited.

Heat modes decay as `exp(−λ t)` with `λ` the eigenvalues of `L`.

---

## 3. Scalar wave on `L` (verified)

The scalar wave `M p̈ = −K p`, integrated by the symplectic leapfrog

```
p^{n+1} = 2 p^n − p^{n−1} − Δt² M⁻¹ K p^n ,    Δt² < 4/λ_max  (stability)
```

inherits `L`'s isotropy
([`unified_scalar_verify.py`](../src/verification3d/unified_scalar_verify.py)):

- **Dispersion anisotropy** of the wave phase speed `c(θ)=√λ(k)/|k|`:
  `Δc/c ≈ 0.048` — the *same* ~5% as heat, i.e. **~10× better than the
  geometry-blind scatter node** (`01e`). The wave anisotropy that blocked Step 1e
  was the scatter node's fault, not the wave equation's: integrate the same `L`
  and the wave is as isotropic as the heat operator.
- **Same eigenvalues:** the modes oscillate at `ω = √λ` with the *identical* `λ`
  that govern heat decay — the concrete sense in which heat and wave are one
  operator, two time-structures.

Honest caveat: a crude time-domain "wavefront radius by direction" measure reads
larger (~1.5) because it is contaminated by threshold and mesh-graph geometry; the
rigorous dispersion measure (5%) is the one to trust. The residual 5% is, again,
the Kuhn mesh symmetry — reducible by a more symmetric mesh, not by changing `L`.

---

## 4. Why this is the right place to stop and consolidate

- It is **one operator**, verified, doing two physics correctly. Adding Step 2
  (pressure) reuses the *same* `L` (the incompressibility projection is a Laplace
  solve), so the foundation transfers directly.
- Every claim is **code-checked** with exact ground truth (π², `k₂/(k₁+k₂)`,
  dispersion). Nothing rests on assertion.
- The honest limits are explicit: not novel (this is classical FE/DEC); ~5%
  residual mesh anisotropy (mesh-quality thread); scalar only (vector momentum is
  Step 3).

---

## 5. Honest scope (unchanged, restated)

- **Not novel.** `L` is the classical cotangent/FE Laplacian. The value of this
  consolidation is a *correct, verified, self-built* scalar foundation — not new
  mathematics. A commercial solver would compute the same `L`; the point of this
  project is to have reached it by hand, with every step checked.
- **Scalar only.** Heat (a scalar field) and scalar waves (pressure-like) are
  covered. Vector transport (momentum, the full Navier–Stokes content) is **not**
  here; it needs the face normals reconstructed into a vector (the `1+3` split),
  which is Step 3.
- **Mesh-limited isotropy.** The ~5% directional residual is the Kuhn mesh; a more
  symmetric tetrahedral mesh (the long-running "equally-spaced tetrahedra" idea)
  is the route to reduce it.

---

## 6. Deliverables of this consolidation

1. **Document** — this file (`03_unified_scalar.md`).
2. **Python** — [`unified_scalar_verify.py`](../src/verification3d/unified_scalar_verify.py):
   one script showing the *same* `L` gives consistent heat (π², exact materials)
   and isotropic scalar waves (dispersion ≈ 0.048), with the shared-eigenvalue
   identity.
3. **Viewer** — a GLSL viewer animating heat diffusion and scalar-wave propagation
   on the same mesh from the same `L`, with material painting.

---

### Roadmap position

```
[Step 1]    S-parameter wave        01, 01b ✓ ; 01c/01d reflectance closed
[Step 1e]   geometric consistency    01e (wave needs full L, not a patched node)
[Step 1-heat] heat route             02 ✓ (contact-area = cotangent conductance)
[Step 1-final] ONE OPERATOR L        03  ← THIS DOCUMENT (heat + scalar wave)
              unified_scalar_verify   ✓   + GLSL viewer
[Step 2]    pressure / tank          (built on the SAME L)
[Step 3]    fluid dynamics           (vector momentum; the 1+3 split)
```

_With Step 1 consolidated on a single verified operator, Step 2 (pressure on a
tank) builds directly on `L`. Theory-first, as always._
