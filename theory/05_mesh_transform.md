# 05 · Mesh Transform: Torus and Sphere (V2.5, theory)

**Status:** theory + **verification (PASS)**. Backed by
`src/verification3d/mesh_transform_verify.py`, which asserts all five targets of
§5. The measured numbers are recorded inline. In keeping with 空論を重ねない, the
central claim — topological deformation preserves `L`, geometric deformation
degrades it — is now a number, not a hope.

Built on `03`/`03b` (the operator `L`) and `04` (its use). This is the
shape-generalisation layer: does the physics survive a change of domain?

---

## 0. The one distinction that decides everything

A mesh has two independent parts: **connectivity** (which vertices form which
tetrahedra) and **geometry** (where the vertices sit). Deformation can change
either:

| kind | what changes | what happens to `L = M⁻¹K` |
|------|--------------|----------------------------|
| **topological** | connectivity only (identify/​glue vertices); vertices do **not** move | tet shapes unchanged → `K`, `M` weights **unchanged** → `L` quality **preserved exactly** |
| **geometric** | vertex coordinates move; connectivity fixed | tet shapes distort → cotangent weights change → `L` **degrades** (anisotropy grows) |

The cotangent / FE stiffness `K_ij` depends on tetrahedron *shape*. So moving
vertices changes `L`; merely re-gluing faces does not. This is the axis V2.5
measures.

**Grounding probe (n=6 Kuhn cube).** Tetrahedron quality `q` (normalised so a
regular tet = 1):

- baseline Kuhn cube: `q_min = q_mean = 0.657` (every tet identical).
- **geometric** cube→ball coordinate map: `q_min = 0.028`, `q_mean = 0.527`
  — the worst tet is ~23× more degenerate; `L` is heavily distorted.
- **topological** opposite-face identification (x,y): `343 → 252` vertices, tet
  shapes **unchanged**; `L` quality preserved, only the spectrum changes to the
  periodic one.

---

## 1. Why the torus is exact and the sphere is not (first principles)

This is the deep reason, and it is pure differential geometry:

- The **torus has a flat metric** — zero Gaussian curvature. A flat torus is
  exactly a **periodic box**: take the cube and identify opposite faces. No vertex
  moves, no tetrahedron distorts. The Kuhn cube *already is* a flat-torus mesh once
  its faces are glued. So the torus admits a **topological** realisation that is
  geometrically exact — `L` is the verified operator, unchanged.

- The **sphere has positive curvature**. By Gauss's *Theorema Egregium*, curvature
  is intrinsic: **no distortion-free map exists** from a flat mesh to a curved
  sphere. Any cube→sphere realisation is necessarily **geometric** and must
  distort tetrahedra. There is no "free" topological sphere the way there is a
  free topological torus.

So the two shapes sit on opposite sides of the V2.5 axis *by their geometry*:
torus = topological (exact), sphere = geometric (distortion unavoidable). That is
itself a result worth recording.

---

## 2. Torus — topological (periodic identification)

**Construction.** Identify the cube's opposite boundary vertices. For a 2-torus
(periodic in x and y, an interval in z), map every vertex with `x = +½` to the one
at `x = −½` with the same `(y,z)`, and likewise for `y`. A 3-torus glues all three
pairs. Build a remap `r: vertex → representative`, then relabel `K`, `M`, and the
tet list through `r`.

**Effect on `L`.** Because no vertex moves, every tetrahedron keeps its Kuhn shape;
the per-tet stiffness and mass are identical to the cube's. The assembled `K`, `M`
on the *identified* index set are the periodic operator — same cotangent weights,
now wired around the seam. **`L` quality is preserved exactly**; what changes is
the spectrum (periodic eigenfunctions `e^{i k·x}` with `k` on the dual lattice)
and the dynamics: a wave that reaches a face **re-enters from the opposite face**.

**Why it's clean.** The seam vertices are genuinely the same points after
identification, so there is no overlap, no gap, no special "junction" geometry —
the "接合面" is handled by relabelling, not by new elements. The grounding probe
(343→252 vertices, shapes unchanged) is exactly this.

---

## 3. Torus / sphere — geometric (move the vertices)

**Construction.** Apply a coordinate map `φ: cube → target`. For a geometric torus,
`φ` wraps the box into a ring (major radius `R`, minor radius `r`); for a ball,
`φ` pushes the cube surface out to a sphere. Connectivity is unchanged; only `V` →
`φ(V)`.

**Effect on `L` (the Jacobian story).** A P1 element's gradient operator is
`G = ∂(barycentric)/∂x`. Under `φ` with Jacobian `J = ∂φ/∂x`, the element maps to
a new shape and its stiffness recomputes from the *deformed* edge vectors. Where
`J` is near-singular or strongly anisotropic, tetrahedra flatten, cotangent weights
blow up or change sign, and `L` loses isotropy. The relevant scalar is the
**tetrahedron quality** `q = c · vol / (rms edge)³` (regular tet → 1, sliver → 0).
The probe's `q_min = 0.028` for the cube→ball map is the warning: a naive
coordinate map makes slivers.

**Honest consequence.** Deforming the cube is a *poor* way to mesh a sphere; a
purpose-built spherical tetrahedralisation would keep `q` high. Within the "reshape
the cube" scope of V2.5, the finding is precisely *how much* `L` degrades, and that
the degradation is concentrated where `φ` is most anisotropic.

---

## 4. The comparison V2.5 makes

For each target shape, put the topological and geometric realisations side by side
and measure what happens to the verified operator:

| shape | topological | geometric |
|-------|-------------|-----------|
| **torus** | periodic box: `q` unchanged, `L` exact, spectrum periodic | wrapped ring: `q` drops, `L` anisotropy grows |
| **sphere** | *does not exist* (Theorema Egregium) | cube→ball: `q` drops (probe: 0.657→0.028), measure `L` |

The headline comparison is the torus: the **same physics** (heat, wave, pressure)
run on the flat-torus mesh (exact `L`) versus a geometrically wrapped torus
(degraded `L`), quantifying the cost of moving vertices vs merely gluing them — a
direct extension of `01e`'s lesson that geometry, not relabelling, is what hurts.

---

## 5. Verification (all five PASS; `mesh_transform_verify.py`)

Run on the Kuhn cube `n=8` (729 vertices, 3072 tets):

1. **Topological torus preserves `L` + periodic spectrum.** Full opposite-face
   identification folds `729 → 512` vertices (`= 8³`, correct). Tet shapes are
   untouched, so quality is identical to baseline (`q_min = q_mean = 0.657`). The
   periodic operator has exactly **1** null eigenvalue (the constant); its lowest
   non-zero eigenvalue `37.49` matches the analytic flat-torus `(2π)² = 39.478` to
   **5.0%** (the expected FE error), with the correct **6-fold degeneracy**
   (modes `(±1,0,0),(0,±1,0),(0,0,±1)`). **PASS.**

2. **Wrap-around is seamless.** After identification every vertex has the **same
   connectivity degree (15)** — there is no boundary, so the seam vertices are
   ordinary interior vertices and a wave crossing a face simply re-enters the
   opposite one. The Laplacian row-sum stays `< 3.3e-16` at the seam. Long-wave
   plane waves (`|k|² ≤ 2`) are stationary periodic modes to **5%**; shorter waves
   show the known FE dispersion of `03b`, not a seam artefact. **PASS.**

3. **Geometric distortion metric (the headline number).**

   | mesh | `q_min` | `q_mean` | sign-flipped cotangents |
   |------|---------|----------|--------------------------|
   | baseline Kuhn cube | 0.657 | 0.657 | **0** |
   | geometric ball (move vertices) | 0.014 | 0.524 | **2844** |
   | geometric torus (wrap) | 0.016 | 0.053 | **2726** |
   | topological torus (glue) | 0.657 | 0.657 | **0** |

   Moving vertices creates **thousands of negative-conductance edges** (positive
   off-diagonal `K`), breaking the discrete maximum principle; gluing creates
   **none**. This is the quantified cost of geometry vs relabelling. **PASS.**

4. **Sphere obstruction is real.** Across several cube→sphere maps the best
   achievable `q_min` is `0.224` — well below the baseline `0.657`; none restores
   quality. The Theorema-Egregium obstruction is concrete: distortion is forced,
   not a bad-map artefact. **PASS.**

5. **Physics still solves on the degraded mesh.** On the geometric ball, heat
   energy decreases monotonically (`1.09e-3 → 8.09e-4`) and the wave leapfrog
   stays bounded (final/initial energy `0.337`). `L` degrades **gracefully** — it
   loses accuracy with quality but does not break. **PASS.**

---

## 5b. Can a different norm remove the negative cotangents?

A natural question (and a good instinct): the sign-flipped cotangents of §5's
geometric meshes are *weights* — can a different weight formula remove them?
`src/verification3d/cotangent_signflip_verify.py` answers with numbers; the answer
is "not for free." Four findings:

1. **Mechanism — sign-flips ⟺ obtuse dihedral angles.** A cotangent weight goes
   negative exactly when the relevant dihedral angle exceeds 90°. The Kuhn cube
   sits *exactly at 90°* (right-angle dihedrals → 0 obtuse, 0 flips), so it is on
   the knife-edge: **any** geometric deformation tips some dihedrals obtuse and
   sign-flips appear. Measured (n=6 ball): 816 flips ↔ 1260 obtuse dihedrals;
   baseline cube 0 ↔ 0. (This is the same right-angle degeneracy that once dropped
   wireframe edges in the viewer.)

2. **Formula tradeoff — removing flips by norm change loses linear-exactness.** On
   one distorted ball mesh:

   | operator | sign-flips | harmonic error (linear `u=x`) |
   |----------|-----------|-------------------------------|
   | cotangent (FE) | 816 | **2.8e-16** (exact) |
   | graph / uniform | **0** | 8.7e-2 |
   | clamp `max(w,0)` | **0** | 2.2e-2 |

   The cotangent weight is the one that reproduces linear fields exactly — the very
   property behind `02`'s interface `T=k₂/(k₁+k₂)` and `04`'s `D∘grad = K`. Uniform
   and clamped weights kill the sign-flips but **lose that exactness**. Changing the
   norm trades the verification foundation for cosmetic positivity.

3. **Mesh improvement helps quality but not flips.** Quality-guarded Laplacian
   smoothing (move a vertex only if local quality doesn't drop) raises
   `q_min 0.401→0.513` and keeps exactness at machine precision — but sign-flips do
   **not** fall (816→838). A flip is an **edge-summed** property of the dihedrals
   around an edge, and globally non-obtuse ("well-centered") tetrahedral meshes are
   a genuinely hard problem in 3D. Volume-quality and dihedral-obtuseness are
   different things.

4. **The clean escape.** Only **geometry preservation** gives both at once: the
   undeformed cube (= the topological torus's tetrahedra) has 0 sign-flips **and**
   reproduces linear fields to 3.3e-16. That is exactly why the **topological**
   torus is the clean path, and why a geometrically deformed mesh is accepted only
   with **graceful degradation** (§5 check 5), never repaired by a cheaper norm.

**Takeaway.** The instinct was right that the weights are where the negativity
lives — but the cotangent weight is load-bearing (it is what makes the operator
*exact*), so it cannot be swapped out without dismantling the verification. The
negative conductances are an honest signature of moving vertices on a right-angle
mesh, removable only by not moving them (topological) or by hard global remeshing.

## 6. What is and isn't claimed

**Verified (§5):**
- Topological deformation (gluing) preserves `L` exactly (q unchanged, 0
  sign-flipped cotangents); geometric deformation (moving vertices) degrades it
  (q_min 0.657→~0.015, thousands of sign-flipped cotangents).
- The torus admits an exact (flat, topological) realisation (verified periodic
  spectrum, 6-fold degeneracy, seamless degree-15 wrap); the sphere does not — a
  curvature obstruction (Theorema Egregium), shown concrete (best q_min 0.224).

**Not claimed:**
- No new mathematics: periodic boundary conditions, the cotangent Laplacian's
  shape dependence, and the curvature obstruction are classical. The contribution
  is to *measure*, on this project's own `L`, the exact cost of each deformation.
- V2.5 reshapes the **existing cube mesh**; it does not build optimal
  torus/sphere meshes. A purpose-built spherical tetrahedralisation (higher `q`)
  is a separate, later option if the degraded ball proves too coarse.
- Curved-space corrections (e.g. the Laplace–Beltrami operator's curvature terms
  on a true curved manifold) are **out of scope**; here the mesh lives in flat ℝ³
  and only its shape changes.

---

## 7. Files and next steps

- Theory: this document, on `03`/`03b`/`04` and the operator `L`.
- Mesh tools available: `src/ksf3d/mesh3d_uniform.py` (`kuhn_cube`),
  `src/ksf3d/fem3d.py` (`fem_laplacian`).
- Verification: `src/verification3d/mesh_transform_verify.py` (all five §5 targets
  PASS). **Next (optional):** torus/sphere viewers reusing the unified/pressure
  rendering on the transformed meshes (topological torus = exact; geometric ball =
  shows the graceful degradation).
