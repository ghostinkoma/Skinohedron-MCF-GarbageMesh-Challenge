# 05 · Mesh Transform: Torus and Sphere (V2.5, theory)

**Status:** theory document (code is the next stage). It asks whether the verified
operator `L = M⁻¹K` carries over when the cube is reshaped into a **torus** or a
**sphere**, and — in keeping with 空論を重ねない — it separates the two ways a mesh
can be deformed, states what each does to `L`, and names the exact check for every
claim. Small grounding probes (reported inline) confirm the direction; full
verification is the next stage.

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

## 5. Verification design (checked next stage)

Each becomes a PASS/FAIL against an exact target when code lands:

1. **Topological torus preserves `L`.** Build the periodic operator; assert tet
   shapes (hence per-tet `K`,`M`) are bitwise the cube's, and that the spectrum
   matches the **analytic flat-torus eigenvalues** `λ = (2π)²(a²+b²+c²)/Lᵢ²` for
   integer modes (periodic box has a closed-form spectrum). *Target:* periodic
   eigenvalues to FE accuracy; `q` identical to baseline.

2. **Wrap-around dynamics.** On the topological torus, a wave packet leaving one
   face re-enters the opposite face; check momentum/energy continuity across the
   seam and that a plane wave `e^{ik·x}` with a dual-lattice `k` is a stationary
   mode (no spurious seam reflection).

3. **Geometric distortion metric.** For cube→ball and wrapped-torus maps, compute
   `q_min`, `q_mean`, the count of inverted/sign-flipped cotangent weights, and the
   wave anisotropy `Δc/c`. *Target:* quantify degradation vs the baseline
   `Δc/c ≈ 0.048` of `03`; record where (which region) it concentrates.

4. **Sphere obstruction is real, not a bug.** Show numerically that **no** vertex
   map from the cube keeps `q` near 1 on the sphere (sampling several maps), making
   the Theorema-Egregium obstruction concrete: distortion is forced, not an
   artefact of a bad map.

5. **Physics still solves on the degraded mesh.** Confirm the Step 1/2 solves
   (heat decay, wave oscillation, pressure Poisson) still run and converge on the
   geometric torus/ball, with accuracy tracking `q` — i.e. `L` degrades gracefully,
   it does not break.

---

## 6. What is and isn't claimed

**Design claims (to verify next stage):**
- Topological deformation (gluing) preserves `L` exactly; geometric deformation
  (moving vertices) degrades it, by an amount set by tetrahedron quality.
- The torus admits an exact (flat, topological) realisation; the sphere does not —
  a curvature obstruction (Theorema Egregium), not a meshing accident.

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
- **Next stage (code):** `src/verification3d/mesh_transform_verify.py` asserting
  the five targets of §5 (periodic spectrum, wrap-around, distortion metrics,
  sphere obstruction, graceful physics); then, optionally, torus/sphere viewers
  reusing the unified/pressure rendering on the transformed meshes.
