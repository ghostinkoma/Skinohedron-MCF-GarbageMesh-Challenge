# Architecture of KSF

This document explains the design decisions behind the Kosaka Skin-o-hedron
Model in plain language. No advanced mathematics required to read this.

---

## The problem KSF is trying to solve

Imagine you want to simulate heat spreading across the surface of a brain,
or airflow over a curved aircraft wing. You need to:

1. **Represent the curved surface** as a mesh of triangles
2. **Run a differential equation** (like the heat equation) on that mesh
3. **Get an accurate answer** even when the triangles are not perfectly shaped

Step 3 is where things have been hard for 30 years.

The standard tool (the "cotangent Laplacian") works beautifully on *nice* meshes
but breaks down on *ugly* ones. And real-world meshes are almost always ugly —
they come from CT scans, 3D scanners, automatic mesh generators, all of which
produce irregular triangles.

---

## The KSF approach: three ideas stacked together

### Idea 1: Layered refinement (§2 of the paper)

Instead of working with one fixed mesh, KSF defines a *sequence* of meshes
`S_0, S_1, S_2, ...` where each level is twice as fine as the previous one.
This is called a *refinement hierarchy* and it is the standard trick used in
multigrid solvers.

The requirement is that the mesh quality stays *bounded* as you refine.
This is the "uniform shape-regularity" condition. It rules out pathological
families where the mesh gets worse as you refine it.

### Idea 2: Trace-free dual tensors (§4–5 of the paper)

For each edge of the mesh, KSF attaches a 2×2 matrix (a "tensor") that captures
how the surface curves near that edge.

The raw tensor has an *isotropic* part (same in all directions) and an
*anisotropic* part (different in different directions). The isotropic part is
the main source of first-order error in the standard approach.

KSF removes the isotropic part by subtracting exactly half the trace:

```
E(edge) = E_raw(edge) - (1/2) * trace(E_raw) * identity
```

The factor 1/2 is exact because surfaces are 2-dimensional
(on an n-dimensional manifold it would be 1/n).

After this projection, the tensor is *trace-free*: `trace(E) = 0`.

This is verified to machine precision in `src/verification/s4_trace_free.py`.

### Idea 3: Honest about what converges

The first draft of the paper claimed the operator is second-order accurate
pointwise (at every single vertex). Running the experiments showed this is false.

What IS true, and what actually matters for most applications, is that:

- **Eigenvalues** (vibration frequencies) converge extremely well: ~4th order
- **Energy** (integrals of squared gradients) converges to 2nd order
- **Pointwise** values converge at 1st order on regular meshes

This distinction matters because:
- PDEs in physics are usually solved variationally (energy-based), not pointwise
- Spectral methods (like computing vibration modes) need eigenvalue accuracy
- Only a few niche applications actually need pointwise 2nd order

---

## Module structure

```
ksf/mesh.py
  └─ IcosphereMesh class: builds the S_k refinement sequence
     Methods: subdivide(), quality_stats(), jitter()

ksf/dec.py
  └─ coboundary(V, F) → (d0, d1)   [discrete exterior derivative]
     laplacian(V, F)  → (L, M)     [stiffness matrix + mass matrix]

ksf/trace_free.py
  └─ trace_free_tensor(V, F, edge) → 2×2 matrix
     check_trace_free(...)          → residual
     check_frame_invariance(...)    → residual

ksf/sph.py
  └─ harmonic(V, name) → (values, eigenvalue)
     [exact spherical harmonics Y10, Y2xy, Y2z2, Y3xyz]
```

Dependencies: only `numpy` and `scipy`. No C extensions, no compiled code.

---

## Sign convention

`L` is symmetric positive semidefinite and approximates `−Δ` (positive operator).

So for an eigenfunction with `Δf = −λf`, the discrete equation is:
```
L @ f ≈ λ * M @ f
```
i.e. `M⁻¹ @ L @ f ≈ λ * f`.

This is the standard convention in finite elements. Some papers use the opposite sign.

---

## Why not just use FEniCS / deal.II / Firedrake?

Those are production finite element frameworks — excellent for serious applications.

KSF is a *research prototype* aimed at understanding one specific mathematical
question: what does the trace-free projection actually buy you, and where are
its limits? For that question a clean, minimal implementation in ~300 lines of
NumPy/SciPy is more useful than a full FEM framework.

If KSF's ideas prove valuable, the next step is to implement them inside a
production FEM framework. The math and the verification suite here give you
exactly what you need to do that.
