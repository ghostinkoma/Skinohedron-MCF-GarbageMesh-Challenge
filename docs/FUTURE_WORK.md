# Future Work

This document sketches directions for extending KSF beyond its current scope.

---

## 1. Navier–Stokes on curved domains

**Why this matters**: Simulating blood flow through a heart valve, airflow over
a curved wing, or ocean currents on the sphere all require solving
Navier–Stokes on curved, irregular meshes.

The Navier–Stokes equations (incompressible, viscous):

```
∂u/∂t + (u·∇)u = −∇p + ν Δu + f
∇·u = 0
```

where `u` is velocity, `p` is pressure, `ν` is viscosity, `f` is body force.

**How KSF could contribute**:

1. **Velocity Laplacian**: Use the KSF operator for `Δu` on curved surfaces.
   Already implemented — this is just `M⁻¹ L` from `ksf/dec.py`.

2. **Pressure–velocity coupling**: Need a compatible discretisation of `∇p`
   and `∇·u`. The DEC framework in `ksf/dec.py` already provides the
   coboundary operators `d0` and `d1` which are the natural discrete gradient
   and divergence. No new code needed — just wire them together.

3. **Convection term**: The `(u·∇)u` term requires upwinding on curved meshes.
   Streamline-upwind Petrov–Galerkin (SUPG) is the standard approach.
   KSF's shape-regular meshes would keep the upwinding stable.

4. **Time integration**: Standard BDF2 or Crank–Nicolson.

**Starting point** (minimal working example):

```python
# Stokes flow on sphere (linearised, no convection)
# ksf already provides L and M; this is the missing piece

from ksf.mesh import IcosphereMesh
from ksf.dec import laplacian, coboundary
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import bmat

mesh = IcosphereMesh(level=3)
V, F = mesh.vertices, mesh.faces

L, M = laplacian(V, F)          # viscous term
d0, d1, _, _ = coboundary(V, F) # gradient / divergence

# Assemble Stokes system [L  d0; d0.T  0] [u; p] = [f; 0]
# (sketch — boundary conditions and tangential projection still needed)
```

**Estimated effort**: 2–4 weeks for a working prototype on the sphere.

---

## 2. Prove Conjecture 1 (spectral super-convergence)

The numerical experiments in §6 show spectral convergence order ~3.7.
This is observed but not proved.

The standard tool is the **Babuška–Osborn spectral approximation theory**:
if the operator is self-adjoint, the mesh family satisfies an approximation
property, and the exact eigenfunction is smooth, then the eigenvalue error
is the *square* of the eigenfunction approximation error.

Since KSF gives ~O(h^1.9) in H¹ (from the variational result), the eigenvalue
error should be ~O(h^3.8) — which matches what we see.

**What's needed**: A rigorous Sobolev estimate for the KSF metric approximation
(axioms 4–5 in Definition 2.1). This is a functional-analysis problem, not
a numerics problem.

---

## 3. Extension to 3D volume meshes

The current KSF is for surfaces (2-manifolds embedded in 3D). The trace-free
projection uses the factor 1/2 because surfaces are 2-dimensional.

For a 3D volume mesh, the factor becomes 1/3, and the DEC operators become:
- `d0: 0-forms → 1-forms` (gradient)
- `d1: 1-forms → 2-forms` (curl)
- `d2: 2-forms → 3-forms` (divergence)

Everything in `ksf/dec.py` generalises straightforwardly. The main new
ingredient is 3D Hodge stars (dual cell volumes, dual face areas).

---

## 4. GPU / parallel implementation

The KSF operator assembly is embarrassingly parallel:
- Each edge is independent (1-ring neighbourhood)
- Each face is independent (area computation)
- The trace-free projection is a rank-1 update per edge

A CUDA or Metal implementation would enable real-time simulation on GPU
for meshes with millions of vertices.

---

## 5. Connection to physics simulations

The S-parameter framework in Part II opens up:

- **Quantum lattice**: S-kernels as region-to-region couplings in
  tensor-network models
- **Photonics**: Rigorous coupled-wave analysis on curved gratings
- **Acoustics**: Room impulse response on non-rectangular geometries
- **Electromagnetics**: FDTD on irregular meshes without staircasing artefacts

---

## 6. arXiv submission

The paper is ready for arXiv preprint submission.

Suggested category: `math.NA` (Numerical Analysis),
cross-list: `cs.CG` (Computational Geometry), `math.DG` (Differential Geometry).

Steps:
1. Register at https://arxiv.org/
2. Upload `paper/kosaka_skin-o-hedron_revised.tex`
3. Link this GitHub repository in the abstract
4. Submit to math.NA

After arXiv submission, update `CITATION.cff` with the arXiv ID.
