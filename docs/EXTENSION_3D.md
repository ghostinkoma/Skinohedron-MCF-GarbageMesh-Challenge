# 3D Extension — Step 1: the combinatorial substrate

This is the **opening act** of the 3D Skin-o-hedron extension. It fixes the
complex (the skeleton), not yet the operator (the Laplacian).

## What is verified here

We fill the solid unit ball with tetrahedra (concentric icosphere shells +
centre, Delaunay-tetrahedralised) and check, at every refinement level:

| Property | 2D (surface) | 3D (solid ball) | Result |
|---|---|---|---|
| Complex closes (1st) | `d1·d0 = 0` | `d1·d0 = 0` | ✅ exactly 0 |
| Complex closes (2nd) | — | `d2·d1 = 0` | ✅ exactly 0 |
| Euler characteristic | `V−E+F = 2` | `V−E+F−T = 1` | ✅ at every level |
| Boundary volume | area → 4π | volume → 4π/3 | ✅ as surface refines |

The jump from 2D to 3D adds a **second exterior derivative** (`d2`) and a
**third simplex type** (tetrahedra). Getting `d2·d1 = 0` and `χ = 1` exactly is
the minimal evidence that the 3D substrate is wired correctly.

`χ = 1` is the correct value for a solid ball (it is contractible), just as
`χ = 2` was correct for the sphere *surface*.

## Why `d·d = 0` is exact (not just small)

The boundary operators are built with the **sorted-simplex, alternating-sign**
convention:

```
∂[v0,…,vp] = Σ_i (−1)^i [v0,…,v̂i,…,vp]    (vertices sorted)
```

With this convention `∂∂ = 0` holds in integer arithmetic, so the coboundary
compositions `d1·d0` and `d2·d1` are machine-zero by construction. Verifying it
numerically still matters: it confirms the edge/face/tet enumeration and
indexing are correct.

## Honest limitations

- **Tetrahedral quality is modest.** Delaunay tetrahedralisation of points lying
  on concentric shells produces some *sliver* tets (nearly flat), so the
  radius-ratio quality `q_min ≈ 0`. This is reported honestly by `tet_quality()`.
- This does **not** affect the combinatorial results above (`d·d = 0`, `χ`),
  which are purely topological.
- It **will** matter for the next step (the Laplacian), and that is precisely
  where the Wardetzky et al. *"no free lunch"* tension is expected to reappear
  in 3D — exactly as it did on irregular surface meshes in §6–§7 of the 2D work.

## How to run

```bash
python3 src/verification3d/s3d_complex.py
```

## Files

```
src/ksf3d/
  ├── mesh3d.py   ← solid-ball tetrahedral mesh + geometry diagnostics
  └── dec3d.py    ← 3D coboundary operators d0, d1, d2 (∂∂ = 0 by construction)
src/verification3d/
  └── s3d_complex.py   ← this verification, with conclusion (結句)
```

## Next step (not yet done)

Place a 3D Laplacian on this complex and test its eigenvalues against the
**exact spherical-Bessel eigenvalues** of the Dirichlet Laplacian on the ball
(`−Δ u = λ u`, `u|∂B = 0`, with `λ = j_{l,k}²` where `j_{l,k}` are zeros of the
spherical Bessel functions). That is where the real convergence question lives.
