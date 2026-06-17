# 3D Extension — Logical Structure & Verification

This document lays out the **logic** of the 3D Skin-o-hedron extension and records
what is now **verified by code** versus what remains open. It continues the
project's discipline: every verdict below is backed by a script that runs and
whose numbers are quoted; nothing is asserted ahead of the evidence.

The 2D programme ended at a clear, honest place (see
[`../review/CONTINUATION_2026-06_round2.md`](../review/CONTINUATION_2026-06_round2.md)):
the operator is **plain P1 FEM**, accuracy is governed by **mesh quality**
(no-free-lunch), and **no new operator** survived scrutiny. The 3D extension is
built to find out whether the same logic holds one dimension up — and it does.

---

## The three acts

| Act | Object fixed | Script | Status |
|---|---|---|---|
| 1 | the **complex** (skeleton): d0,d1,d2 with `d·d=0`, `χ=1` | [`s3d_complex.py`](../src/verification3d/s3d_complex.py) | ✅ done |
| 2 | the **operator**: P1 FEM Laplacian, convergence vs exact spectra | [`s3d_laplacian.py`](../src/verification3d/s3d_laplacian.py) | ✅ done (this round) |
| 3 | the **S-parameter tetrahedron**: equal-spaced 4-port junction | [`s3d_sparameter.py`](../src/verification3d/s3d_sparameter.py) | ✅ first object (this round) |

Act 1 was the original "opening act". Acts 2 and 3 are the continuation.

---

## Act 1 — the combinatorial substrate (recap)

`solid_ball` fills the unit ball with tetrahedra (icosphere shells + centre,
Delaunay). The discrete de Rham complex closes **exactly** (`d1·d0 = 0` and
`d2·d1 = 0`, machine-zero by the sorted-simplex sign convention) and the Euler
characteristic is `V−E+F−T = 1` at every level — the correct value for a
contractible solid ball, the 3D analogue of the surface's `χ = 2`. The one
blemish, reported honestly, is **tetrahedral quality `q_min = 0`** (Delaunay
slivers between shells). Act 2 shows why that blemish is decisive.

## Act 2 — the operator and the 3D no-free-lunch law

We place the **P1 FEM Laplacian** on the tets and test the lowest Dirichlet
eigenvalue against exact ground truth: `λ₁ = π²` on the ball (the spherical-
Bessel ground state `j_{0,1}²`) and `λ₁ = 3π²` on the cube. Three mesh families,
increasingly "equally spaced":

| Family | quality | fitted order | reading |
|---|---|---|---|
| **A** `solid_ball` (Delaunay shells) | `q_min = 0` | `O(h^−0.7)` | spectrum **structure** correct (1‑3‑5 multiplicities) but **does not converge** |
| **B** `graded_ball` (geometric + rotated shells) | `q_min → 0` | `O(h^0.83)` | converges, but only **first order** — family not uniformly shape‑regular |
| **C** `kuhn_cube` (congruent uniform lattice) | `q ≡ 0.717` | `O(h^2.02)` | **optimal second order** |

The law is the same as 2D, now in 3D: **convergence is bought with mesh quality,
not with the operator.** Slivers stall it; merely-better meshes give first order;
only a **uniformly shape-regular** family (constant `q`, i.e. congruent tets)
recovers the optimal `O(h²)`. This is exactly why **"as equally-spaced as possible
tetrahedra" matters** — it is the precise lever that restores convergence.

Honest verdict (carrying the 2D Tier-3 result up): the operator is **ordinary P1
tetrahedral FEM**. No new operator. What family C demonstrates is the textbook
fact that structured, congruent meshes are second-order accurate — correct and
important, **not novel** over existing FEM.

## Act 3 — the equally-spaced S-parameter tetrahedron

The most equally-spaced placement of 4 ports is the **regular tetrahedron**: its
4 vertices are the unique 4-point set on a sphere maximising the minimum pairwise
angle — every pair subtends the same `arccos(−1/3) = 109.47°`
(verified: `|⟨vᵢ,vⱼ⟩ + 1/3| ≤ 1.1×10⁻¹⁶`). On these ports we build the
**fully tetrahedral-symmetric (T_d) lossless reciprocal junction**

```
S(a,b) = e^{ia} P₀ + e^{ib}(I − P₀),   P₀ = ¼ J  (all-ones projector).
```

Verified to machine precision for arbitrary phases: **unitary** (`‖SᴴS − I‖ ~ 1e-16`),
**reciprocal** (`S = Sᵀ` exactly), **passive** (every `|eigenvalue| = 1`). Its
spectrum splits **1 + 3**: one symmetric mode `e^{ia}` and a 3-fold degenerate
`e^{ib}` — the trivial + standard 3-D irreps of `T_d`. That 3-fold degeneracy is
the **scattering twin of the Laplacian's l=1 triplet** in Act 2: the same
tetrahedral/rotational symmetry showing up in two different operators. A clean
exact corollary falls out: the symmetric lossless tetrahedral junction is **never
reflectionless** (`S_ii = 0` would need `e^{ia} = −3 e^{ib}`, impossible for unit
phases).

This is a finite, exact **symmetry object** — a tidy "S-parameter tetrahedron".
It is not yet a field solver. Two readings of the idea, and how they connect:

1. **Single junction (built here).** A 4-port `T_d`-symmetric S-matrix on the
   maximally equal-spaced ports. Done and verified.
2. **Scattering lattice (next).** Put one such node on **every cell of the
   uniform Kuhn/BCC tetrahedral lattice** (Act 2, family C), with ports on the 4
   shared faces — a TLM-like (transmission-line-matrix) isotropic 3D scattering
   network. Because the cells are congruent and equally spaced, **every node is
   identical and the lattice is isotropic by construction**. This is the natural
   fusion of "equally-spaced tetrahedra" + "S-parameter", and it is the concrete
   next build.

---

## The honest bottom line for 3D

- **Substrate:** correct (`d·d=0`, `χ=1`). ✅
- **Operator:** plain P1 FEM; convergence is **quality-limited** exactly as in 2D;
  uniform tetrahedra restore `O(h²)`. Correct, but **no new operator**. ⚠️
- **S-parameter tetrahedron:** a real, exact, symmetric lossless 4-port; its
  degeneracy mirrors the Laplacian's. A clean object, **not yet a solver**. ✅(object)

So the 3D extension is **not hopeless** — it is sound, it works, and the
"equally-spaced tetrahedra" instinct is mathematically the right one (it is what
buys back convergence). What it is **not** is a source of novelty over FEM / TLM /
finite-element exterior calculus. That remains to be demonstrated, not assumed —
same burden as 2D.

## Open constructions (a Round-3 menu)

1. **Boundary-conforming uniform ball mesh.** Family C is uniform but on a cube;
   `graded_ball` conforms to the sphere but is not uniformly regular. The missing
   piece is a *uniform interior + exact sphere boundary at once*
   (isosurface-stuffing / BCC-snap, à la Labelle–Shewchuk). Then re-run Act 2 on
   the ball and confirm `O(h²)` against the spherical-Bessel spectrum directly.
2. **The assembled scattering lattice.** Build reading (2) above, then verify
   global energy conservation across the lattice and the numerical **dispersion
   relation** against the exact wave speed. This is where any content beyond
   standard TLM would have to appear.
3. **Tie the anisotropy to 3D geometry.** The 2D `fix02` showed the sphere is
   umbilic (no trace-free curvature); a 3D volume mesh of a non-spherical body is
   where a curvature/material-director anisotropy `G` could finally be
   non-trivial. Test the anisotropic 3D FEM operator there.

### Reproduce

```bash
pip install -r requirements.txt
python3 src/verification3d/s3d_complex.py      # Act 1: substrate (d·d=0, χ=1)
python3 src/verification3d/s3d_laplacian.py    # Act 2: operator + mesh-quality law
python3 src/verification3d/s3d_sparameter.py   # Act 3: equal-spaced S-parameter tetrahedron
```

_Last updated: 2026-06 — Acts 2 & 3 added; operator = P1 FEM (no new operator),
convergence quality-limited (uniform tets give O(h²)), tetrahedral 4-port built
and verified. Open items 1–3 above._
