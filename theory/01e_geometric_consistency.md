# 01e · Geometric Consistency — Why the Scatter Node Was Rejected

> **Reconstructed (2026-11).** The original `01e` was cited by `01d`, the prologue
> (`00a`) and the theory index but never committed. This note was rewritten from those
> citations; every number below is produced by `src/verification3d/sparam_geometry_verify.py`
> (also reconstructed). Where a cited number differs, the difference is stated.
>
> **Known theory.** A consistent discretisation of the wave equation has an
> anisotropy that vanishes under refinement; a scheme whose anisotropy does not vanish
> is not consistent with the isotropic equation. The FE (cotangent) operator's
> residual `~5%` is its own `O(h²)`-type dispersion on the Kuhn mesh.

**Status:** **verified** (reconstructed), `src/verification3d/sparam_geometry_verify.py`.
The decisive negative result of the S-parameter route (`01`–`01d`).

---

## 1. The question

The scatter node (`01`, `01b`, `ksf3d.sparam_wave`) moves energy between the four
faces of every tetrahedron with the same unitary `S` in every cell and equal impedance
everywhere. It uses the mesh's **connectivity** but not its **geometry** (no face areas,
no distances). Does a wave on it travel at the same speed in every direction?

The measure, as in the prologue:

$$\Delta c / c = \frac{\max_i c(\theta_i) - \min_i c(\theta_i)}{\operatorname{mean}_i c(\theta_i)} .$$

## 2. Scatter node: large, and not shrinking

A monopole pulse is injected at the centre cell of the Kuhn cube; the front speed is
fitted along six lattice directions (distance vs first-arrival step, cells within a
narrow cone around each direction).

| n | `<100>` | `<010>` | `<001>` | `<110>` | `<111>` | `<-111>` | `Δc/c` |
|---|---|---|---|---|---|---|---|
| 8  | 0.199 | 0.201 | 0.226 | 0.177 | 0.289 | 0.214 | **0.517** |
| 12 | 0.194 | 0.203 | 0.207 | 0.191 | 0.314 | 0.218 | **0.555** |
| 16 | 0.192 | 0.195 | 0.195 | 0.192 | 0.324 | 0.208 | **0.606** |

(speeds in mesh widths `h` per step)

- The anisotropy is about **50–60%** and does **not decrease** with refinement — it is a
  property of the node, not of resolution. The cited value "`~50%`, h-independent"
  (and the `~0.49` used as a reference in `unified_scalar_verify.py` /
  `heat_conduction_verify.py`) is consistent with this; the exact figure depends on the
  arrival threshold and cone (stated in the script) and drifts slightly upward with `n`.
- The fast direction is `<111>`, the main diagonal that every Kuhn tetrahedron
  shares: a geometry-blind node inherits the mesh's combinatorial bias.

## 3. FE operator on the same mesh: about ten times more isotropic

Plane-wave phase speed from the Rayleigh quotient of `L = M⁻¹K` (as in `03`):

| n | `c<100>` | `c<110>` | `c<111>` | `Δc/c` |
|---|---|---|---|---|
| 12 | 0.9886 | 0.9520 | 0.9981 | **0.0470** |
| 20 | 0.9959 | 0.9533 | 1.0007 | **0.0482** |

`Δc/c ≈ 0.048`, about **12×** smaller than the scatter node — the "`~10×`" of `00` and `03`.

## 4. What this decided

The scatter node was **rejected** as the wave operator, rather than patched. The route
moved to conduction done with the cotangent / FE conductance (`02`), which became the
one operator `L = M⁻¹K` (`03`). The S-parameter view is kept for what it does well —
losslessness, reciprocity, interfaces (`01`–`01d`) — but not as the propagation operator.

## 5. What is and isn't claimed

**Verified (reconstructed script):** the scatter node's `Δc/c ≈ 0.52–0.61` at n = 8–16,
not shrinking; the FE operator's `Δc/c ≈ 0.048` on the same mesh.

**Not claimed:** that the original 01e used exactly this measurement; that `0.49` is a
precise value (it is a threshold-dependent figure of the same size); any statement about
geometry-aware scatter nodes (with face-area weights), which were not built.

## 6. Files

- Verification: `src/verification3d/sparam_geometry_verify.py` (checks A–B; ~10 s).
- Engines: `src/ksf3d/sparam_wave.py` (scatter node), `src/ksf3d/fem3d.py` (FE operator).
