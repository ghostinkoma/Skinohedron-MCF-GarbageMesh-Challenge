# Review & Answer — Continuation (Round 2, 2026-06)

A continuation of [`REVIEW_AND_ANSWER.md`](REVIEW_AND_ANSWER.md). Round 1 settled
the four criticisms of Review #01 (C1, C2 **conceded**; C3 **addressed**; C4
**conceded**) and produced one constructive Tier‑1 result
([`fix01`](verbatim/fix01_tracefree_anisotropic.py)): used as an anisotropic
conductivity, the trace‑free tensor yields an operator that is **distinct from
the isotropic FEM Laplacian**, SPD, and O(h²)‑convergent.

That immediately exposed the two open burdens written into
[`REVISION_PROPOSAL.md`](REVISION_PROPOSAL.md):

- **Tier 3 (the real burden).** Distinct from *isotropic* FEM is not enough.
  Novelty must be proven against **anisotropic** FEM, the correct baseline — or
  conceded as a negative result.
- **Tier 2.1 (geometric meaning).** The `E` in `fix01` was an arbitrary lab‑frame
  field. Can it come from real geometry, and is the resulting operator a
  *consistent* discretisation of a genuine anisotropic surface energy?

Round 2 settles both, with code. **Same discipline as before: nothing is marked
resolved without a script that runs and the numbers in — including when the
answer is "no".**

---

## Outcome summary (Round 2)

| # | Open item | Verdict | Evidence |
|---|---|---|---|
| T3 | Is the anisotropic operator anything beyond anisotropic P1 FEM? | ⚠️ **conceded** (it is *exactly* anisotropic FEM) | [`fix03_aniso_is_fem.py`](verbatim/fix03_aniso_is_fem.py) |
| T2.1a | Can curvature supply the trace‑free anisotropy on S²? | ⚠️ **no** (the sphere is umbilic) | [`fix02_geometric_anisotropy.py`](verbatim/fix02_geometric_anisotropy.py) PART 1 |
| T2.1b | With an external director, is `L(G)` a consistent anisotropic operator? | ✅ **yes** (O(h²) to the exact energy) | [`fix02_geometric_anisotropy.py`](verbatim/fix02_geometric_anisotropy.py) PART 2 |

---

## T3 · The anisotropic operator **is** anisotropic FEM — conceded

`fix03` repeats the decisive `c1` experiment one level up. The same operator is
assembled two genuinely independent ways:

- **Route A** — gradient quadrature with `G` in the middle,
  `K_e = Area · ∇φ · G · ∇φᵀ` (the `fix01` / "KSF" assembly).
- **Route B** — a classical fact about anisotropic P1 FEM: a constant SPD
  conductivity is removed by a linear change of tangent coordinates
  `y = B⁻¹x` with `G|_tangent = B Bᵀ`. In the new coordinates the form is
  **isotropic**, so the element matrix is just `|det B|` times the ordinary
  **cotangent** stiffness of the remapped triangle. Route B never multiplies by
  `G`; it only moves vertices and calls the repo's own isotropic cotan formula.

Tested with three different trace‑free tangential fields `E` (the fixed `fix01`
field, a geometric director field, and a *random* one), at `eps ∈ {0, 0.1, 0.3}`:

| `E` field | worst `relerr(A vs B)` | `relerr(A vs isotropic FEM)` at eps=0.3 |
|---|---|---|
| fix01 fixed | 4.8 × 10⁻¹⁶ | 9.9 × 10⁻² |
| geometric director | 4.9 × 10⁻¹⁶ | 1.2 × 10⁻¹ |
| random trace‑free | 5.2 × 10⁻¹⁶ | 2.0 × 10⁻¹ |

Routes A and B agree to **machine precision for every field**, while both differ
from the isotropic Laplacian by 10–20%. **Conclusion (conceded):** the
anisotropic KSF operator is, identically, the isotropic cotangent (Dziuk P1 FEM)
Laplacian of a per‑face `G`‑remapped mesh — i.e. **textbook anisotropic P1 FEM**.
Just as `c1` showed at the isotropic level, **no new *operator* is produced at the
anisotropic level either.** The revised paper may **not** claim a novel operator.

What is *not* refuted: the layered / trace‑free *formalism* might still be a
convenient, frame‑invariant **bookkeeping** for the anisotropy carrier. But that
is a small, expository claim that must be weighed against finite‑element exterior
calculus and material Hodge stars, which already do this — and it is **not
demonstrated here**, so it stays **open and unclaimed**.

## T2.1 · Geometric anisotropy — a clean negative and a modest positive

`fix02` PART 1 — **curvature cannot supply `E` on S².** The dimensionless shape‑
operator anisotropy `|κ₁−κ₂|/|κ₁+κ₂|`, measured with the *exact* Gauss map of the
unit sphere, is `~10⁻¹⁶` at every level: the sphere is **umbilic** (`κ₁=κ₂=1`),
so the **deviatoric (trace‑free) part of the shape operator is exactly zero**.
(With merely *approximate* normals one sees `~10⁻²`, but that is the estimator's
own discretisation noise — it shrinks under refinement and is not a real
geometric signal.) So on S² the trace‑free tensor's content is an **external
director field**, not an intrinsic curvature invariant. A non‑umbilic surface
(e.g. an ellipsoid) is where a curvature‑derived `E` would become non‑trivial —
a concrete next test, not done here.

`fix02` PART 2 — **with an external director, `L(G)` is consistent.** Take the
tangential projection of the global x‑axis as a director `d`, set
`E = d dᵀ − ½(I − n nᵀ)` and `G = (I − n nᵀ) + eps·E`. Against the exact
continuous energy `a(u) = ∫_S ∇u · G · ∇u` (computed by high‑order spherical
quadrature) for `u = x·y`, the discrete energy `uᵀ L(G) u` converges:

| eps | fitted order | reference `a(u)` |
|---|---|---|
| 0.0 | **O(h²·⁰⁰)** | 5.026548 ( = analytic `8π/5` ✓ ) |
| 0.3 | **O(h²·⁰²)** | 5.152131 (quadrature) |

So `L(G)` is a sound, O(h²)‑consistent discretisation of a genuine anisotropic
surface‑diffusion energy (and reproduces the analytic `8π/5` at `eps=0`).

**Honest scope.** O(h²) consistency is exactly what anisotropic P1 FEM already
delivers — which, by T3, is what this operator *is*. PART 2 therefore confirms
the construction is **correct and usable**; it does **not** show any advantage
over anisotropic FEM. Consistency ≠ novelty.

---

## Updated honest bottom line

After two rounds, the standing of the whole programme — stated without
inflation — is:

> The trace‑free dual tensor, once actually wired into an operator, gives a
> **genuine, well‑posed, SPD, O(h²)‑consistent anisotropic surface‑diffusion
> operator**. That operator is **identically anisotropic P1 FEM** (T3). On the
> sphere its anisotropy must be carried by an **external director**, because the
> surface is umbilic (T2.1a). Nothing here is a new *operator*, a new
> *convergence rate*, or a demonstrated *advantage* over the mature anisotropic‑
> FEM / FEEC literature.

This is a smaller result than Draft v1 or even v2 asserted — and it is the one
the numbers actually support.

## What remains genuinely open (for a Round 3, if pursued)

1. **Non‑umbilic surfaces.** Repeat `fix02` on an ellipsoid / torus, where the
   deviatoric shape operator is non‑zero, and define `E` *from curvature*. Only
   there can "geometric anisotropy" mean something intrinsic.
2. **Any advantage over anisotropic FEM / FEEC.** The Tier‑3 burden is conceded
   as **unmet**, not disproven in general. A concrete problem where the layered /
   trace‑free formalism measurably helps (accuracy, conditioning, assembly cost,
   or solver structure) would be the only thing that revives a novelty claim. If
   honest testing finds none, **that is the publishable result** and should be
   stated as such.
3. **The layered (multi‑shell) structure** is still unused by the operator. Its
   one plausible job — a multigrid hierarchy for solving `L(G)u=b` — remains
   untested (needs a solver + timings).

_Last updated: 2026‑06 (Round 2) — T3 conceded (`fix03`); T2.1a negative,
T2.1b positive (`fix02`). No novelty claim stands; the open programme is items
1–3 above._

---

### Reproduce

```bash
pip install -r requirements.txt
python3 review/verbatim/fix03_aniso_is_fem.py          # T3: operator identity test
python3 review/verbatim/fix02_geometric_anisotropy.py  # T2.1: geometry + O(h^2) energy
```
