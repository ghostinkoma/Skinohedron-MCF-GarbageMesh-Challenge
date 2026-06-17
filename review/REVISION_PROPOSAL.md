# Revision Proposal — how to make the substrate stand

Both operator-level criticisms were **conceded** (see
[`REVIEW_AND_ANSWER.md`](REVIEW_AND_ANSWER.md)):

- **C1**: the tested operator is exactly the linear FEM (cotangent) Laplacian;
  the trace-free tensor did not enter it. *Verified:*
  [`verification/c1_tracefree_vs_fem.py`](verification/c1_tracefree_vs_fem.py).
- **C2**: the `O(h^3.7)` super-convergence is icosahedral symmetry, not KSF;
  breaking symmetry at fixed shape-regularity drops the order toward `O(h^2)`.
  *Verified:*
  [`verification/c2_random_regular_mesh.py`](verification/c2_random_regular_mesh.py).

Conceding is not the end. This document proposes how the substrate could be
rebuilt so that it actually has content — and is **strict** about separating what
is already verified from what is still only a hypothesis. Asserting the latter as
the former is exactly the mistake that started this whole saga; we do not repeat
it.

---

## The single root cause

Both criticisms reduce to one defect:

> **The trace-free dual tensor was never wired into the operator. It was a
> standalone object verified on the side, while the operator that did the work
> was plain isotropic FEM.**

So there is exactly one honest way forward: **make the trace-free tensor act
inside the operator, and produce something FEM does not.**

---

## Tier 1 — what is VERIFIED now (`fix01_tracefree_anisotropic.py`)

Use the trace-free tensor `E` as an **anisotropic conductivity** in the operator:

```
G = (I - n nᵀ)  +  ε · E ,      tr E = 0 ,   E = Eᵀ ,   E tangential
```

assembled into a P1 anisotropic stiffness `L(G)`. Measured facts:

| Property | Result |
|---|---|
| (A) `ε = 0` reproduces cotangent FEM | rel. diff ~1e-16 (validates assembly) |
| (B) `ε > 0` is **distinct** from FEM | 9% at ε=0.1, 27% at ε=0.3 (entrywise + spectrum) |
| (C) the `ε > 0` operator is **well-posed** | SPD; eigenvalue converges to a definite limit at **O(h²)** |

**Conclusion (verified):** the trace-free tensor *need not be decorative*. There
is a concrete, working construction in which it changes the operator into a
genuine, distinct, SPD, convergent anisotropic operator.

**This is the only thing Tier 1 claims.** It does **not** claim usefulness,
physical meaning, or any advantage over existing methods.

---

## Tier 2 — PLAUSIBLE, not yet shown

These are reasonable next steps, explicitly **unproven**:

1. **Tie `E` to real geometry.** Here `E` was an arbitrary tangential trace-free
   field. To be meaningful it should encode something — surface curvature
   (shape operator), a material director field, or a fibre orientation. Then
   `L(G)` would model **anisotropic diffusion on a surface**, a real and useful
   PDE class.
   *Needed:* a definition mapping geometry → `E`, and a convergence test against
   a problem with a known anisotropic answer.

2. **A convergence theorem for the anisotropic operator** under uniform
   shape-regularity (a Strang-type / Babuška–Osborn argument for the
   `G`-weighted bilinear form). Tier 1 shows O(h²) numerically for one `E`; a
   proof would generalise it.

3. **The layered (multi-shell) structure** of the Skin-o-hedron — so far unused
   by the operator — might earn its keep as a multilevel / multigrid hierarchy
   for solving `L(G) u = b` efficiently. *Needed:* an actual solver + timing.

---

## Tier 3 — the BURDEN the revised paper must accept

If KSF is repositioned around the anisotropic operator, then **novelty must be
proven against the right baseline**:

- Not against the isotropic Laplacian (that comparison is lost — C1).
- But against **anisotropic FEM**, which is a mature field (anisotropic
  diffusion, DEC with material Hodge stars, finite-element exterior calculus).

The revised paper may only claim novelty for whatever the layered mesh +
trace-free dual-tensor formalism adds *on top of* anisotropic FEM — and must
**demonstrate** it, not assert it. If, after honest testing, it adds nothing,
that is itself a publishable negative result and should be stated plainly.

---

## Concrete rewrite checklist for the manuscript

1. **Abstract / title:** drop "new discrete Laplacian." Reposition as "a
   frame-invariant assembly of *anisotropic* operators on layered polyhedral
   complexes, with the trace-free dual tensor as the anisotropy carrier."
2. **Theorem 5.x / Prop 6.x:** delete the isotropic-operator and pointwise
   claims (conceded). Replace with the Tier-1 *verified* statement (distinct,
   SPD, O(h²)-convergent anisotropic operator) — clearly labelled empirical.
3. **Conjecture 6.2 (super-convergence):** withdraw as general; restate only as
   "symmetric meshes superconverge" (textbook), with `c2` as evidence.
4. **§7 (non-Delaunay):** keep, reworded as already agreed ("non-degenerate but
   not maximum-principle preserving").
5. **Applications / Part II:** move out of the main paper; mark preliminary;
   note that PDEs are solved in **weak form** (this answers C3 — nodal values of
   a Galerkin solution stay bounded, so first-order pointwise error does not
   cause NaN blow-up).
6. **Add** `fix01` and the honest Tier-2/Tier-3 open problems as the paper's
   actual research programme.

---

## What a reader arriving at the same idea should take away

If you independently invented a "trace-free dual-tensor Laplacian" and are
excited: the isotropic version is just FEM (C1), and any super-convergence you
see on a pretty mesh is symmetry, not your construction (C2). The **live**
question — open as of this writing — is whether a *geometrically meaningful*
trace-free anisotropy gives an operator that is both distinct from and better
than anisotropic FEM for some real problem. Tier 1 shows the operator is at
least real and well-posed. The rest is unclaimed. Good luck — and please add
your result here, positive or negative.

_Last updated: 2026-06 — Tier 1 verified (`fix01`); Tiers 2–3 open._
