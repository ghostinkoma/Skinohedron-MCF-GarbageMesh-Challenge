# Peer Review & Response

This directory keeps the **adversarial peer review** of the Kosaka Skin-o-hedron
Model and our point-by-point response — including the experiments we ran to
settle each objection, **and any objection that turned out to be correct.**

## Why this is public

Most authors hide negative reviews. We keep them, in the open, on purpose:

1. **Honesty as evidence.** The first draft of this work made a false claim
   (pointwise `O(h²)` on general meshes). We found it ourselves, by running the
   experiment, and retracted it. The same discipline applies here: a criticism is
   answered with code, not rhetoric — and if the criticism is right, we say so.
2. **Reproducibility taken to its limit.** Every objection that can be tested is
   tested by a script in `review/verification/`. A reader can re-run them.
3. **Pre-empting the next reviewer.** If you arrive with the same objection, the
   experiment is already here.

## Ground rule

> **No claim of "rebutted" without a verification script behind it.**

Writing "we refuted this" *without* evidence would repeat the exact mistake of
the first draft (asserting before verifying). So each response below is tagged
with its current status, and a response is only marked **resolved** once its
script exists, runs, and the numbers are in.

---

## The review

The full review text is in [`2026-06_review_01.md`](2026-06_review_01.md).
It raises four criticisms. Summary and **current** status:

| # | Criticism (reviewer's words, paraphrased) | Our read | Status |
|---|---|---|---|
| **1** | The operator *is* the linear FEM Laplacian; the trace-free tensor is decorative, so there is no novelty. | Possibly correct — this is the most serious one. | ⬜ **testing** |
| **2** | The `O(h^3.68)` super-convergence is just icosahedron symmetry, not a property of KSF. | Likely partly correct. | ⬜ **testing** |
| **3** | A pointwise-divergent operator cannot be used for Navier–Stokes / heat (it will blow up to NaN). | Misunderstanding — FEM solves weak forms. | ⬜ to address |
| **4** | Part II (S-parameters) is an empty sketch; machine-precision unitarity is trivial for a lossless toy. | Largely correct. | ⬜ to address |

No status above is "resolved" yet. This file will be updated as each is settled.

---

## How each criticism will be tested

- **Criticism 1 — `c1_tracefree_vs_fem.py`** (the decisive one).
  Assemble the trace-free KSF operator and the plain linear-FEM (cotangent)
  operator on the *same* mesh and compare them numerically. If they agree to
  machine precision, the reviewer is right and we will say so plainly, and
  re-frame the paper's contribution honestly (a unifying viewpoint, not a new
  operator). If they differ, that difference is the paper's actual core, and we
  characterise exactly where and why.

- **Criticism 2 — `c2_random_regular_mesh.py`.**
  The reviewer's challenge is concrete and fair: show the high-order convergence
  on a *shape-regular but symmetry-broken* mesh family (e.g. Lloyd-relaxed
  random Voronoi/Delaunay), not the highly structured icosphere. If the order
  drops toward `O(h²)`, the reviewer is right and Conjecture 6.2 is withdrawn /
  weakened. If it stays near 4, that is a strong rebuttal. (Note: the 3D cube
  test already gave `O(h²)`, i.e. *not* super-convergent, so we expect the
  reviewer is at least partly correct.)

- **Criticism 3.**
  Clarify in the manuscript that the applications are solved in **weak
  (variational) form**, as every production FEM solver does; nodal values of a
  finite-element solution stay bounded, so "blow up to NaN" does not follow from
  first-order pointwise consistency. This is an exposition fix plus, if useful, a
  small demonstration.

- **Criticism 4.**
  Separate Part II from the main paper and label it explicitly as a preliminary
  framework, narrowing the claims so the target surface is smaller and honest.

---

## Files

```
review/
├── README.md                       ← this file
├── 2026-06_review_01.md            ← the review, verbatim
├── RESPONSE_01.md                  ← point-by-point response (updated as tests land)
└── verification/
    ├── c1_tracefree_vs_fem.py      ← Criticism 1 (operator identity?)
    └── c2_random_regular_mesh.py   ← Criticism 2 (super-convergence real?)
```

## Status legend

- ⬜ **testing / to address** — not yet settled; no conclusion claimed.
- ✅ **resolved (rebutted)** — verification script exists and supports the response.
- ⚠️ **resolved (conceded)** — the criticism was correct; the paper is corrected.

_Last updated: 2026-06 — initial structure; verifications pending._
