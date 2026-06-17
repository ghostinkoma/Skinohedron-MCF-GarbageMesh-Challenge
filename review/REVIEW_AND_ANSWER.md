# Review and Answer — a running log

A chronological, honest record of the peer review this project received and how
each point was addressed. This is a **living document**: each criticism is
answered with a verification script, and the verdict is written **only after the
numbers are in** — including when the criticism turns out to be correct.

The review verbatim: [`2026-06_review_01.md`](2026-06_review_01.md).
The forward path: [`REVISION_PROPOSAL.md`](REVISION_PROPOSAL.md).

> **Discipline.** The first draft of this work failed precisely by *claiming
> before verifying*. We do not repeat that here. Nothing is marked "rebutted"
> without a script in `verification/` that a reader can re-run. A criticism that
> is right is marked **conceded**, with the evidence, not buried.

## Status legend
- ⬜ **open** — not yet settled; no conclusion asserted.
- ✅ **rebutted** — verification supports our response.
- ⚠️ **conceded** — the criticism was correct; manuscript corrected accordingly.
- 🛠 **addressed** — exposition / structural fix (no numerical dispute).

---

## Outcome summary

| # | Criticism | Verdict | Evidence |
|---|---|---|---|
| 1 | Operator = linear FEM; trace-free tensor is decorative | ⚠️ **conceded** | [`c1_tracefree_vs_fem.py`](verification/c1_tracefree_vs_fem.py) |
| 2 | `O(h^3.7)` super-convergence is icosahedron symmetry | ⚠️ **conceded** | [`c2_random_regular_mesh.py`](verification/c2_random_regular_mesh.py) |
| 3 | Pointwise divergence ⇒ unusable (NaN blow-up) | 🛠 **addressed** | weak-form clarification (below) |
| 4 | Part II is an empty sketch | ⚠️ **conceded** | restructure (below) |

**And a constructive result:** a verified way the trace-free tensor *can* matter —
[`fix01_tracefree_anisotropic.py`](verification/fix01_tracefree_anisotropic.py) —
written up in [`REVISION_PROPOSAL.md`](REVISION_PROPOSAL.md).

---

## Timeline

### 2026-06 · Review #01 received
Verdict *Reject / Major Redirection*, four criticisms. We logged them and
committed to settling the testable ones with code.

### 2026-06 · Criticism 1 tested → CONCEDED
`c1_tracefree_vs_fem.py`: the cotangent Laplacian (the operator the paper's
convergence test actually uses) equals an independently-assembled linear P1 FEM
stiffness matrix to **rel. diff ~1e-16** — bit-for-bit, across levels. The
trace-free tensor does not enter that operator. **The reviewer is right:** at the
operator level the paper rediscovered the classical Dziuk FEM Laplacian. No new
*operator* may be claimed.

### 2026-06 · Criticism 2 tested → CONCEDED
`c2_random_regular_mesh.py`: on the perfect icosphere the eigenvalue error
converges at ~`h^3.65`. Breaking the icosahedral symmetry while keeping the mesh
**shape-regular** (q_min ≈ 0.32, bounded away from 0) drops the order to
~`h^2.77` (mean over 5 random seeds, ±0.03). **The reviewer is right:** the
super-convergence is the classical symmetry-cancellation effect of structured
meshes, not a property of KSF. *Why:* on a symmetric mesh the leading O(h²)
truncation terms lie in directions related by the icosahedral group and largely
cancel in the symmetric average, exposing the higher-order residual; breaking
the symmetry removes the cancellation and the generic O(h²) reappears.
**Conjecture 6.2 is withdrawn** as a general claim.

### 2026-06 · Criticism 3 → ADDRESSED (exposition)
Substance: a Galerkin finite-element solution solves a well-posed linear system;
its **nodal values stay bounded**. First-order *pointwise* consistency error does
not imply nodal blow-up, and every production solver (FEniCS, deal.II, OpenFOAM)
runs on the **weak form** for exactly this reason. The reviewer's "NaN blow-up"
does not follow. **Fix:** state explicitly in the applications section that
problems are posed in weak form. (No numerical dispute; this is a wording defect
in the manuscript, now flagged.)

### 2026-06 · Criticism 4 → CONCEDED (restructure)
Agreed: Part II is a framework, not a theorem, and machine-precision unitarity of
a lossless toy is expected, not a discovery. **Fix:** split Part II out of the
main paper, label it preliminary, narrow its claims.

### 2026-06 · Constructive response → REVISION PROPOSAL
Conceding C1/C2 raised the obvious question: can the trace-free tensor be made to
matter at all? `fix01_tracefree_anisotropic.py` shows that, used as an
**anisotropic conductivity** `G = (I - n nᵀ) + ε E`, it yields an operator that
(A) reduces to FEM at ε=0, (B) is genuinely **distinct** from FEM for ε>0
(9–27%), and (C) is **SPD and O(h²)-convergent** to its own limit. So the tensor
*need not* be decorative — there is a concrete construction in which it changes
the operator. **That is all that is claimed.** Usefulness, physical meaning, and
any advantage over the mature field of anisotropic FEM remain **open** and must
be proven, not asserted ([`REVISION_PROPOSAL.md`](REVISION_PROPOSAL.md), Tiers 2–3).

---

## Honest bottom line

The original operator-level claims did not survive contact with their own
verification scripts. We conceded both, on the record, with reproducible
evidence. What remains is **not** the grand theory of the first draft, but
something smaller and real: a frame-invariant way to assemble *anisotropic*
operators on layered meshes, whose value (if any) is now an explicit, testable
research programme rather than an assertion.

_Last updated: 2026-06 — C1 & C2 conceded with evidence; C3 addressed; C4
restructured; fix01 verified; Tiers 2–3 open._
