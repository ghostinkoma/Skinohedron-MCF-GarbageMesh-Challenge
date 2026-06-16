# Review and Answer — a running log

A chronological, honest record of the peer review this project received and how
each point is being addressed. This is a **living document**: each criticism is
answered with a verification script, and the verdict is written **only after the
numbers are in** — including when the criticism turns out to be correct.

The review verbatim: [`2026-06_review_01.md`](2026-06_review_01.md).

> **Discipline.** The first draft of this work failed precisely by *claiming
> before verifying*. We will not repeat that here. Nothing below is marked
> "rebutted" without a script in `verification/` that a reader can re-run. A
> criticism that is right will be marked **conceded**, not buried.

## Status legend
- ⬜ **open** — not yet settled; no conclusion asserted.
- 🔬 **under test** — script being written / run.
- ✅ **rebutted** — verification supports our response.
- ⚠️ **conceded** — the criticism was correct; manuscript corrected accordingly.

---

## Timeline

### 2026-06 · Review #01 received
An adversarial review arrived with verdict *Reject / Major Redirection* and four
criticisms. Rather than rebut rhetorically, we logged them here and committed to
settling the two testable ones with code.

### 2026-06 · Triage and plan
Our honest first read of each criticism, before testing:

| # | Criticism | First read | Status |
|---|---|---|---|
| 1 | Operator = linear FEM; trace-free tensor is decorative ⇒ no novelty | **Possibly correct.** The manuscript itself says Δ_h equals the FEM Laplacian. Must be tested by direct numerical comparison. | ⬜ open |
| 2 | `O(h^3.68)` super-convergence is icosahedron symmetry, not KSF | **Likely partly correct.** Our own 3D cube test already gave plain `O(h²)`, *not* super-convergence — consistent with the reviewer. Needs a symmetry-broken shape-regular mesh test in 2D. | ⬜ open |
| 3 | Pointwise divergence ⇒ unusable for Navier–Stokes (NaN blow-up) | **Misunderstanding.** Production FEM solves *weak* forms; finite-element nodal values stay bounded. But the manuscript never said "weak form," so the fix is exposition. | ⬜ open |
| 4 | Part II is an empty sketch; lossless-toy unitarity is trivial | **Largely correct.** Part II is a framework, not a theorem; machine-precision unitarity for a lossless toy is expected. | ⬜ open |

### 2026-06 · Experiments (pending)
- `verification/c1_tracefree_vs_fem.py` — does the trace-free KSF operator equal
  the linear-FEM cotangent operator on the *same* mesh, to machine precision?
  **This single experiment decides whether the paper has an operator-level
  contribution at all.**
- `verification/c2_random_regular_mesh.py` — measure the eigenvalue convergence
  order on a *shape-regular but symmetry-broken* mesh family (Lloyd-relaxed
  random Delaunay). Does the high order survive, or collapse toward `O(h²)`?

_Results will be filled in here as each script lands._

---

## Point-by-point response (updated as tests land)

### Response to Criticism 1 — "it's just FEM" ⬜
Planned test: assemble `L_tracefree` (KSF) and `L_FEM` (cotangent) on identical
meshes; compare entrywise and compare spectra.
- **If identical:** the reviewer is correct at the operator level. We will
  **concede** and re-frame the contribution honestly — the value, if any, is a
  *unifying viewpoint / assembly formalism*, not a new operator, and the paper
  must say so in the abstract.
- **If different:** characterise exactly where (anisotropy? non-flat metric?
  boundary handling?) — that difference becomes the paper's real core.

*No verdict yet — `c1` not yet run.*

### Response to Criticism 2 — "super-convergence is icosahedron magic" ⬜
Planned test: `c2` on random shape-regular meshes. We note up front that the **3D
cube result in this very repo already shows `O(h²)`**, i.e. no super-convergence
off the icosphere — so we expect the reviewer is at least partly right and that
Conjecture 6.2 will be **weakened or withdrawn**, restricted to highly
symmetric meshes only.

*No verdict yet — `c2` not yet run.*

### Response to Criticism 3 — "pointwise divergence ⇒ NaN" ⬜→ exposition
Substance: a Galerkin finite-element solution minimises energy over a
finite-dimensional space; its nodal coefficients solve a well-posed linear
system and stay bounded. First-order *pointwise* consistency error does not imply
*nodal* blow-up. Every production solver (FEniCS, deal.II, OpenFOAM) runs on the
weak form for exactly this reason. **Fix:** state explicitly in the applications
section that problems are posed in weak form; optionally add a small demo.

### Response to Criticism 4 — "Part II is a poem" ⚠️ (leaning concede)
Agreed in substance. **Fix:** split Part II out of the main paper, label it a
*preliminary framework*, and narrow its claims; do not present machine-precision
unitarity of a lossless toy as a result.

---

_Last updated: 2026-06 — review logged, plan set, experiments c1/c2 pending._
