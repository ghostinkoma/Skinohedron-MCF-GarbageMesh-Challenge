# 08b · Higher Reynolds — Re=400 & 1000, Both Elements, and the "Rev"

**Status:** **verified**, `src/verification3d/cavity_highre_verify.py`. Pushes the
cross-validated cavity (`08`, `08a`) to `Re=400` and `Re=1000` with **both** the PSPG and
the Taylor–Hood elements, and searches for the **"rev"** — the Reynolds-driven change in
flow structure, and the point where the discretisation begins to strain. Steady laminar
throughout (true unsteady transition is `~Re 8000`, far beyond here).

---

## 1. Both elements vs Ghia (1982)

| Re | PSPG (n=64) centreline RMS | Taylor–Hood (n=32) RMS |
|----|------|------|
| 400 | 0.0311 | 0.0139 |
| 1000 | 0.0555 | 0.0265 |

Both converge to the published benchmark. Taylor–Hood is markedly more accurate per node
(its P2 velocity is the `06f` enrichment). `Re=1000` is genuinely harder — thinner
boundary layers — so the error is larger but still shrinking under refinement, reported
honestly rather than hidden.

---

## 2. The cross-check holds, and loosens — the numerical "rev"

PSPG (with `τ`) and Taylor–Hood (no parameter) still agree at higher Re, but the
agreement **loosens as Re rises**:

| Re | PSPG vs Taylor–Hood mutual RMS |
|----|------|
| 400 | 0.0062 |
| 1000 | 0.0161 |

Two independent discretisations agreeing is the evidence the physics is right (`08a`).
The *loosening* (`0.006 → 0.016`) is the **numerical rev**: at `Re=1000` the schemes begin
to strain against the boundary layers, and they start to part company — a signal, visible
in the data, that more resolution or higher order is needed before pushing further.

---

## 3. The physical "rev" — corner vortices and centre migration

The Reynolds-driven structural change is real and captured by both elements:

| Re | bottom-right reverse flow `u_min` | bottom-left `u_max` | primary-vortex centre |
|----|------|------|------|
| 400 | −0.0705 | +0.0003 | (0.56, 0.60) |
| 1000 | −0.1413 | +0.0030 | (0.53, 0.56) |

- **Secondary corner vortices grow with Re:** the bottom-right reverse flow roughly
  **doubles** (`−0.07 → −0.14`) and the bottom-left vortex strengthens ~`10×`. This is the
  emergence of the well-known Ghia secondary vortices — weak at `Re=400`, clear at
  `Re=1000`.
- **The primary vortex centre migrates toward the geometric centre:** `(0.56,0.60) →
  (0.53,0.56)`, matching Ghia's published centres `(0.5547,0.6055)` and `(0.5313,0.5625)`
  to two digits. This migration *is* the canonical cavity "rev".

So the "rev" the search was for appears on two levels: the **flow** reorganises (corner
vortices grow, centre moves) and the **discretisation** strains (method agreement
loosens). Both are measured, both are consistent with the literature.

---

## 4. What is and isn't claimed

**Verified:**
- Both PSPG and Taylor–Hood reproduce Ghia at `Re=400` and `Re=1000` (converging; TH more
  accurate per node).
- The two methods still agree (mutual RMS `0.006`/`0.016`), with agreement loosening as Re
  rises — the discretisation's strain made visible.
- The physical rev is captured: secondary corner vortices grow with Re and the primary
  centre migrates to the Ghia-published locations.

**Not claimed:**
- **Still steady laminar.** No unsteady/turbulent transition (that is `~Re 8000`, and
  would need time integration and far finer meshes). The "rev" here is *within* the steady
  regime plus the numerical strain — not a Hopf bifurcation to time-dependence.
- `Re=1000` is **not** fully mesh-converged at these resolutions (RMS `~0.03` for TH);
  the trend is convergent but a finer mesh / higher order would tighten it.
- Corner-vortex strength is measured by a simple reverse-flow proxy, not a full
  stream-function topology; the trend is robust, the absolute numbers are indicative.
- 2-D only; no new mathematics — Ghia, PSPG, Taylor–Hood are classical. The contribution
  is the dual-element validation and the explicit, measured account of the rev.

---

## 5. Files

- Theory: this note, on `08`/`08a`.
- Verification: `src/verification3d/cavity_highre_verify.py` (checks A–C; imports both
  element solvers).
- Next rung: a true time-dependent run toward unsteadiness (much higher Re, or a
  deliberately destabilised setup), and/or mesh-convergence at `Re=1000` to tighten the
  numbers — both larger efforts (see `07a`).
