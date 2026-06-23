# 11 · Two Forward Markers — Pour-Vortices (GLSL) and "Does Precision Really Improve?"

**Status:** viewers + note (a *forward marker* and a *summary*, explicitly qualitative
where it is qualitative). Two visual artifacts that set up later chapters and gather the
work so far. Both follow the project rule: any *number* shown is verified solver output;
where a viewer is **illustrative** (the GLSL pour), it is labelled as such and grounded in
already-verified physics, not presented as a solved field.

---

## 1. Pour-vortices — cylindrical cup vs square *masu* (`viewer/viewer_pour.html`)

A GLSL (WebGL) viewer of the secondary vortices that appear when water is poured into a
**round cup** versus a **square masu**, shown with the governing equations. It is a
**qualitative** illustration grounded in the verified result of `10`:

- **Equations on screen:** steady incompressible NS `(u·∇)u = −∇p/ρ + ν∇²u`, `∇·u = 0`,
  `ν = UL/Re`; vorticity `ω = ∂v/∂x − ∂u/∂y`; the Moffatt corner-eddy form
  `ψ ∝ r·e^{−λθ}` (counter-rotating nested eddies in a corner).
- **The physical point (verified in `10`):** the **round cup has no corners**, so no
  Moffatt corner eddies form — only the smooth central primary vortex. The **square masu
  has four corners**, so counter-rotating Moffatt eddies appear there, and they
  **strengthen with Re** (exactly the `α`-scaling mechanism of `10 §3`: corner vortices
  grow as nonlinear advection turns on). The viewer encodes this: `shape=cup` zeroes the
  corner term; `shape=masu` adds decaying corner oscillations whose amplitude rises with
  Re.

This is a *forward marker* for a later chapter that will solve the actual pouring flow on
the operator; for now it visualises, with the right equations, the corner-vortex physics
already established. **It is not a solved NS field** — the stream function is an
illustrative encoding, labelled qualitative.

---

## 2. "Does precision really improve?" (`viewer/viewer_precision.html`)

A data viewer — **verified numbers only** — gathering the accuracy story as a summary of
the rev investigation and the enrichment route. Five panels:

- **Enrichment P1→P8** (`06f`): the field-intrinsic `(3)` sinusoid error falls spectrally
  `2.3e-2 → 2.35e-13` by raising local order. *Yes, high-order really reaches machine
  precision.*
- **Refinement `O(h²)`** (`06d`/`06f`): the approximation `(1)` error falls `5.0→1.3%`
  with mesh; mass `M` is exact at every `n`.
- **Real benchmark** (`08b`/`10`): Ghia `Re=1000` centreline RMS falls under both
  refinement (PSPG) and enrichment (Taylor–Hood, more accurate per node).
- **Rev separation** (`10`): the numerical rev shrinks (`(1)`), the physical corner vortex
  converges to a fixed value (`(3)`) — the two plotted together.
- **Summary:** precision *does* improve with `(1)`+`(3)`; what improves is the
  approximation; the field-intrinsic structure itself remains (correctly).

This answers the question the chapter title poses — *if we do `(1)`, `(3)`, and
enrichment, does accuracy really rise?* — with the verified curves: **yes**, and it makes
explicit *which* part rises (approximation) and which part is real physics that should
**not** vanish (the field-intrinsic structure).

---

## 3. What is and isn't claimed

**Verified (the precision viewer):** every number is from `06d`/`06f`/`08b`/`10`; the
viewer renders presets, no in-browser physics.

**Illustrative (the pour viewer):** the pour vortices are a **qualitative** GLSL encoding
of verified corner-eddy physics, with correct equations — **not** a solved NS field. The
later chapter that solves the actual pouring flow on the operator is the forward marker
this sets up.

No new mathematics; this is a consolidation-and-setup step.

---

## 4. Files

- Viewers: `viewer/viewer_pour.html` (GLSL, illustrative), `viewer/viewer_precision.html`
  + `viewer/data_precision.js` (verified data).
- Theory: this note, on `10`/`06f`/`08b`.
- Forward marker: a later chapter solving the real pouring/free-surface-ish flow on the
  operator (a larger build — free surface / inflow BCs).
