# 13 · Surface Tension — a New (5th) Domain the Field Operators Can't Reach

**Status:** **verified**, `src/verification3d/surface_tension_verify.py`. Adds surface
tension as a new domain and tests the intuition that it **cannot** be described by heat or
Coulomb. The answer: the intuition is **correct** — and there is a twist that makes it
beautiful. Surface tension is genuinely distinct from the scalar-field operators, yet it
still uses the **same Laplacian**, applied to the interface **geometry** (the project's
mean-curvature-flow core, `H = Δ_S x`).

---

## 1. Why heat and Coulomb can't reach it

Heat (`02`/`03`) and Coulomb (`06e`) are the **same** machine: a linear PDE for a
**scalar field** `φ` on a **fixed** domain, `K φ = source`. Surface tension is different
in kind on three verified counts:

- **It acts on geometry, not a field.** The curvature comes from the Laplacian applied to
  the interface **position** `X`: `κ n = −Δ X` (mean curvature flow). Verified exact for a
  circle: `κ = 1/R` at `R=1,2,4`. The operator is the same; its *argument* is the shape,
  not a field living on the shape.
- **It makes a pressure jump.** Young–Laplace: the pressure is **discontinuous** across
  the interface, `Δp = σκ` (verified `σ/R` for `R=0.5,1,2`). A smooth `Kφ` field cannot
  produce a jump — this is structurally outside heat/Coulomb.
- **It minimises area (free boundary).** Under `Ẋ = Δ_S X` a non-circular loop's area
  shrinks (`6.32→3.81`) and its isoperimetric ratio rises to a circle (`0.41→0.93`). Heat
  smooths a *field*; this smooths the *shape* — the domain itself evolves.

Net surface-tension force on a closed interface is zero to machine precision (`6.6e-16`):
it drives internal pressure and shape, not translation — a consistency check.

---

## 2. The twist — it is still the same Laplacian, on geometry

The unifying thread does not break; it deepens. `κ n = −Δ X` is the **same Laplacian
operator** that carries heat, waves, pressure, Coulomb and temperature — but applied to
the interface **position** rather than to a scalar field. This is exactly the project's
**mean-curvature-flow** identity `H = Δ_S x` (the repository's namesake). So surface
tension is not outside the one-operator story; it is the one operator **acting on the
geometry instead of on a field**. Heat/Coulomb are `K` on `φ`; surface tension is `Δ` on
`X`. Same `Δ`, different argument — and that difference is precisely what heat and Coulomb
cannot express.

---

## 3. Accuracy — and a caught bug, kept honestly

Curvature converges **`O(h²)`** with the correct operator: RMS error
`1.3e-3 → 3.3e-4 → 8e-5 → 2e-5` for an ellipse, error ratio `~4.00` per refinement. So the
`06f` enrichment framework applies to surface tension too — enrich to sharpen curvature.

**Honest record:** the first attempt used a naive average-`ds` second difference
`(X_{i-1}−2X_i+X_{i+1})/ds²`. It **did not converge** — it plateaued at `~0.099` — because
on a non-uniformly-parametrised curve that stencil carries an `O(1)` parametrisation
error. The fix is the proper non-uniform Laplace–Beltrami stencil
`2/(h_l+h_r)·((X_r−X_i)/h_r − (X_i−X_l)/h_l)`, which restores `O(h²)`. The plateau was
caught by checking convergence against the exact ellipse curvature — the discipline
working as intended.

---

## 4. What is and isn't claimed

**Verified:**
- Curvature `κ n = −Δ X` (exact for circles; `O(h²)` for the ellipse with the correct
  operator).
- Young–Laplace pressure jump `Δp = σκ`.
- Mean curvature flow shrinks area and rounds shapes to circles (area-minimising).
- Net force on a closed interface is zero (machine precision).
- Surface tension shares the Laplacian but acts on geometry — not reducible to heat/Coulomb.

**Not claimed:**
- This is the **surface-tension domain in isolation** (curvature, Young–Laplace, MCF). A
  **fully coupled** two-phase NS solve — surface-tension force `σκ n δ_S` in the momentum
  equation driving a real flow with a moving interface — is **not** done here; that needs
  interface tracking / a free surface (the biggest gap noted in `12 §4`). This note
  establishes the domain and its operator; the coupling is the next build.
- 2-D curves only (curvature as `1/R`); the 3-D mean-curvature surface (`H = κ₁+κ₂`) is
  the same identity but not implemented here.
- No claim that surface tension explains the cavity rev (it doesn't — `10`); it is a new
  domain for free-surface/multiphase physics, where it is essential.

---

## 5. Files

- Theory: this note, on `02`/`03`/`06e`/`12`.
- Verification: `src/verification3d/surface_tension_verify.py` (checks A–E).
- Next build: couple `σκ n δ_S` into the NS momentum equation with a tracked interface
  (free surface) — the pouring-flow goal of `11`/`12 §4`.
