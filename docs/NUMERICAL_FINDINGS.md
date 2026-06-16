# Numerical Findings

All numbers here are produced by `python3 src/verification/run_all.py`.
Raw output is in `results/results.json`.

---

## §2: The mesh family is well-behaved

The icosphere refinement sequence satisfies the KSF axioms:

| Level | Vertices | Mesh size h | Shape quality q_min | Area error |
|-------|----------|-------------|---------------------|------------|
| 1 | 42 | 0.582 | 0.985 | 7.2% |
| 2 | 162 | 0.299 | 0.977 | 1.9% |
| 3 | 642 | 0.151 | 0.975 | 0.47% |
| 4 | 2,562 | 0.075 | 0.974 | 0.12% |
| 5 | 10,242 | 0.038 | 0.974 | 0.030% |
| 6 | 40,962 | 0.019 | 0.974 | 0.0075% |

- Shape quality stays above **0.974** at all levels → uniform shape-regularity ✅
- Area error decays at order **2.00** in h → geometric approximation is second order ✅

---

## §3: The discrete algebra is exact

The discrete exterior calculus complex closes:

- `d₁ ∘ d₀ = 0` exactly (floating point zero, not just small)
- Euler characteristic `V − E + F = 2` on every level

This means the combinatorial structure is correct.

---

## §4–5: The trace-free projection works perfectly

For Y₂ₓᵧ (a smooth test function on the sphere):

| Level | Trace residual | Frame-invariance residual |
|-------|---------------|--------------------------|
| 2 | 1.08 × 10⁻¹⁹ | 9.76 × 10⁻¹⁸ |
| 3 | 1.19 × 10⁻²⁰ | 2.78 × 10⁻¹⁸ |
| 4 | 7.41 × 10⁻²² | 8.84 × 10⁻¹⁹ |

These are essentially machine precision (≈ 10⁻¹⁶ for float64, but the
cancellations here are exact). The projection is correct and coordinate-invariant.

---

## §6: The central finding — pointwise vs. spectral

This is the most important section. We test the Laplacian against the
exact spherical harmonic Y₂ₓᵧ, whose exact eigenvalue is λ = 6.

**On a regular (icosphere) mesh:**

| h | Pointwise error | Spectral error |
|---|----------------|----------------|
| 0.299 | 4.6% | 0.0046% |
| 0.151 | 2.1% | 0.00041% |
| 0.075 | 1.0% | 0.000032% |
| 0.038 | 0.50% | 0.0000024% |
| 0.019 | 0.25% | 0.00000018% |

- Pointwise error halves when h halves → **first order** (O(h^1.06))
- Spectral error drops by factor ~12 when h halves → **near fourth order** (O(h^3.68))

**On a persistently irregular mesh (jittered so quality never improves):**

| h | Pointwise error | Spectral error |
|---|----------------|----------------|
| 0.317 | 33% | 0.31% |
| 0.165 | 152% | 0.83% |
| 0.083 | 260% | 0.88% |
| 0.041 | 509% | 0.88% |
| 0.021 | 944% | 0.84% |

- Pointwise error **gets worse** as mesh refines → **diverges**
- Spectral error **stagnates** around 0.88% → does not improve

### What this means

The operator is **not** second-order accurate pointwise. The first draft's claim
was too strong.

But for most real applications:
- You care about **eigenvalues** (vibration frequencies, buckling loads) → excellent
- You care about **energy** (heat, electric potential, fluid flow) → good
- You rarely need pointwise accuracy at every single vertex

So the useful result stands. It just needs to be stated honestly.

---

## §7: Non-Delaunay meshes

As we introduce random distortion (jitter) to the mesh:

| Jitter | Worst shape quality | Negative-weight edges | Min Voronoi area |
|--------|--------------------|-----------------------|------------------|
| 0.00 | 0.974 | 0.0% | 0.00379 |
| 0.10 | 0.624 | 0.2% | 0.00346 |
| 0.20 | 0.012 | 5.1% | 0.00247 |
| 0.30 | 0.005 | 13.8% | 0.00171 |
| 0.45 | 0.008 | 20.7% | 0.00082 |

- Negative edge weights → **maximum principle is lost** (heat can flow "uphill")
- But minimum Voronoi area stays **strictly positive** → operator never becomes singular

The operator is **non-degenerate** but **not maximum-principle preserving** on
bad meshes. This is forced by the no-free-lunch theorem: you cannot have both.

---

## §10–12: Part II S-parameter operator

Tested on a lossless, reciprocal scattering kernel on the sphere:

| Level | Unitarity residual | Reciprocity residual | Max gain |
|-------|--------------------|---------------------|----------|
| 2 (162 cells) | 1.06 × 10⁻¹⁵ | 4.31 × 10⁻¹⁷ | 1.000000 |
| 3 (642 cells) | 7.51 × 10⁻¹⁶ | 3.43 × 10⁻¹⁷ | 1.000000 |
| 4 (2562 cells) | 6.90 × 10⁻¹⁶ | 4.03 × 10⁻¹⁷ | 1.000000 |

All three physical invariants satisfied to machine precision.
Part II is implementable as written.

---

## Summary table (paper §A)

| Section | What was checked | Pass/Fail |
|---------|-----------------|-----------|
| §2 | Shape quality bounded, area defect O(h²) | ✅ |
| §3 | d∘d = 0, Euler = 2 | ✅ |
| §4–5 | Trace residual ≲ 1e-19, frame-invariance ≲ 1e-18 | ✅ |
| §6 | Pointwise order 1.06 (not 2), spectral order 3.68 | ⚠️ corrected |
| §7 | Non-degenerate, loses max principle on bad meshes | ⚠️ corrected |
| §10 | Unitarity / reciprocity / passivity to machine eps | ✅ |
