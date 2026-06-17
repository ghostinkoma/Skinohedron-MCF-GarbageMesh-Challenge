# 01c · Plane-wave Interface Test (Step 1.5)

**Status:** theory / experiment design (no code yet). Follows
[`01b_sparameter_model.md`](01b_sparameter_model.md) in the workflow
theory → model → code → viewer. This is **Step 1.5**: it adds *no new physics*.
It is the controlled experiment that promotes the deferred item of
[`sparam_wave_verify.py`](../src/verification3d/sparam_wave_verify.py) —
"quantitative reflected-energy fraction == R²" — from DEFERRED to a measured,
asserted result, so the material model `(R,T)` of `01b §4` is verified as a
*physical observable*, not only as the algebraic identity `R²+T²=1`.

**Why this is worth a step of its own.** Step 2 (pressure) and Step 3 (fluid)
both rely on impedance boundaries behaving correctly. If we do not pin down `(R,T)`
quantitatively now, any later discrepancy could be blamed on the interface and we
would not know. Fixing it here is the "foundation-first" choice.

---

## 1. The problem with a point source (why 3D needs a plane wave)

The coefficients of `01b §4`,

```
R = (Z₂ − Z₁)/(Z₂ + Z₁) ,   T = 2√(Z₁Z₂)/(Z₁ + Z₂) ,                (1)
```

are defined for a wave at **normal incidence** on a planar interface. A point
(monopole) source emits a spherical wave that meets the plane `x = 0` over a
continuum of incidence angles `φ`. The angle-dependent reflectance `R(φ) ≠ R(0)`,
so a spherically-spreading pulse measured against (1) will not match — exactly the
failure seen in the Step-1 suite (measured fraction ≠ R²). The fix is to make the
incident wave **planar and normal** to the interface.

---

## 2. Experiment design: a quasi-1D channel

**Domain.** A rectangular slab `Ω = [−Lx, Lx] × [0, w] × [0, w]`, meshed by the
uniform Kuhn tetrahedra (`kuhn_cube` generalised to a box; same congruent cells).
`w` (transverse width) is kept small — a few cells — because the wave is uniform
there by construction.

**Two media.** `Z(x<0) = Z₁`, `Z(x≥0) = Z₂`. One planar interface at `x = 0`.

**Transverse periodicity.** On the four side walls (`y=0,y=w,z=0,z=w`) impose
**periodic** boundary conditions: a port leaving through `y=w` re-enters at `y=0`,
etc. This removes transverse reflections and forces the only spatial variation to
be along `x`. The result is a genuinely **1-D wave in a 3-D mesh**: a plane wave at
normal incidence. (Periodic pairing is just an extra case of the partner map `Π`
of `01b §1`: side-wall ports are paired with their opposite-wall counterparts
instead of being boundary ports.)

**Ends.** The two `x = ±Lx` ends are **absorbing** (`ρ = 0`, `01b §4 (C-bnd)`), so
each wave packet crosses the interface essentially once and leaves; no end echoes
contaminate the measurement.

```
   absorbing                interface                 absorbing
     end                      x=0                        end
   x=-Lx                       |                        x=+Lx
     |   Z1 (medium 1)         |        Z2 (medium 2)      |
     |        --> incident     |                           |
     |        <-- reflected    |   --> transmitted         |
     |_________________________|___________________________|
        (periodic in y,z : plane wave, normal incidence)
```

---

## 3. What is measured

Launch a plane wave packet in medium 1, travelling +x. Account energy by where it
finally leaves the domain (through the absorbing ends), using the exact per-port
absorbed energy `b_f²` at `ρ=0` ports (the lossless interior conserves energy, so
all injected energy eventually exits one end):

```
E_refl  = Σ_t Σ_{f ∈ left end}  b_f(t)²          (left out  = reflected)
E_trans = Σ_t Σ_{f ∈ right end} b_f(t)²          (right out = transmitted)
E_total = E_refl + E_trans            (= E_injected, lossless interior)         (2)
```

**Predictions (normal incidence, from (1)):**

```
E_refl / E_total  →  R²  ,        E_trans / E_total  →  T²  .          (3)
```

Note: with the `√Z` energy normalisation of `01b §2`, the *power* split is exactly
`R²` and `T²` (this is the energy-correct form; the naive amplitude `T=2Z₂/(Z₁+Z₂)`
would not split energy as `R²:T²`). The test therefore also confirms the `√Z`
normalisation choice.

---

## 4. Acceptance criteria (promotes the DEFERRED item)

The plane-wave interface test passes iff, for each tested pair `(Z₁,Z₂)` ∈
{(1,0.01),(1,0.3),(1,0.5)} and sufficient transverse periodicity:

1. **Energy accounting:** `|E_total − E_injected| / E_injected < 10⁻³`
   (all energy leaves through the ends; interior is lossless).
2. **Reflection match:** `|E_refl/E_total − R²| < ε(h)`, with `ε(h)` the
   discretisation error, shrinking under refinement (test at 2-3 widths/lengths).
3. **Transmission match:** `|E_trans/E_total − T²| < ε(h)`.
4. **Limits:** `Z₂ → Z₁` gives `E_refl/E_total → 0` (no interface); `Z₂ ≪ Z₁`
   (metal) gives `E_refl/E_total → 1` (near-total reflection).

On success, the Step-1 suite's DEFERRED line is replaced by "VERIFIED (01c)".

---

## 5. Honest scope

- This is **not** new physics or a new method: it is the textbook normal-incidence
  reflection/transmission test, run on the tetrahedral TLM engine of `01b` to
  confirm the material model behaves as a physical observable.
- The only modelling addition is **periodic transverse pairing** in `Π`, which is
  a boundary-condition variant, not a new operator. The update `𝒰 = 𝒞𝒮` is
  unchanged.
- A residual discretisation error `ε(h)` is expected and must *decrease* under
  refinement; we report the trend rather than claiming exactness. (Unlike the
  energy identities of Step 1, which are exact, this is a continuum-limit match.)

---

### Roadmap position

```
[Step 1]   S-parameter wave        01, 01b, verify ✓
[Step 1.5] plane-wave interface     01c  ← THIS DOCUMENT (theory)
                                    01c-code  ← next, on approval
[Step 2]   pressure / tank
[Step 3]   fluid dynamics
```

_Next artifact, on approval: the code `src/verification3d/sparam_interface_verify.py`
implementing the periodic-transverse channel, the energy accounting (2), and the
acceptance asserts (4). Reuses the `01b` engine unchanged except for periodic
pairing in `PortComplex`. No viewer for this step — it is a measurement, not a
display._
