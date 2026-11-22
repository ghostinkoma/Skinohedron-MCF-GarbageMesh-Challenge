# 14 · Gravity as a Field on the Same Operator — and the Field Sets the Body Force

> **Reading guide.** This note mixes verified numbers with interpretation. See `theory/README.md`, *How to read the claims*, for which results are implementation checks [I], approximation-accuracy results [A], and interpretation [N]. "Coulomb" as a name for the stiffness `K` means *the shared Laplacian*, not that viscosity is electrostatic.

**Status:** **verified**, `src/verification3d/field_coupling_verify.py`. A direct
extension of `06e` (Coulomb): Newtonian gravity is the **same Poisson operator** as
electrostatics, and the gravitational field — like the Coulomb field — enters Navier–Stokes
as a body force. The operator claims (§1, §3) are checked against an exact value or an
invariant with the project's own operator; §2 gives closed-form reference magnitudes only.
Nothing beyond the operator structure is asserted.

> **Correction (issue #3, A2).** An earlier version of this note and its script (a) solved
> §1 with a 1-D radial finite-difference scheme that did not use the project's `K` and
> reported the match as exact (`0.000%`); (b) presented §2's closed-form values as an NS
> result; and (c) computed the §3 superposition error as `‖x − x‖`, which is identically
> zero. All three are corrected below; the numbers are from the corrected script.

---

## 1. Gravity rides the same Poisson operator as Coulomb

Newtonian gravity and electrostatics are the **same `K`** (`K ≈ −∇²`) with different
constants:

- **gravity:** `K Φ = −4πG·b(ρ)`, body force `f = −ρ∇Φ`
- **Coulomb:** `K φ = b(ρ_q)/ε`, body force `f = −ρ_q∇φ`

where `b` is the P1 load vector of the density. The self-gravity potential of a uniform
water ball (`R = 1 m`) is solved with the **project's operator itself** —
`fem_laplacian` on the Kuhn cube, mapped to the box `[−2, 2]³`, with the exact exterior
potential `−GM/r` on the box boundary. The hydrostatic central pressure
`P(0) = ρ(Φ(R) − Φ(0))` is compared with the analytic

`P(0) = (2/3)πGρ²R² = 1.3978×10⁻⁴ Pa`:

| n | `P(0)` [Pa] | error | `‖Φ−Φ_exact‖_∞ / ‖Φ‖_∞` |
|---|---|---|---|
| 8 | `1.3566×10⁻⁴` | `2.95%` | `4.3×10⁻²` |
| 16 | `1.3845×10⁻⁴` | `0.95%` | `2.0×10⁻²` |
| 24 | `1.3890×10⁻⁴` | `0.63%` | `9.5×10⁻³` |
| 32 | `1.3919×10⁻⁴` | `0.42%` | `8.1×10⁻³` |

The error decreases monotonically, at roughly first order — limited by the staircase
approximation of the ball's surface on the cube mesh (the density jump is not aligned
with element faces). It is **not** exact; it converges. So gravity joins the domain map
on the same operator `L = M⁻¹K`.

Same Poisson operator, different coupling constant. This is a classical fact, now explicit
in the project's operator.

---

## 2. The field sets the body-force scale (reference magnitudes, not an NS solve)

A field enters NS as a body force. As an **illustration with closed-form magnitudes**
(no Navier–Stokes or free-surface solve is performed), the same water droplet
(`R = 1 m`) has very different characteristic pressure under three fields:

| field | characteristic pressure | source | symmetry |
|----|------|------|------|
| zero-g (surface tension only) | `0.144 Pa` | Young–Laplace `2σ/R` | isotropic |
| earth-g (uniform, hydrostatic) | `19.6 kPa` | `ρg·2R` | top–bottom asymmetric |
| self-gravity | `1.39×10⁻⁴ Pa` | FE value from §1 (n=32) | spherically symmetric |

In free fall the droplet's interior pressure scale is the surface-tension / self-gravity
value, not the terrestrial hydrostatic one. Whether and how this changes the *flow* is
not computed here; it needs the coupled two-phase NS solve listed as open in `12 §4`.

---

## 3. Several fields ride one operator and superpose

Multiple fields (gravity + Coulomb) ride the same `K` with different sources. On the
Kuhn cube (n=12, `fem_laplacian`), solving the **combined** source once gives the sum of
the separately solved potentials to `3.5×10⁻¹⁶` (relative). Since the body force `−ρ∇φ` is
linear in `φ`, the forces superpose too. This is the expected linearity of the operator
(a consistency check, not a deep result). And the operator is **dimension-
independent**: the same P1 cotangent/gradient operator is exact for a linear field in 2-D
(triangle) and 3-D (tetrahedron), `uᵀKu = 1` for `u = x`. The n-simplex is the minimal
cell whose `n+1` vertices fix an affine field uniquely — the natural carrier for one
scalar field per vertex in any dimension.

---

## 4. What is and isn't claimed

**Verified:**
- gravity = Coulomb = the same Poisson operator `K`; with the project's operator the
  self-gravity central pressure converges to the analytic value (`2.95% → 0.42%`,
  n=8→32);
- potentials from several sources superpose linearly on one operator (`3.5×10⁻¹⁶`);
- the operator is dimension-independent (the n-simplex is the minimal affine cell).

**Not claimed:** that the field changes the *flow* — §2 gives closed-form pressure
magnitudes only, with no NS solve. Nothing beyond this operator structure. No unification of physical
theories, and no statement about the microscopic (quantum) origin of the coefficients —
those are out of scope here.

---

## 5. Files

- Theory: this note, on `06e`/`13`.
- Verification: `src/verification3d/field_coupling_verify.py` (checks A, C, D assert;
  B prints closed-form reference magnitudes only).
