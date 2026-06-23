# 14 · Gravity as a Field on the Same Operator — and the Field Changes NS

**Status:** **verified**, `src/verification3d/field_coupling_verify.py`. A direct
extension of `06e` (Coulomb): Newtonian gravity is the **same Poisson operator** as
electrostatics, and the gravitational field — like the Coulomb field — enters Navier–Stokes
as a body force, so it changes the flow and the pressure. Every claim here is checked
against an exact value or an invariant; nothing beyond the operator structure is asserted.

---

## 1. Gravity rides the same Poisson operator as Coulomb

Newtonian gravity `K Φ = 4πG·Mρ` and electrostatics `K φ = M ρ_q/ε` are the **same `K`**
with different constants. Solving the self-gravity potential of a uniform ball with this
operator and integrating `ρ·g(r)` gives the central pressure

`P(0) = (2/3)πGρ²R²`,

matched **exactly** (`0.000%`). So gravity joins the domain map on the same operator
`L = M⁻¹K`:

- **gravity:** `K Φ = 4πG·Mρ`, body force `f = −ρ∇Φ`
- **Coulomb:** `K φ = M ρ_q/ε`, body force `f = −ρ_q∇φ`

Same Poisson operator, different coupling constant. This is a classical fact, now explicit
in the project's operator.

---

## 2. The field changes the NS behaviour

A field enters NS as a body force, so it sets the flow and the pressure. The same water
droplet (`R = 1 m`) has totally different pressure under three fields:

| field | central / characteristic pressure | symmetry |
|----|------|------|
| zero-g (surface tension only) | `0.144 Pa` | isotropic |
| earth-g (uniform, hydrostatic) | `19.6 kPa` | top–bottom asymmetric |
| self-gravity | `1.4×10⁻⁴ Pa` | spherically symmetric |

In free fall the droplet's interior pressure is the surface-tension / self-gravity value,
not the terrestrial hydrostatic one, and the NS behaviour differs accordingly — the
"field parameter" sets the body force, hence the flow.

---

## 3. Several fields ride one operator and superpose

Multiple fields (gravity + Coulomb) ride the same `K` with different sources; their NS
body forces **superpose linearly** (error `0`). And the operator is **dimension-
independent**: the same P1 cotangent/gradient operator is exact for a linear field in 2-D
(triangle) and 3-D (tetrahedron), `uᵀKu = 1` for `u = x`. The n-simplex is the minimal
cell whose `n+1` vertices fix an affine field uniquely — the natural carrier for one
scalar field per vertex in any dimension.

---

## 4. What is and isn't claimed

**Verified:**
- gravity = Coulomb = the same Poisson operator `K` (self-gravity central pressure exact);
- the field sets the NS body force, so it changes the flow and pressure (water droplet);
- multiple fields ride one operator and superpose linearly;
- the operator is dimension-independent (the n-simplex is the minimal affine cell).

**Not claimed:** nothing beyond this operator structure. No unification of physical
theories, and no statement about the microscopic (quantum) origin of the coefficients —
those are out of scope here.

---

## 5. Files

- Theory: this note, on `06e`/`13`.
- Verification: `src/verification3d/field_coupling_verify.py` (checks A–D; operator
  structure only).
