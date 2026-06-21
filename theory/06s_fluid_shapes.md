# 06s · Fluid on Shapes — does the V3 fluid solve on cube, torus, and sphere?

**Status:** theory + **verification (PASS)**. Backed by
`src/verification3d/fluid_shapes_verify.py` (checks A–D, measured inline). Before a fluid viewer shows the three
solids, this chapter establishes — in keeping with 空論を重ねない — *what* fluid
problem each shape supports and *what exact target or invariant* verifies it. It is
the fluid analogue of `05` (which carried heat/wave/pressure onto shapes) and uses
the Stage A/B operators (`06a`, `06b`) and the verified projection (`04`).

The honest headline, from the probes below: **incompressibility is preserved to
machine precision on all three shapes — even the degraded sphere — because the
projection is a linear, exact constraint; mesh degradation costs accuracy in the
dynamics, not the divergence-free property.**

---

## 0. Why shapes differ for a fluid (more than for heat)

For heat and waves (`05`), a shape change was a change of operator `L`. For a
**fluid** three extra things change with shape:

1. **What is a wall.** The cube's Couette/Poiseuille assume plane walls. A torus is
   naturally *periodic* (a ring, no wall in the flow direction); a solid sphere has
   only an outer surface. So each shape gets the flow that is natural to it, not a
   forced channel.
2. **Mesh degradation bites harder.** The sign-flipped cotangents of `5b` (geometric
   ball: ~2844; geometric torus: ~2726) affect advection + projection more than a
   single heat solve. Stability, not just accuracy, must be checked.
3. **Periodicity of advection.** A torus is the *most natural* stage for advection —
   a scalar is carried around the ring and returns. A sphere carries a scalar by
   rotation. The cube needs an artificial periodic direction.

So the model assigns each shape its natural flow and its own verifiable target.

---

## 1. The flow and the target, per shape

| shape | natural flow (prescribed `u`) | what is verified | yardstick |
|---|---|---|---|
| **cube** (reference) | uniform `u=(1,0,0)`, x-periodic | mass conservation, projection `∇·u=0`, variance `2κt` | machine / FE |
| **torus** (x-periodic ring) | uniform around the ring | one-lap **mass conservation**, skew energy, variance | machine / FE |
| **sphere** (solid ball) | solid-body rotation `u = ω×r` (divergence-free) | projection `∇·u=0` on the **degraded** mesh, bounded energy (**graceful**) | machine / bounded |

Two properties are required of *every* shape (shape-independent):

- **Incompressibility.** `‖D u‖ → ` machine precision after the verified projection,
  because the constraint is linear and solved exactly — *independent of mesh
  quality*. (Probe: degraded ball gives `1.2e-14`.)
- **Mass conservation.** The skew-symmetric advection conserves `∫φ = 1ᵀMφ` to
  machine precision on a closed/periodic domain. (Probe: torus `0.20043 → 0.20043`.)

What *does* depend on shape is the **accuracy** of the transported field: on the
flat cube and topological torus it is FE-accurate; on the geometric sphere it
degrades gracefully (the dynamics lose accuracy with tet quality, but the solve
stays bounded and the constraints stay exact) — exactly the `05` §5 lesson, now for
a fluid.

---

## 2. The natural flows, defined

- **Cube / torus — uniform advection.** `u = (U,0,0)`, constant; on the torus the
  x-faces are identified so the flow circulates. Divergence-free trivially.

- **Sphere — solid-body rotation.** `u(x) = ω × (x − x_c)` about the centre. This is
  exactly divergence-free (`∇·(ω×r) = 0`) and tangent to spherical shells, so a
  scalar blob is carried around without leaving the ball — the natural closed flow
  on a sphere, needing no wall.

Each is a *prescribed* velocity (Stage B style): the scalar is passive, so the
problem is linear in `φ` and verifiable, while still exercising advection +
projection + diffusion on the shaped mesh.

---

## 3. Verification design (`fluid_shapes_verify.py`)

All four PASS — cube `||Du||=2.2e-14` / mass `2.2e-16`; torus skew `1.3e-17` /
one-lap mass `2.8e-16`; sphere (2844 sign-flips) `||Du||=1.2e-14` / rotation energy
ratio `1.000`; cross-shape incompressibility machine on all. The checks:

- **A. Cube (reference).** Uniform advection–diffusion: mass conserved (machine),
  projection `‖Du‖` machine, Gaussian variance `2κt` (FE). Confirms the baseline.

- **B. Torus (ring).** Uniform advection around the identified ring: **one-lap mass
  conservation** to machine precision, skew energy `φᵀC_skew φ = 0`, variance growth
  FE-accurate. The flat (topological) ring matches the cube's accuracy.

- **C. Sphere (ball, degraded).** Solid-body rotation advection on the geometric
  ball: the **projection keeps `∇·u=0` to machine precision despite ~2844
  sign-flipped cotangents** (the headline); the advection stays **bounded/stable**
  (graceful degradation, energy ratio bounded); mass behaviour reported.

- **D. Cross-shape invariant.** State the unifying result: incompressibility holds
  to machine precision on *all three* shapes, and mass is conserved on the closed
  ones — the fluid framework transforms, with accuracy (not constraints) tracking
  mesh quality.

---

## 4. What is and isn't claimed

**Design claims (to verify in code):**
- The V3 fluid (advection + diffusion + projection) runs on cube, torus, and sphere.
- Incompressibility and mass conservation are shape-independent (machine precision),
  including on the degraded sphere; transported-field accuracy tracks tet quality.

**Not claimed:**
- `u` is **prescribed** (passive scalar) here, as in Stage B; coupled nonlinear flow
  on shapes is later (after Stage C on the cube).
- The sphere/torus use the **geometric** (degraded) meshes; the topological torus is
  exact, the geometric ball degrades gracefully. No claim that the geometric ball is
  high-accuracy — only that it solves stably with exact constraints.
- No turbulence, no walls-with-no-slip on curved shapes (the natural closed flows are
  used instead); those are out of scope.

---

## 5. Files and next steps

- Theory: this document, on `05`/`06a`/`06b`/`04`.
- **Verification:** `src/verification3d/fluid_shapes_verify.py` (checks A–D).
- Then: extend the fluid viewer to cube/torus/sphere (now justified by this
  verification), reusing the solid-viewer rendering; and Stage C (full NS) on the
  cube.
