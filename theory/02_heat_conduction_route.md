# 02 · The Heat-Conduction Route (the integrating foundation)

**Status:** exploration findings + decision (backed by code). This document
records where the scatter form and the FE/DEC operator **integrate**, prompted by
two leads: *(a) look hard for the place the scatter form is consistent with
FE/DEC; (b) treat it as heat — conduction can be computed from contact area.* The
exploration was done with runnable code
([`heat_conduction_verify.py`](../src/verification3d/heat_conduction_verify.py));
every claim below is a measured number, not an assertion.

**One-line result.** The place they integrate is **heat (diffusion)**.
"Conduction from contact area" is correct — *as the cotangent/FE conductance* —
and because diffusion needs only the metric (areas/lengths), not face directions,
the heat route **sidesteps the spurious wave anisotropy** that 01e showed could
not be patched on the scatter node.

---

## 1. Why heat, not waves, is where it integrates

01e established, with numbers, that the scalar **wave** scatter node carries a
real ~50% directional anisotropy that area-weighting does not remove; only the
full FE/DEC operator is ~10× more isotropic. The obstruction was that a wave
carries *directional* information through faces, and the geometry-blind /
area-only node loses the face normals.

**Diffusion is different.** Fourier's law `q = −k∇T` produces a flux whose
inter-cell exchange is governed by a scalar **conductance**: how much heat crosses
a contact, not which way it points. Heat between two regions is

```
Q = c · (T_j − T_i) ,     c = conductance of the shared connection.       (1)
```

The user's intuition — "as heat, conduction can be computed from contact area" —
is exactly (1). The only question is what `c` is, geometrically.

---

## 2. Contact area, done two ways (one wrong, one right)

### 2.1 The naive way (inconsistent on tetrahedra)
Cell-centred two-point flux: `c = A_face / d(centre_i, centre_j)` — contact area
over centroid distance. Measured cube Neumann spectrum
([`heat_conduction_verify.py`](../src/verification3d/heat_conduction_verify.py) A):

```
exact lowest nonzero eigenvalue = π² = 9.8696
naive area/distance FV →  ~7.59   (wrong by ~23%, and anisotropic)
```

It *converges* — but to the **wrong operator**. Reason: the two-point flux
approximation is only consistent on **orthogonal** meshes (centre-to-centre line
⟂ shared face). Tetrahedra are non-orthogonal, so centroid distance is the wrong
length and the naive conductance is inconsistent. (This is the finite-volume
"TPFA fails on non-orthogonal grids" fact.)

### 2.2 The right way (the cotangent / FE conductance)
Use the **correct geometric conductance** — contact area over the *dual*
(circumcentric) length, i.e. the 3-D cotangent weight. That conductance is exactly
the off-diagonal of the P1 FE stiffness matrix, `c_ij = −K_ij`. Measured:

```
cotangent / FE heat operator →  π² (correct),  Δc/c ≈ 0.048 anisotropy
                                 (~10× more isotropic than the wave node)
```

So **"contact-area conduction" is right** — provided the area is divided by the
*right* (dual) length. With that, the heat operator **is** the FE/DEC Laplacian.
This is the integration point: the scatter-form material idea and FE/DEC meet, as
a **thermal conductance network**.

---

## 3. Materials are exact in the heat route

The whole motivation for the scatter form was material contrast (metal/wood/water
reflecting energy differently). In the heat route this becomes **conductivity**
`k`, and the two-material steady state is exact
([`heat_conduction_verify.py`](../src/verification3d/heat_conduction_verify.py) B):

```
interface temperature  T_iface = k₂/(k₁+k₂)   (series thermal resistance)
measured error  ≤ 1e-14   for k-ratios 1, 3, 10, 0.1
```

Unlike the wave reflectance (which needed a geometry-robust total-energy trick and
was anisotropy-limited), the heat material law is recovered to **machine
precision**. Materials "just work" in the heat route.

---

## 4. The conductance network is physical (mostly)

For (1) to be a real heat network it needs **non-negative** conductances `c_ij ≥ 0`
(discrete maximum principle: heat flows hot→cold, never the reverse). Measured on
the Kuhn mesh
([`heat_conduction_verify.py`](../src/verification3d/heat_conduction_verify.py) C):

```
negative-conductance fraction ≈ 0% on the Kuhn cube
```

So the Kuhn cotangent network is physical. Where a mesh *does* produce negative
cotangent conductances (poorly-shaped / non-Delaunay / non-well-centred tets), the
network violates the maximum principle — and that is, once again, the
**mesh-quality / "equally-spaced tetrahedra"** thread that recurs throughout this
project. Heat conduction gives that thread a crisp physical meaning: *a good mesh
is one whose dual conductances are all non-negative.*

---

## 5. What this resolves and what it opens

**Resolved (with code):**
- The scatter form integrates with FE/DEC **as heat**: contact-area conduction =
  cotangent/FE conductance.
- The 01e wave-anisotropy obstacle is **avoided** in the heat route (diffusion
  needs only the metric).
- Material contrast is **exact** as conductivity (T_iface = k₂/(k₁+k₂)).
- Mesh quality acquires a physical criterion: **non-negative dual conductances.**

**Honest scope / not claimed:**
- This is still classical (FE/DEC heat equation, finite-volume on the dual).
  Nothing novel; it is the *correct* and *consistent* foundation, which the naive
  area/distance form was not.
- The heat route is **scalar diffusion**, not yet fluid momentum. But it is the
  right rung: Step 2 (pressure/tank) is a diffusion-like elliptic solve, and the
  heat operator is precisely the consistent Laplacian that pressure projection
  needs. So this also de-risks Step 2.

---

## 6. Proposed next step (for approval, theory-first as always)

Adopt the **thermal conductance network** (cotangent/FE) as the working operator
and build Step 2 (pressure/tank) on it, since the incompressibility projection is
a Laplace solve that this operator gets right. The material/`S`-parameter idea
lives on as **per-edge conductivity**. Concretely, the next theory doc would be
`theory/03_pressure_on_heat_operator.md`: pressure as the scalar field whose
Laplacian (the verified heat operator) enforces `∇·u = 0`, with a tank (closed
domain) as the first exact test (eigenmodes / equilibration).

---

### Roadmap position

```
[Step 1]   S-parameter wave        01, 01b (verified); 01c/01d reflectance closed
[Step 1e]  geometric consistency    01e (wave anisotropy: needs full FE operator)
[Step 2pre] HEAT route              02  ← THIS DOCUMENT (integration point found)
            heat_conduction_verify  ✓ (consistency, exact materials, valid network)
[Step 2]   pressure / tank          (build on the verified heat/cotangent operator)
[Step 3]   fluid dynamics
```

_The heat route is the consistent foundation the project was circling toward:
the scatter-form material idea, made geometrically correct, is a thermal
conductance network = FE/DEC. Next, on approval, pressure is built on it._
