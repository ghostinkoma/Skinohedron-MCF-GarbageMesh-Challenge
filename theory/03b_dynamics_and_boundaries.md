# 03b · Dynamics, Boundaries, and the Parametrised Instrument (Step 1.5)

**Status:** consolidation document (theory), backed by code and a runnable
viewer. It continues `03` (the single operator `L = M⁻¹K`) by adding the two
things that turn `L` from a static spectrum into a *physics you can drive and
watch*: (1) the **time structures** on `L` with their parameters (damping,
forcing, dispersion limit), and (2) the **boundary conditions** that couple the
cube to an outside world (convection and radiation). Everything here is grounded
in numbers measured from the same `L` used in `viewer/viewer_unified.html`. This
is the Step 1.5 foundation for V2.

Mesh used throughout: the equally-spaced Kuhn cube, `n = 8`
(729 vertices, 3072 tetrahedra), unit cube centred at the origin.

Key measured constants for this mesh (from `data_unified.js`):

```
maxKM   = max_i (Kdiag_i / M_i)  ≈ 512
lamMax  = 2.2 · maxKM            ≈ 1126.5      (estimate of largest eigenvalue of L)
dtWave  = 0.4 · sqrt(4 / lamMax) ≈ 2.384e-2    (leapfrog step, wave)
dtHeat  = 0.35 / maxKM           ≈ 6.835e-4    (explicit Euler step, heat)
```

---

## 1. The two time structures (recap) and their stability

On the one operator `L`:

```
heat:        M ṗ  = −K p              explicit Euler:   p ← p − dtHeat · M⁻¹ K p
scalar wave: M p̈  = −K p              leapfrog:         p⁺ = 2p − p⁻ − dtWave² · M⁻¹ K p
```

Both are conditionally stable; the steps above are chosen below the CFL limit of
this mesh. The wave leapfrog conserves a discrete energy to ~1% over thousands of
steps (verified in `03`); heat decays monotonically as expected.

---

## 2. Damped wave: `M p̈ = −K p − γ M ṗ`

Adding viscous damping `−γ M ṗ` gives the leapfrog update (let `a = γ·dtWave/2`):

```
p⁺ = ( 2p − (1−a) p⁻ − dtWave² · M⁻¹ K p ) / (1 + a)
```

**Behaviour.** Each eigenmode `λ` of `L` becomes a damped oscillator: it
oscillates at `cos(√λ · t)` under an envelope `exp(−γ t / 2)`.

- `γ = 0` → no decay; waves persist indefinitely. This is the "vacuum / nothing
  to dissipate into" picture — the closest thing here to a wave spreading freely
  through space.
- `γ > 0` → amplitude follows `exp(−γ t / 2)`; half-life `t½ = ln 2 / (γ/2)`.

**Verification** (lowest non-zero eigenmode seeded, 400 steps, tail/initial
amplitude ratio vs. predicted envelope `exp(−γ·dt·N/2)`):

| γ   | measured ratio | predicted envelope |
|-----|----------------|--------------------|
| 0.0 | 1.011          | 1.000  (persists)  |
| 0.5 | 0.047          | 0.044              |
| 2.0 | ~0.000         | ~0.000             |

The measured decay tracks the analytic envelope. `γ = 0` is genuinely
non-dissipative within discrete-energy fluctuation. **PASS.**

---

## 3. Continuous source and the dispersion limit

A continuous (held) source drives one vertex with `s(t) = A · sin(2π f · step)`,
where `f` is in **cycles per step**. A single click instead deposits one Gaussian
impulse (the `pulse` mode).

**Dispersion / propagation cutoff.** A discrete wave operator supports real
travelling waves only up to a maximum angular frequency `ω_max = √λ_max`. Above
it there is no real wavevector and a forced source is **evanescent** — energy
stays pinned near the source instead of radiating. In cycles per step:

```
f_cutoff = ω_max · dtWave / (2π) = sqrt(lamMax) · dtWave / (2π) ≈ 0.127 cyc/step
```

The viewer's frequency slider spans `0.003 … 0.120 cyc/step`, i.e. deliberately
**below** this cutoff, so the whole slider range propagates. The readout reports
an approximate wavelength `λ ≈ c / f` (with phase speed `c ≈ 2·f_cutoff ≈ 0.254`
cells/step) and flips to "evanescent" if `f` exceeds the cutoff. Lower `f` → long
waves (wavelength can exceed the box, ~12–13 cells at `f = 0.02`); higher `f` →
shorter waves and visible interference, until the cutoff where propagation stops.

This is the honest statement: **the mesh has a built-in shortest wavelength
(~2 cells, Nyquist), and the cutoff frequency is its image in time.**

---

## 4. Boundary conditions: coupling the cube to an outside world

Until now the cube was closed (natural/insulating boundary: no flux leaves). To
let heat exchange with an environment at temperature `T_env`, we add a boundary
term carried only by surface nodes.

**Boundary mass.** Let `B_i = M_i` if vertex `i` is on the cube surface, else `0`
(surface nodes: 386 of 729, i.e. 9³−7³). `B` weights the exchange by how much surface each
node represents. Because `M⁻¹ B = 1` on the boundary and `0` inside, the boundary
term acts cleanly only where the cube meets the outside.

### 4a. Convection — Robin boundary: `M ṗ = −K p − h·B·(T − T_env)`

Surface nodes relax toward `T_env` at a rate set by the heat-transfer
coefficient `h`. This is a Robin (mixed) boundary condition.

**Behaviour.** The whole body exponentially approaches `T_env`; the relaxation
half-life is inversely proportional to `h`.

**Verification** (body initialised uniform `T=1`, `T_env=0`, half-life =
steps to mean `T = 0.5`):

| h    | half-life (steps) | stable |
|------|-------------------|--------|
| 0.2  | 15371             | yes    |
| 0.5  | 6154              | yes    |
| 1.0  | 3082              | yes    |
| 2.0  | 1546              | yes    |
| 4.0  | 778               | yes    |
| 8.0  | 393               | yes    |
| 16.0 | 200               | yes    |
| 32.0 | 102               | yes    |

Half-life `× h` is constant (≈ 3082 step·units): a clean `t½ ∝ 1/h` law.
Stability limit is `h ≲ 2/dtHeat ≈ 2930`, far above the usable range. **PASS.**

Viewer presets: **water** `h = 4.0` (strong convection), **air** `h = 0.8`
(weak). The readout shows `half-life ≈ 3082/h steps`.

### 4b. Radiation (vacuum) — Stefan–Boltzmann: `M ṗ = −K p − εσ·B·(T⁴ − T_env⁴)`

In vacuum there is no convection; the only way to lose heat is to radiate it. The
flux scales as `T⁴`, so the cooling term is non-linear.

**Behaviour.** Cooling is **fast when hot, slow when cool** — qualitatively
different from convection's single time-constant.

**Verification** (εσ = 1, `T_env = 0`, mean temperature every 4000 steps):

| start T₀ | t=0    | 4000   | 8000   | 12000  | 16000  |
|----------|--------|--------|--------|--------|--------|
| 2.0 (hot)| 1.9942 | 0.7133 | 0.5685 | 0.4971 | 0.4519 |
| 1.0      | 0.9996 | 0.6478 | 0.5394 | 0.4795 | 0.4396 |
| 0.5 (cool)| 0.5000| 0.4537 | 0.4209 | 0.3959 | 0.3760 |

The hot body sheds heat dramatically faster (Δ ≈ 1.28 in the first window) than
the cool body (Δ ≈ 0.05): the `T⁴` signature. **PASS.**

Viewer preset: **vacuum** `εσ = 1.5`, boundary mode = radiation. The readout
reads "T⁴ radiative: fast when hot, slow when cool".

---

## 5. The parametrised instrument (`viewer/viewer_unified.html`)

The viewer makes the above touchable. It renders the field on the tetrahedral
mesh as a vertex-interpolated **heatmap** (colormap blue→green→yellow→orange→
red→purple→white; intensity auto-normalised to the current peak so the field is
always visible), with:

- **mode**: wave / heat — the same `L`, two time structures.
- **wave controls**: source = continuous / pulse; **frequency** slider (with
  wavelength + propagating/evanescent readout); **damping γ** slider (with
  half-life readout; γ=0 = persistent).
- **heat controls**: outside **environment** = water / air / vacuum (selects
  convection vs. radiation and a default coefficient); **environment temperature**
  `T_env`; **exchange rate** (`h` or `εσ`, 0–20) with half-life / T⁴ readout.
- **interaction**: hold the mouse on the cube to excite/heat while held; release
  to let it evolve and cool (no reset — new excitation interacts with the
  existing field). Drag on the background to rotate.
- **view**: auto-rotate; **wireframe** (the complete geometric edge set — see
  below); **slice** (clip plane to see the interior) with position; **transparency**.

### Two viewer facts worth recording (hard-won)

- **Wireframe edges ≠ conductance edges.** The conductance graph (`K` off-
  diagonals) drops edges whose cotangent weight is ~0 at right-angle dihedrals —
  for `n=8` that is **2240 of 4184** edges missing. The wireframe must use the
  full geometric edge set `wedges` (4184 unique edges, all 729 vertices present,
  verified complete). Drawing the wireframe from conductance edges silently
  yields a broken-looking mesh that is not actually broken data.
- **Wireframe as indexed lines in the same program.** The robust pattern (matching
  the repo's working `viewer.html`) is to draw the wireframe with `drawElements
  (gl.LINES)` into the *same* vertex buffers as the faces, lifted slightly toward
  the camera (`uLift`) and flat-coloured (`uFlat`), rather than a second program
  with its own buffers. Separate-program + `drawArrays` is fragile: leftover
  enabled vertex-attribute arrays sized for the faces fail WebGL's array-size
  validation on the larger line draw, discarding the whole call.

---

## 6. What is and isn't claimed (Step 1.5 honesty)

**Established and verified:**

- One operator `L` drives heat (decay) and waves (oscillation) on the same
  spectrum (`03`).
- Damping `γ` produces the analytic envelope `exp(−γt/2)`; `γ=0` is
  non-dissipative.
- The discrete dispersion cutoff `f_cutoff ≈ 0.127 cyc/step` is the temporal
  image of the mesh's ~2-cell Nyquist wavelength.
- Robin convection gives `t½ ∝ 1/h`; Stefan–Boltzmann radiation gives `T⁴`
  cooling. Both verified numerically and stable in the usable range.

**Not claimed / honest limits:**

- No new mathematics: damped waves, Robin and Stefan–Boltzmann boundaries are
  classical. The contribution is a *verified, self-consistent instrument* on this
  project's own `L`.
- Residual spatial anisotropy (~5%, Δc/c ≈ 0.048 from `01e`/`03`) is inherited
  from the Kuhn mesh's broken symmetry; it is unchanged here.
- All coefficients (`h`, `εσ`, `γ`, `f`) are in **simulation units** tied to the
  mesh step sizes above; they are physically *meaningful in ratio and scaling*
  (the verified laws), not calibrated to SI.

---

## 7. Files

- Theory: this document, on top of `01`–`01e`, `02`, `03`.
- Operator / data: `src/verification3d/export_unified_viewer.py` →
  `viewer/data_unified.js` (emits `verts`, conductance `edges`/`cond`, lumped
  `mass`, `kdiag`, boundary triangles `btris`, all triangles `atris`, and the
  full wireframe edge set `wedges`).
- Instrument: `viewer/viewer_unified.html`.
- Verification of every number above: re-runnable against `data_unified.js`
  (damped-wave envelope, dispersion cutoff, Robin `1/h` half-life and stability,
  Stefan–Boltzmann `T⁴` cooling).

**Next (Step 2):** pressure on the same `L` — pressure as the scalar field whose
Laplacian (this verified `L`) enforces `∇·u = 0`, with a closed tank as the first
exact test (`04_pressure_tank.md`).
