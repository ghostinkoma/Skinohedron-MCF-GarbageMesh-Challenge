# 09 · Capstone — From a Garbage Mesh to a Classical Benchmark

A short integrating note across the whole project: where it started, what one operator
turned out to carry, and the moment the work met an external standard. It does not
re-prove anything (each claim points to its document and script); it states the shape of
what was built and, honestly, its bounds.

---

## 1. The arc

The project began with a blunt question — can a *garbage* tetrahedral mesh still do
physics? — and answered it by refusing to stack theory on unverified theory. The path:

- **Foundation (`01`–`05`).** Reject the geometry-blind scatter node (`01e`); do
  conduction correctly on the cotangent/FE operator (`02`, `π²`, exact two-material
  interface); bet that **one** operator `L = M⁻¹K` carries heat *and* waves (`03`); add
  pressure/incompressibility (`04`) and shape transforms (`05`).
- **V3 fluid (`06`–`06h`).** Stokes → advection → full Navier–Stokes (`06a`–`06c`), on
  cube/torus/sphere (`06s`); then a domain map — inertial mass exact, Coulomb/stiffness
  `O(h²)` (`06d`–`06f`), temperature both a fused energy domain and a multiplicative glue
  (`06g`–`06h`).
- **Audit (`06z`, `07`, `07a`).** Treat it all as hypothesis: re-run every script (all
  reproduce), fix one framing imprecision, concede the phantom-agent critique in full
  (pressure unreliable, `‖Bu‖` trivial, benchmarks untested), and recommend: consolidate,
  then validate against an external standard before extending.
- **P0 (`08`–`08b`).** Do exactly that. A **trustworthy pressure** (PSPG, then
  cross-checked by parameter-free Taylor–Hood) on the **lid-driven cavity**, matching the
  **Ghia (1982)** benchmark at `Re=100, 400, 1000`, with the secondary-vortex "rev"
  reproduced and the primary-vortex centre landing on Ghia's published coordinates.

---

## 2. What one operator turned out to carry

`L = M⁻¹K` is read, with the *same* `K`, as heat, scalar waves, pressure/incompressibility,
electrostatics (Coulomb), viscous momentum, temperature diffusion, and — stabilised — the
cavity pressure. The error structure is uniform and now well-understood: **linear fields
exact, polynomial-source nodally exact, sinusoidal eigenmodes `O(h²)` → machine precision
by local enrichment**. Domains **fuse** (each adjoint-consistent coupling adds an exact
invariant) and one of them, **temperature**, also **mediates** the others. This is the
"one tetrahedron carries many physics" idea, made concrete and checked.

---

## 3. The moment that matters

For most of its life the project validated against *its own* exact solutions — rigorous,
but internal. `08`–`08b` is the first time the output was held against numbers the project
did **not** generate: the Ghia cavity tables. It matched them — at three Reynolds numbers,
with two independent elements that agree with each other more tightly than with the table.
That an amateur-driven, self-built operator reaches the classical benchmark numbers is a
real result; it is the difference between an elegant internal story and a flow that the
wider world would recognise as correct.

It is worth saying plainly, without overstating it: this does not make the notebook a
competitive solver (see `07`/`06z` — no turbulence, limited geometry, small scale). It
makes the *foundation* earned.

---

## 4. Honest bounds (unchanged from the audit)

Steady laminar only; `Re≤1000`; 2-D for the benchmark (3-D demonstrated structurally);
pressure trustworthy *when stabilised*; the unifying narrative remains interpretation on
idealised tests, broader than any single proof. No new mathematics anywhere — the
contribution is verification discipline and conceptual unity, not novelty.

---

## 5. The open frontier — why does the "rev" exist?

`08b` found the **rev**: the flow reorganises with Re (corner vortices grow, the primary
centre migrates) and the two discretisations begin to part company as Re rises. The next
research thread asks *why* such structural change exists, and what — if anything — is
**missing** from the model that a fuller physics would supply. Candidate sources, to be
separated by experiment (see the rev-investigation plan):

1. **Approximation error** — the rev as an artifact of finite `h` / `O(h²)` (the
   discretisation straining), removable by enrichment.
2. **A missing domain** — e.g. Brownian/electrostatic effects not represented.
3. **A field-intrinsic structure** — the rev as a property of the continuous field, not
   the mesh.
4. **Neglected chemistry** — reaction terms absent from a purely mechanical model.

The disciplined next step is to design tests that tell these apart — the same
theory→model→code→viewer loop, now aimed at the origin of the rev.

---

## 6. Index

Foundation `01`–`05` · V3 fluid `06`–`06h` · audit `06z`/`07`/`07a` · P0 cavity
`08`/`08a`/`08b` · this capstone `09`. Viewers unified responsive under `viewer/`;
structure documented in `STRUCTURE.md`.
