# Repository Structure & Conventions

This file codifies the existing organisation so the system is explicit and uniform.
It documents conventions rather than moving files (the project archives and appends
rather than deletes).

## Top-level layout

```
theory/      theory documents (numbered, see convention below) + README index
src/         source: core operators, verification scripts, exporters
viewer/      self-contained HTML viewers + their data_*.js (verified presets)
results/     saved numerical outputs
paper/       write-up material
docs/        supplementary documentation
tests/       test scaffolding
archive/     superseded material (kept, not deleted)
review/      review notes
examples/    examples
STRUCTURE.md this file
README.md    project entry point
```

## src/ layout

```
src/ksf/            V1/V2 core (2-D / early operators)
src/ksf3d/          V3 core: mesh3d_uniform.py (kuhn_cube), fem3d.py (fem_laplacian -> K, lumped M), dec3d.py
src/verification/   V1/V2 verification scripts
src/verification3d/ V3 verification + viewer-export scripts
src/viewer/         (reserved)
```

## Naming conventions

**Theory documents** — `NN[a-z]_short_name.md`, where `NN` is the stage:
- `00` synthesis/prologue · `01`–`05` foundation · `06`–`06h` V3 fluid · `06z` audit ·
  `07`/`07a` synthesis + open problems · `08`–`08b` P0 cavity · `09` capstone.
- Letter suffixes (`a`,`b`,…) are sub-stages of the same number; `z` is reserved for an
  audit of that stage. Bilingual docs use `.en.md` / `.ja.md`.
- Every theory doc is listed in `theory/README.md` with a one-line status and its script.

**Verification scripts** — `src/verification3d/<topic>_verify.py`. Each is runnable
stand-alone, prints `PASS`/`FAIL` per check, and asserts against an exact solution, a
conserved invariant, or an external benchmark. Demonstrations (not benchmarks) end in
`_demo.py`.

**Viewer-export scripts** — `src/verification3d/export_<topic>_viewer.py`: produce a
`viewer/data_<topic>.js` of *verified* solver output.

**Viewers** — `viewer/viewer_<topic>.html` + `viewer/data_<topic>.js`. No physics runs in
the browser; viewers render verified presets only. All viewers are responsive (a unified
`@media (max-width:760px)` block; container pattern `#wrap` flex + `#side`/`#sidebar` +
`#stage`).

## The one operator

All physics is read from a single discrete operator `L = M⁻¹K` on the Kuhn-cube P1 mesh:
`fem_laplacian(V,T)` in `src/ksf3d/fem3d.py` returns `(K, M)` with **M a lumped-mass
vector** (use `M*x` elementwise, `np.sum(M*x)` for integrals, `1/M` for the inverse).
`kuhn_cube(n)` in `src/ksf3d/mesh3d_uniform.py` builds the mesh.

## Conventions for new work

1. Follow theory→model→code→viewer; no claim without a passing script.
2. New theory doc → add a `theory/README.md` row; new viewer → `data_*.js` from a verified
   exporter; keep the responsive block.
3. Negative results are preserved, not deleted.
4. Each stage may carry its own `..z`-style audit; caveats live inline so the narrative
   cannot outrun the proof.
