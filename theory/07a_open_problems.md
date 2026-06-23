# 07a · Open Problems — Remaining Tasks for V3 and Beyond

Enumerated and prioritised. Each item states *what*, *why it matters*, and *how it could
fail* — because, per `06z`, the tasks that can fail against an external standard are the
valuable ones. Priority: **P0** = load-bearing (do before extending), **P1** = important,
**P2** = enrichment / nice-to-have.

---

## P0 — load-bearing  [STEP 1/2 DONE: cavity matches Ghia Re=100, see 08]

**1. Trustworthy pressure (inf-sup-stable / stabilised).**
*What:* replace equal-order P1/P1 with a Taylor–Hood (P2–P1) or MINI pair, or adopt PSPG
properly (not just as a diagnostic). *Why:* the pressure is currently **not trustworthy**
(`06c`, `06z`); without it, "incompressible NS" holds for velocity only. *Failure mode:*
the stable pressure may expose velocity errors currently masked by the projection, or
break the machine-precision incompressibility — either would be informative.

**2. Lid-driven cavity benchmark.**
*What:* the standard wall-bounded test (Re = 100, 400, 1000) against published Ghia et al.
data. *Why:* converts "periodic, analytic, velocity-only" into a **real, externally
checkable** flow; this is the single most decisive validation. *Failure mode:* corner
singularities, wall boundary layers, or the unstable pressure could make the centreline
profiles miss the published data — exactly the reality-contact the project needs.

---

## P1 — important (breadth and honesty)

**3. Flow past a cylinder (Strouhal number).**
*What:* unsteady wake, vortex shedding, measure St(Re). *Why:* tests genuinely unsteady,
non-symmetric dynamics with an experimental number to hit. *Failure mode:* the moderate-Re
laminar scheme may not capture shedding frequency without proper outflow BCs.

**4. Nonlinear Rayleigh–Bénard convection.**
*What:* full (not linearised) thermal convection; measure the critical Rayleigh number
`Ra_c ≈ 1708`. *Why:* `06g`/`06h` are linear/steady; this is the nonlinear feedback where
flow advects temperature that controls flow. *Failure mode:* the symplectic/projection
machinery may not preserve the convective instability threshold accurately.

**5. Complex / non-periodic boundary conditions.**
*What:* inflow/outflow, slip, Robin. *Why:* current tests are periodic + simple Dirichlet;
real problems need more. *Failure mode:* mass conservation (which `06b` showed depends on
periodicity) may degrade at open boundaries.

**6. Long-time enstrophy / energy-cascade behaviour.**
*What:* run Taylor–Green to long times, track enstrophy `∫|ω|²`. *Why:* the agent flagged
this; short runs can hide slow drift. *Failure mode:* secular enstrophy growth would
reveal accumulating high-frequency contamination not seen in the short PSPG test.

---

## P2 — enrichment and structure

**7. 3-D high-order (`Pn`) stiffness assembly.**
*What:* implement the enrichment of `06f` in 3-D (it is currently demonstrated in 1-D).
*Why:* the verified route to machine precision for smooth modes. *Failure mode:* 3-D `Pn`
assembly cost / conditioning may not deliver the clean 1-D spectral rate.

**8. Mesh-quality cure for sign-flips / obtuse dihedrals.**
*What:* the negative-cotangent problem on deformed/sphere meshes (`05`, `cotangent_signflip`).
*Why:* the sphere is still "degraded (2844 sign-flips)"; a well-centred / intrinsic
triangulation would fix it. *Failure mode:* may require changing the meshing approach,
not a local patch.

**9. First-principles `ν` (molecular Coulomb → Green–Kubo).**
*What:* the bridge mapped in `06e` but not built — a small MD + stress-autocorrelation
study. *Why:* would close the "honey vs water" chain from first principles. *Failure mode:*
a different engine entirely (statistical mechanics); may stay out of this notebook's scope.

**10. Thermo-viscous nonlinear feedback.**
*What:* time-dependent flow advecting a temperature field that controls `ν(T)` which
controls the flow (`06h §5`). *Why:* the genuine fusion the glue enables. *Failure mode:*
coupled nonlinear stiffness may lose the machine-precision interface results.

---

## Deferred by intent (not tasks yet)

- **The coupling-setting-field resemblance** (Higgs-like) noted in `06h` is **not** to be
  formalised until there is something to verify — recorded only as a direction.

---

## Cross-cutting (engineering, not physics)

- **Viewer parity:** keep `viewer/viewer_model.html` and the per-stage viewers showing
  *verified* numbers only (presets from solver output, no in-browser physics).
- **Reproducibility:** the 20 verification scripts should remain a single re-runnable
  suite (the `06z` audit re-ran them all; keep that easy).
- **Synthesis-drift guard:** any new `06x`/`07x` note must carry its scope caveats inline,
  so the narrative cannot outrun the proof (the lesson of `06z`).

---

## The one-line recommendation (from `06z`)

**Consolidate (this synthesis), then GO narrowly: P0 first — a trustworthy pressure and
the lid-driven cavity.** If the cavity matches published data, the foundation is earned
and broad GO is justified; if not, we learn where the edifice meets reality before
building on top.
