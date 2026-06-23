# 12 · Residual Domains — Coulomb (2) and Reaction (4), and What NS Still Lacks

**Status:** **verified**, `src/verification3d/residual_domains_verify.py`. Implements the
method `10` pointed to: a candidate missing physics is confirmed when the purely
mechanical model leaves a **residual** that the candidate term removes. Builds the two
domains the rev investigation flagged — **Coulomb (2)** and **reaction (4)** — shows their
fingerprints, shows they compose with enrichment **(1)**, and uses the exercise to **enumerate
what the present Navier–Stokes model still lacks**.

---

## 1. Hypothesis (2) — Coulomb / electrohydrodynamics

A charged fluid: charge density `ρ_q` sets an electrostatic potential through
`K φ = M ρ_q` — the **same `K`** verified for electrostatics in `06e` — and the Coulomb
body force `f = −ρ_q ∇φ` enters the momentum equation. Measured (n=48, Re=200, dipolar
charge, *no wall driving*):

- **Coulomb OFF (pure mechanical NS): max speed `0`** — with no walls moving, the
  mechanical model produces no flow at all.
- **Coulomb ON (EHD): max speed `0.007`**, scaling **linearly** with charge
  (`×0.5→0.003, ×1→0.007, ×2→0.013`).

The mechanical model's residual (zero flow) is filled entirely by the Coulomb domain.
This is electrohydrodynamics, built from the operator the project already has.

---

## 2. Hypothesis (4) — reaction

A species `c` advected–diffused by the cavity flow, with a Fisher–KPP reaction
`R = k·c(1−c)`. Measured (n=48):

- **No reaction (transport only):** total `c` conserved (`0.060→0.060`), it merely
  dilutes; `c>0.5` region = 92 nodes.
- **With reaction (autocatalytic):** total `c` **grows** (`0.060→0.125`), saturates to
  `c→1`, and a **front propagates** (`c>0.5` region = 307 nodes).
- residual `‖c_react − c_noreact‖ = 0.142`.

Pure transport cannot create the front; the reaction term fills that residual.

---

## 3. (1)+(2) compose — enrichment lifts the Coulomb-driven flow

Because `φ` rides the **same `K`**, refinement and high-order enrichment improve the EHD
flow just as they improve the mechanical flow. Measured: the Coulomb-driven max speed
converges under refinement (`0.00660 → 0.00657 → 0.00658`, Cauchy differences
`0.00003 → 0.00001`), and `φ`'s `O(h²)` error → machine precision by the `06f` enrichment.
So **(1) refinement/enrichment and (2) the Coulomb domain compose** — the accuracy route
and the missing-domain route stack.

This is exactly the comparison the existing precision viewer is built for: the same
"refine/enrich → error down" curves now apply to the Coulomb-driven flow, not only the
mechanical one.

---

## 4. What the present NS model still lacks — enumerated

The residual method, applied here, makes the **remaining gaps in the NS modelling**
concrete. Each is a residual the current mechanical, single-phase, isothermal,
non-reacting, non-charged incompressible NS cannot reproduce, with the term that fills it:

1. **Body-force / multiphysics coupling** — Coulomb (`§1`, EHD), buoyancy (`06g`),
   magnetic (MHD). *Fingerprint:* a flow with no mechanical drive. *Status:* Coulomb &
   buoyancy demonstrated; MHD not.
2. **Reactions / source terms** — `§2`. *Fingerprint:* species growth/fronts. *Status:*
   demonstrated (passive scalar); reaction→flow feedback (heat release → buoyancy) not.
3. **Free surface / multiphase** — pouring, splashing, bubbles. *Fingerprint:* a moving
   interface. *Status:* **not modelled** (the pour viewer `11` is illustrative only). This
   is the single biggest gap for the "pouring water" goal.
4. **A trustworthy pressure beyond Re=1000 / unsteadiness** — vortex shedding, transition
   (`~Re 8000`). *Fingerprint:* time-dependent wake. *Status:* steady laminar only;
   time integration of the full NS not yet built.
5. **Non-Newtonian / variable viscosity feedback** — `ν(T)` advected by the flow it
   controls (`06h` is steady). *Fingerprint:* shear-thinning structure. *Status:*
   interface verified, the dynamic feedback not.
6. **Compressibility / acoustics** — the wave route (`03`) exists but is not coupled into
   the incompressible solver. *Fingerprint:* density waves. *Status:* separate, not fused.
7. **3-D at scale** — the benchmark is 2-D; 3-D is structural-demo only (`08`). *Status:*
   real 3-D high-`Re` needs the HPC machinery the project does not have.

The first three are the live frontier for the project's own goals (pouring, charged/
reacting flows); 4–7 are the standard CFD distance noted in `06z`/`07`.

---

## 5. What is and isn't claimed

**Verified:** the Coulomb (EHD) and reaction fingerprints exist and are removed by their
terms; the Coulomb domain reuses `K` (`06e`); (1)+(2) compose (EHD flow converges under
refinement).

**Not claimed:** these are **constructed** settings showing the *method*, not evidence
that the cavity rev needs them (`10` showed it does not). No claim that the demonstrated
EHD/reaction match any specific experiment — the magnitudes are illustrative; the
*structure* (residual filled by the term) is the verified content. Free surface, MHD,
reaction→flow feedback, and unsteady high-Re remain unbuilt (`§4`).

---

## 6. Files

- Theory: this note, on `06e`/`06f`/`10`.
- Verification: `src/verification3d/residual_domains_verify.py` (checks (2),(4),(1)+(2)).
- Viewer: the precision viewer (`11`) framework compares the refine/enrich curves,
  now applicable to the Coulomb-driven flow too.
