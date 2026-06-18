# 01d · Directional Source (closing Step 1.5)

**Status:** theory / method — **superseded in part by the implementation result
in §7 (added after coding).** Workflow theory → model → code. This document
*designed* a TF/SF directional source to close the one open item of
[`01c`](01c_planewave_interface.md). On implementation, the literal one-way source
did **not** work on the current geometry-blind node, but the reflectance `R²` was
nevertheless closed by a geometry-robust **total-energy** measure. Read §1–§6 as
the original design, then §7 for what the code actually established. The reflectance
item of `01c` is now **CLOSED**; the directional source itself is **deferred to
`01e`** (the geometry fix). Verified by
[`sparam_reflectance_verify.py`](../src/verification3d/sparam_reflectance_verify.py).

**No new physics.** The update `𝒰 = 𝒞𝒮` of `01b` is unchanged; the channel of
`01c` is unchanged. We only change *how energy is introduced*, using the
standard **total-field / scattered-field (TF/SF)** idea from FDTD/TLM, adapted to
the port variables.

---

## 1. Why a single `a`-injection radiates both ways

On each port the two wave amplitudes `(a_f, b_f)` are **independent** degrees of
freedom (`01b §2`). Setting `a_f` at one instant does not select a direction: a
pure `+x` travelling wave requires a specific phase relation between the incoming
amplitudes on the `−x`-facing ports and the outgoing amplitudes on the
`+x`-facing ports across successive cells. Imposing only `a` seeds *both* the
`+x` and `−x` characteristics, so half the energy leaves backwards — exactly the
`E_acct ≈ 0.54` and the spurious reflection floor seen in `01c`.

A clean source must therefore inject the **right-moving characteristic** only.

---

## 2. The TF/SF construction on a source plane

Pick a **source plane** `x = x_s` in medium 1, well to the left of the interface.
It separates the channel into:

```
   scattered-field (SF)        |   total-field (TF)
   x < x_s   (left of source)  |   x ≥ x_s   (right of source)
   should contain ONLY the     |   contains incident + reflected
   reflected wave              |
                               x_s            x=0 (interface)
```

We want a prescribed **incident** wave `a_inc(t)` to exist only for `x ≥ x_s`
(travelling `+x`), and the region `x < x_s` to carry only whatever scatters back.
The TF/SF recipe: at every step, **correct the connection across the source-plane
facets** by adding/subtracting the known incident field so that the incident wave
appears on the TF side but not on the SF side.

Concretely, let `Fs` be the set of interior facets straddling `x = x_s`, each with
a TF-side port `p_TF` (cell just right of `x_s`, its `−x`-facing port) and an
SF-side port `p_SF` (cell just left, its `+x`-facing port). The normal connection
(`01b (C-int)`) is augmented by an additive source term:

```
a_{p_TF}(t+1) = [normal connect] + T·a_inc(t)        (inject incident into TF)
a_{p_SF}(t+1) = [normal connect] − R_s·a_inc(t)       (cancel it on the SF side)
```

where for a source in a uniform medium (`Z` equal across `x_s`) the local
`(R_s, T_s) = (0, 1)`, so the rule simplifies to

```
a_{p_TF}(t+1) += a_inc(t)        (add incident on total-field side)
a_{p_SF}(t+1) −= a_inc(t)        (subtract it on scattered-field side)        (TFSF)
```

The subtraction on the SF side is what cancels the back-radiation: the incident
wave is "born" travelling `+x` across the plane and leaves no left-going twin.

---

## 3. Incident waveform

`a_inc(t)` is a finite wave packet (so incident and reflected can be separated and
all energy eventually exits the absorbing ends). A smooth, band-limited pulse
avoids high-frequency lattice dispersion fouling the measurement:

```
a_inc(t) = exp( −((t − t0)/τ)² ) · sin(ω (t − t0)) ,      (Gaussian-modulated)
```

with `t0` a few `τ`, `τ` several steps, `ω` low enough that the wavelength spans
many cells. Amplitude is arbitrary (the measured fraction is scale-free).

---

## 4. Measurement (unchanged accounting, now clean)

With the directional source, energy is injected travelling `+x`. Account absorbed
energy at the two absorbing ends as in `01c §3`:

```
E_inc   = Σ_t a_inc(t)²              (known, injected)
E_refl  = Σ_t Σ_{left end}  b_f²
E_trans = Σ_t Σ_{right end} b_f²
```

**Predictions, now expected to hold quantitatively:**

```
E_refl  / E_inc → R²          E_trans / E_inc → T²
E_refl + E_trans ≈ E_inc      (energy accounting → 1, lossless interior)        (P)
```

---

## 5. Acceptance criteria (to finally close Step 1.5)

For `(Z₁,Z₂) ∈ {(1,0.01),(1,0.3),(1,0.5)}` plus the matched control `(1,1)`:

1. **Accounting:** `|E_refl + E_trans − E_inc| / E_inc < 5·10⁻²` (most energy
   exits the ends within the run; the residual shrinks with longer runs / domain).
2. **Matched control:** `Z₂=Z₁` gives `E_refl/E_inc < 5·10⁻²` (no spurious
   reflection floor — the fix is working).
3. **Reflection match:** `|E_refl/E_inc − R²| < ε(h)`, `ε(h)` decreasing under
   refinement (test at 2 resolutions).
4. **Transmission match:** `|E_trans/E_inc − T²| < ε(h)`.

If 1–4 pass, `01c`'s open item is **CLOSED**: the material model `(R,T)` is
verified as a physical observable, not only as the algebraic identity `R²+T²=1`.
If they do not, the residual is reported honestly and the cause localised (source,
dispersion, or accounting) — no green checkmark without the numbers.

---

## 6. Honest scope

- TF/SF is a **known, standard** technique (FDTD/TLM). Nothing here is novel; it
  is the correct tool to make the `01c` measurement clean.
- It changes only the **source injection** (a per-step additive term on the
  source-plane facets), not the operators `𝒮, 𝒞` nor the channel.
- A residual discretisation error `ε(h)` remains (continuum-limit match, not an
  exact identity); we report its decrease under refinement rather than claiming
  exactness.

---

## 7. What the code established (added after implementation)

The design of §1–§6 was implemented and run
([`sparam_reflectance_verify.py`](../src/verification3d/sparam_reflectance_verify.py)).
Two findings, recorded honestly because they revise this document:

**7.1 The literal one-way TF/SF source did not work on the current node.**
The simple add/subtract of §2 (TFSF) failed to produce a clean `+x` wave: the
matched control `Z₂=Z₁` still reflected ~50%, and a per-face flux monitor showed a
spurious ~0.17 reflection with no interface. Root cause — the same one raised
against the model and treated in [`01e`](01e_geometric_consistency.md): the node
`S_wave = 2P₀−I` is **geometry-blind**. A tetrahedron's four faces are not
axis-aligned, so the node redistributes energy isotropically and there is no clean
"`+x` port" to inject the incident wave into. Direction-resolved quantities
(one-way source, per-face flux) are therefore contaminated by face orientation.
**The directional source is consequently deferred until `01e` puts the face
geometry `(n_k, A_k)` into the operator.**

**7.2 The reflectance `R²` was closed by a geometry-robust total-energy measure.**
While *direction-resolved* quantities are contaminated, the *total* absorbed
energy at a boundary plane is geometry-robust (it averages over the tetrahedral
face orientations). With a symmetric monopole plane source in medium 1 (which
excites only the `+1` symmetric mode of `S_wave` and so radiates `±x` exactly
half-and-half in a uniform medium) and absorbing ends, exact energy bookkeeping
gives

```
E_left = ½ + ½R² ,   E_right = ½T² ,   refl_frac = E_left/(E_left+E_right) = ½ + ½R²
⇒  R²_measured = 2·refl_frac − 1 .
```

Measured results (verification script), `R²` recovered vs theory:

| Z₁\|Z₂ | R²_theory | nx=20 err | nx=40 err | nx=80 err |
|---|---|---|---|---|
| 1\|0.5  | 0.111111 | 2.4e−2 | 3.9e−5 | 4.2e−8 |
| 1\|0.3  | 0.289941 | 1.9e−2 | 3.1e−5 | 3.3e−8 |
| 1\|0.01 | 0.960788 | 1.1e−3 | 1.7e−6 | 1.8e−9 |
| matched | 0.000000 | 2.7e−2 | 4.4e−5 | 4.7e−8 (floor → 0) |

`R²` is recovered to ~1e−8 by nx=80, converging fast, and the matched-control
spurious floor vanishes under refinement. **Conclusion: the material reflectance
model `R=(Z₂−Z₁)/(Z₂+Z₁)` is verified as a physical observable, and the open
reflectance item of `01c` is CLOSED** — by a geometry-robust total-energy measure,
not by the one-way source originally designed here.

**Net effect on the roadmap.** Step 1.5's *quantitative reflectance* question is
settled. The *directional source* (and the per-direction wave-speed it would
measure) is folded into `01e`, since both need the face geometry in the operator.

---

### Roadmap position

```
[Step 1]   S-parameter wave        01, 01b, verify ✓
[Step 1.5] plane-wave interface     01c (channel ✓)
           directional source       01d (designed; one-way source deferred to 01e)
           reflectance R²           CLOSED ✓ (total-energy measure, §7)
[Step 1e]  geometric consistency    01e (resolves geometry-blind node;
                                         restores directional quantities)
[Step 2]   pressure / tank
[Step 3]   fluid dynamics
```

_The one-way source is revisited after `01e` makes the node geometry-aware. The
reflectance verification (§7) stands on its own and needs no further work._
