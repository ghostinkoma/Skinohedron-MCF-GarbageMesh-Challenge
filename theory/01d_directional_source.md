# 01d · Directional Source (closing Step 1.5)

**Status:** theory / method (no code yet). Workflow theory → model → code. This
is the thin addition that closes the one open item of
[`01c`](01c_planewave_interface.md): the simple port injection was not a clean
`+x` travelling wave (spurious back-radiation, energy accounting < 1), so the
reflected-energy fraction did not match `R²`. Here we specify a **directional
source** that injects a wave moving in `+x` only, so reflection/transmission can
be measured cleanly against `R²/T²`.

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

### Roadmap position

```
[Step 1]   S-parameter wave        01, 01b, verify ✓
[Step 1.5] plane-wave interface     01c (channel ✓, R² open)
           directional source       01d  ← THIS DOCUMENT (theory)
           01d-code                 ← next, on approval: closes R² measurement
[Step 2]   pressure / tank
[Step 3]   fluid dynamics
```

_Next artifact, on approval: code adding the TF/SF source term to the channel and
re-measuring `E_refl/E_inc` vs `R²`, with the §5 acceptance asserts. Reuses the
`01b` engine and `01c` channel unchanged; only the source is added._
