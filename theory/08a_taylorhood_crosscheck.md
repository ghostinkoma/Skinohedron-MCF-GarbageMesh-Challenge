# 08a · Taylor–Hood Cross-Check — the Cavity Match Is Not a Tuning Artifact

**Status:** **verified**, `src/verification3d/cavity_taylorhood_verify.py`. A direct
follow-up to `08`. P0 validated the lid-driven cavity with **PSPG-stabilised** equal-order
P1 — but PSPG carries a stabilisation parameter `τ`, so a sceptic can ask whether the
Ghia match was *bought by tuning `τ`*. This note removes that doubt by re-solving the same
cavity with a genuine **inf-sup-stable** element that has **no parameter at all**.

---

## 1. Why this check, and in this order

The honest weakness of `08` was the free parameter `τ`. The disciplined way to retire it
is **independent replication**: solve the identical problem with a method that shares
*none* of PSPG's tuning freedom, and see whether it lands on the same answer. The
**Taylor–Hood** element (P2 velocity / P1 pressure) is exactly that — inf-sup (LBB)
stable by construction, **no stabilisation term**. Doing this *before* pushing to higher
Reynolds number is deliberate: if PSPG and Taylor–Hood disagreed even at `Re=100`, any
high-`Re` discrepancy could not be diagnosed. Fix the foundation, then raise the load.

---

## 2. Result — two independent elements, one benchmark

**Taylor–Hood (no `τ`) vs Ghia (1982), Re=100:**

| `n` | centreline u RMS vs Ghia | `u(0.5)` (Ghia −0.2058) |
|----|------|------|
| 16 | 0.0096 | −0.1859 |
| 24 | 0.0062 | −0.1935 |

The parameter-free element converges to the published benchmark — and reaches the same
accuracy at `n=24` (P2) that PSPG reached near `n=64` (P1), because P2 is the velocity
enrichment of `06f`.

**Head-to-head, PSPG vs Taylor–Hood (the decisive line):**

> PSPG (`n=64`, `τ=0.002`) vs Taylor–Hood (`n=24`, no `τ`): centreline RMS **0.0015**.

The two independent discretisations agree with **each other** more tightly (`0.0015`)
than either agrees with Ghia's tabulated points (`~0.006`, which carry their own
interpolation error). So the Ghia match is **not** an artifact of choosing `τ`: a
stabilised method and a parameter-free inf-sup-stable method converge to the *same*
flow. P0 is doubly confirmed.

---

## 3. Consistency with the one-operator philosophy

Taylor–Hood enriches the **velocity** space to P2 while keeping pressure P1 on the same
vertices. That is precisely the **local high-order enrichment of `06f`** ("enrich `K`,
keep the low-order companion") — here applied to make the pressure inf-sup stable. So the
cross-check is not a departure from the single-operator program; it is that program's own
enrichment route, used to earn a trustworthy pressure two independent ways.

---

## 4. What is and isn't claimed

**Verified:**
- Taylor–Hood P2-P1 (no stabilisation parameter) converges to Ghia `Re=100`
  (RMS `0.0062` at `n=24`).
- PSPG and Taylor–Hood agree with each other to RMS `0.0015` — the benchmark match is
  parameter-independent, i.e. correct physics.

**Not claimed:**
- Still **only `Re=100`**. The cross-check licenses the move to higher `Re`, it does not
  perform it (next: `Re=400/1000`, where corner vortices test both elements harder).
- Taylor–Hood here is **2-D**; the 3-D Taylor–Hood (P2 tetrahedra) is a heavier build,
  not done.
- No new mathematics: Taylor–Hood and the Ghia benchmark are classical; the contribution
  is retiring the `τ` doubt by independent replication, inside this project's discipline.

---

## 5. Files

- Theory: this note, on `08`/`06f`/`06z`.
- Verification: `src/verification3d/cavity_taylorhood_verify.py` (checks A–B; imports the
  PSPG solver from `cavity_pspg_verify.py` for the head-to-head).
- Next rung: `Re=400` and `Re=1000` with **both** elements (consistent cross-check at
  higher Reynolds number).
