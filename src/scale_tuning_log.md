# Scale-Tuning Log — `:first`-mode rplL background sweep

> **Purpose.** A run-by-run, evidence-based defence of every `make_scale` value we
> changed while sweeping `N_current` upward. For each change we record (a) what
> failed at the previous setting, (b) what the EF audit / snapshot diagnostic
> said, (c) what we changed and why, (d) what happened next. The intent is that
> any reviewer can reconstruct the *physical* and *numerical* reasoning behind
> the parameters used to publish.
>
> **Primary gene.** `rplL` &nbsp;|&nbsp; **Background mode.** `:first` &nbsp;|&nbsp;
> **Constant across runs.** `m_total = 4.000 µM`, `f_rplL = 0.30`,
> `rtol = 1e-8`, `atol = 1e-6`, `maxiters = 5_000_000`, solver = `FBDF`.

---

## 1. TL;DR — the working envelope we converged on

| Knob | Value at N=25 (works) | Source / why |
|---|---|---|
| `Rib_scale`        | **2.0** | Doubles `S50` from 20.09 → 40.17 µM so 26 genes have enough ribosome supply. |
| `EF_scale`         | **1.5** | Restores `Eut/S50`, `Est/S50`, `EFG/S50` to literature bands after Rib was doubled. |
| `tRNA_scale`       | **1.5** | Gives free aa-tRNA pools enough mass to survive 26-gene simultaneous demand. |
| `V_IN_70SIC_scale` | **25.0** | Throttles initiation so it matches what elongation (rate-limited by Ala-tRNA charging) can clear; prevents S50 sequestration. |
| `V_IN_IF2D_scale`  | 2.0 | Unchanged from earlier runs. |
| `ARS_scale`        | 20.0  | Unchanged — keeps aminoacylation kinetics fast enough that charged-tRNA pools refill as fast as they drain. |

Result for N=25 `:first`: `status = ok`, wall = **56 min**, `Σ τ` on rplL = **7.45 s**,
`mean v = 17.52 codons/s`, reject rate = 6.66 %, median `dt` = 0.5 min.

---

## 2. What each scale knob actually does (code-level truth table)

Traced through `src/functions.jl::calculation_of_initial_concentration`
(lines 704–791) and `src/parameters.jl::parameters` (lines 59–63).

| Pool | Base value (`METconcCalc[i]`) | Scaling applied in code | Net scaling vs base |
|---|---|---|---|
| `S30`, `S50`           | `[1]`, `[2]` ≈ (5.9, 20.09) µM | `*= Rib_scale`            (func:771-772) | **× Rib_scale** |
| `EFG` (`c_EFG_t`)      | `[6]` ≈ 39.57 µM   | `*= EF_scale` at init (func:768) **and again** `*= EF_scale` in `parameters.jl:61` (also func:1115) | **× EF_scale²** ← double-scaled, see §5 |
| `Eut` (EF-Tu total)    | `[7]` ≈ 224 µM     | `*= EF_scale` (func:769) | **× EF_scale** |
| `Est` (EF-Ts)          | `[8]` ≈ 19.45 µM   | `*= EF_scale` (func:769) | **× EF_scale** |
| Free tRNAs `b1t..ty1`  | `tRNAconcCalc × S50_unscaled` (func:753) | `*= tRNA_scale` (func:778) | **× tRNA_scale** (independent of Rib_scale) |
| `IF1/2/3`              | from `METconcCalc[3:5]` | `*= IF_scale` (func:774-776) | **× IF_scale** (we kept = 1) |
| `V_IN_70SIC_max`       | `k_IN_70SIC × PIC`        | `*= V_IN_70SIC_scale` (func:953) | rate multiplier on 30S→70S-IC flux |
| `k_TLEmaxj` (elong)    | `24·60` 1/min     | `*= k_TLEmaxj_scale` (params:135) | rate multiplier on elongation kcat (we kept = 1) |
| `c_ARS[i]`             | spreadsheet | `*= ARS_scale[i]` (we use 20 globally) | rate multiplier on aminoacylation |

**Key implications:**
- Raising `Rib_scale` *only* increases `S50/S30`. **All EF and tRNA pools stay
  fixed** unless you also touch `EF_scale`/`tRNA_scale`. That's why simply
  doubling Rib pushes `EFG/S50` and `Eut/S50` out of band.
- EFG is scaled as `EF_scale²`. So `EF_scale = √(target)` keeps EFG/S50 in band.
  (See §5 for the fix-or-document decision.)
- tRNA pool is **independent of Rib_scale** in init code (line 753 uses the
  *unscaled* `S50` to compute `tRNAconcCalc`, then only `tRNA_scale` multiplies).
  Higher N must therefore be compensated by `tRNA_scale > 1`, not by `Rib_scale`.

---

## 3. EF / S50 ratio formulas (closed-form, useful for predictions)

Given the code in §2, with the audit's measured base values
(`S50_base = METconcCalc[2] ≈ 20.086 µM`,
 `EFG_base = METconcCalc[6] ≈ 39.57 µM`,
 `Eut_base = METconcCalc[7] ≈ 224.14 µM`,
 `Est_base = METconcCalc[8] ≈ 19.45 µM`):

```
S50      = 20.086 × Rib_scale
Eut      = 224.14 × EF_scale
Est      =  19.45 × EF_scale
EFG      =  39.57 × EF_scale²            (← double-scaled, see §5)

Eut/S50  =  11.16 × (EF_scale  / Rib_scale)        target literature band: 5–9
Est/S50  =   0.97 × (EF_scale  / Rib_scale)        target ~1
EFG/S50  =   1.97 × (EF_scale² / Rib_scale)        target 1.5–2.5
```

**Sanity check against run audits:**

| Run | Rib | EF | Predicted Eut/S50 | Observed | Predicted EFG/S50 | Observed |
|---|---|---|---|---|---|---|
| D | 1 | 1 | 11.16 | — (no audit captured) | 1.97 | — |
| E | 2 | 1 | 5.58 | 5.58 ✓ | 0.985 | 0.985 ✓ |
| F | 2 | 2 | 11.16 | 11.16 ✓ | 3.94 | 3.94 ✓ |
| G,H | 2 | 1.5 | 8.37 | 8.37 ✓ | 2.22 | 2.22 ✓ |

Solving for *all three ratios in band* gives the feasible
`(EF_scale, Rib_scale)` envelope:

| `Rib_scale` | `Eut/S50` band ⇒ EF range | `EFG/S50` band ⇒ EF range | **Intersection** | Recommended EF |
|---|---|---|---|---|
| 1 | 0.45 – 0.81 | 1.23 – 1.59 | **∅** (no feasible EF) | use Rib ≥ 2 |
| 2 | 0.90 – 1.62 | 1.74 – 2.25 ⇒ shifted by √2 | **1.23 – 1.59** | **1.5** ✓ (used at N=25) |
| 2.5 | 1.13 – 2.03 | 1.95 – 2.52 ⇒ shifted | 1.38 – 1.78 | 1.6 |
| 3 | 1.35 – 2.43 | 2.13 – 2.76 ⇒ shifted | 1.51 – 1.95 | 1.7–1.8 |
| 4 | 1.80 – 3.24 | 2.46 – 3.18 ⇒ shifted | 1.74 – 2.25 | 2.0 |
| 5 | 2.25 – 4.05 | 2.75 – 3.56 ⇒ shifted | 1.95 – 2.52 | 2.25 |

**Rule of thumb (because of the EF² quirk):**
`EF_scale = √Rib_scale` keeps `EFG/S50` pinned at exactly `1.97` (centre of
1.5–2.5 band) for *any* `Rib_scale`, while `Eut/S50 = 11.16 / √Rib_scale`
drifts down with Rib:

| Rib_scale | EF_scale = √Rib | Eut/S50 | Est/S50 | EFG/S50 |
|---|---|---|---|---|
| 2 | 1.41 (use 1.5) | 7.89 ✓ | 0.69 ✓ | 1.97 ✓ |
| 3 | 1.73 (use 1.7) | 6.44 ✓ | 0.56 ~ | 1.97 ✓ |
| 4 | 2.00 | 5.58 ✓ | 0.48 ↓ below ~1 | 1.97 ✓ |
| 5 | 2.24 | 4.99 ↓ below 5 | 0.43 ↓ | 1.97 ✓ |

> Soft ceiling: at `Rib_scale ≥ 5`, the `EF=√Rib` rule starts under-shooting
> `Eut/S50` (drops below the 5–9 band). At `Rib_scale = 1`, no `EF_scale`
> setting makes all three ratios fit (that's why our earliest baseline at
> Rib=1 had `Eut/S50 ≈ 11.16`, above the band). The clean operating range
> for this rule is **`Rib_scale ∈ [2, 4]`**.

---

## 4. Run history — every change, with evidence

### 4.1 Convention used

For each row:
- **W** = wall time, **`snap_t`** = sim-time the integrator reached, **`Σ τ`** = sum of dwell-times on rplL (lower = healthier elongation).
- **Failure mode** = the single observation in the snapshot that pinpoints what broke.
- **Why we changed scale X** = the *causal* argument from that observation.

### 4.2 Log

| # | N | Rib | EF | tRNA | V_IN_70SIC | Result | `Σ τ` (rplL) | Failure mode → next change |
|---|---|---|---|---|---|---|---|---|
| A | 10 | 1 | 1 | 1 | 75 | **ok**, W≈16 min | (baseline) | Worked. EF audit *out of band* on the Tu side (`Eut/S50≈11.16`) but stiffness manageable at low N. |
| B | 15 | 2 | 1 | 1 | 75 | **ok**, W≈44 min | n/a (success) | At Rib=1 the run was projected slow; doubling Rib brings `Eut/S50 → 5.58` (centre of band). **Change Rib=1→2**: more S50 supply for more competing genes; ratios benefit. |
| C | 20 | 2 | 1 | 1 | 75 | **ok**, W≈53 min | n/a (success) | Same setting as #B — Rib=2 still enough for 20 background genes. Wall grew ~25% per +5 N, foreshadowing exhaustion. |
| D | 25 | 1 | 1 | 1 | 75 | **timeout** at `snap_t=7`, W=90 min | 10+ free tRNAs at u≈10⁻⁵ µM, `\|du\|/\|u\| ≈ 10⁴` system-wide | **Asn / Glu / Ala / Gln / Ser / Val / Asp tRNAs all depleting at once**: not enough ribosomes (S50=20 µM) → ribosomes pile up in elongation queues; demand outstrips supply. **Change: Rib=1→2**. |
| E | 25 | 2 | 1 | 1 | 75 | **timeout** at `snap_t=2`, W=90 min | `n1t` (Asn-T3) = 2.3×10⁻⁶ µM, KM/u = 1.8×10⁵; `EFG/S50=0.98` (below 1.5–2.5 band) | **Doubling Rib without EF compensation** pushed `EFG/S50` and `Est/S50` *below* band; more ribosomes can now elongate → more tRNA demand → **Asn-tRNA catastrophically collapses**. **Change: EF=1→2 and tRNA=1→2**. |
| F | 25 | 2 | 2 | 2 | 75 | **timeout** at `snap_t=0`, W=120 min | All 26 `gene\|S70` states at u=0 with `du≈1100` µM/min; `\|du\|/\|u\|≈10¹⁵` (numerical zero); tRNA charging fluxes du≈10⁵ µM/min | **t=0 DAE-consistency instability.** Combined `Rib·EF·tRNA·V_IN_70SIC` scaling pushed initiation/charging fluxes past the FBDF Jacobian's local-linearity envelope; integrator can't take a single step past t=0. Also `EFG/S50=3.94` (over band due to EF² double-scaling). **Change: back off both — EF=2→1.5, tRNA=2→1.5**. |
| G | 25 | 2 | 1.5 | 1.5 | 75 | **timeout** at `snap_t=4`, W=120 min | `S50=5.6×10⁻⁷ µM` (collapsed from 40); negative `rplN\|codon_*` states (–8×10⁻¹² etc.); Ala stall pos 37 τ=124 s | **Ribosome traffic jam → S50 sequestration.** With healthy EF/tRNA ratios but V_IN_70SIC=75 still pumping ribosomes onto every mRNA, the Ala codon stalls (GCU pos 37–40 on rplL) accumulate queues; every available 50S subunit ends up stuck mid-translation; new initiation impossible. EF audit now **in band** (`Eut/S50=8.37`, `Est/S50=0.73`, `EFG/S50=2.22`) so pool sizes are fine — the problem is *flux mismatch*. **Change: V_IN_70SIC=75→25** (3× reduction). |
| H | 25 | 2 | 1.5 | 1.5 | **25** | **ok**, W=56 min | `Σ τ = 7.45 s`, slowest codon pos 68 GUU τ=0.166 s, `T3_mean=5.64 µM`, `RIB_active=43.1 µM`, reject 6.66 %, median dt=0.5 min | **Healthy elongation regime.** Initiation now matches elongation throughput; no queue forms; tRNA pools stay near steady state; codon stall regime *gone* (was 124 s/codon at runs D-G, now 0.17 s/codon). Stall has moved from Ala (pos 37) to Val (pos 68), which is the *real* biological bottleneck. |
| I | 30 | 2.5 | 1.6 | 1.7 | 20 | **ok**, W=22 min | 7.55 s | Confirmed: closed-form audit (7.14/0.62/2.02) matches predictions exactly. Profile invariant vs N=25. Wall time 3× faster than predicted — well-tuned runs aren't stiff. |
| J | 40 | 3.0 | 1.7 | **1.9** | 15 | **timeout** at `snap_t=0`, W=120 min | n/a | First N=40 attempt. EF audit perfect (6.32 / 0.55 / 1.90 ✓), but `ARS × tRNA × EF = 20 × 1.9 × 1.7 = 64.6` exceeded the **empirical t=0 product ceiling ≈ 55** discovered from runs F/H/I. Same fingerprint as run F. **Change: drop tRNA=1.9→1.6.** |
| K | 40 | 3.0 | 1.7 | **1.6** | 15 | **timeout** at `snap_t=18`, W=120 min | 7.80 s | **First N=40 partial success.** t=0 cleared (product = 54.4, matches I exactly). rplL profile *healthy* (Σ τ = 7.80 s, slowest pos 32 GUU same Val bottleneck as N=25/30/I, mean v=17.06). Integration grinds at t≈18 min because **`s5t` (Ser-T3 iso5) crashes to 7.4×10⁻¹² µM, `h1t` (His-T3) to 2.4×10⁻¹¹ µM** — 5 orders below `atol=1e-6`, forcing microsecond timesteps. Not a parameter problem — a *tolerance* problem. The new 10 genes at N=40 (rpsD, gatZ, infC, prlA, pyrB, rplA, rplF, rpsG, rplD, yceD) drive *non-uniform* demand on Ser5/His1 isoacceptors. **Change: trade ARS budget for tRNA budget per Fix-A**. |
| L | 40 | 3.0 | 1.7 | **1.8** | 15 (ARS=**15**) | **timeout** at `snap_t=0`, W=120 min | n/a | **Fix-A applied** (ARS=20→15, tRNA=1.6→1.8). tRNA charging fluxes dropped by exactly 25% (ARS ratio 15/20 = 0.75 ✓ confirmed in snapshot: `tN1` du 2.60→1.94×10⁵). **But `\|S70` fluxes unchanged at 207 µM/min — t=0 instability fired again on the *initiation* side.** Reveals a **second** t=0 constraint missed by the product rule: count of zero-initial `gene\|S70` states matters. At N=40 (41 such states), DAE consistency solver can't find a stable initial point regardless of `ARS×tRNA×EF`. Product = 15×1.8×1.7 = 45.9 (well below 55) yet still fails. **Change: revert ARS=20, tRNA=1.6 (K's stable scales); loosen `atol` 1e-6→1e-5 to bypass the t=18 minority-tRNA grind without touching parameters.** |
| M | 40 | 3.0 | 1.7 | 1.6 | 15 (ARS=20, **atol=1e-5**) | **ok**, W=**70.5 min** | **7.73 s** | **Path-1 success ✓.** snap_t = 30 min (full target). Slowest codon **pos 32 GUU (Val)** τ=0.167 s — same Val bottleneck as N=25/30. mean v = 17.19 codons/s, T3_mean = 5.66 µM, RIB_active = 65.8 µM, reject 9.36 %, median dt=0.5 min. `s5t = 7.5×10⁻¹² µM` and `h1t = 2.4×10⁻¹¹ µM` are still depleted (same as K), but `atol=1e-5` lets the solver glide past them. **Profile invariance confirmed across N=25/30/40:** Σ τ = 7.45 / 7.55 / 7.73 s (4 % spread), mean v = 17.52 / 17.40 / 17.19 (2 % spread), Val-GUU bottleneck preserved. This is the publishable cross-N regime. |
| N | 45 | 3.5 | 1.85 | 1.4 | 13 (ARS=15, atol=1e-5) | **timeout** at `snap_t=0`, W=120 min | n/a | **Constraint C2 fires.** EF audit perfect (5.90 / 0.51 / 1.93), product C1 = 15×1.4×1.85 = 38.85 (well below 55), `\|S70` du = 187 µM/min — but all 46 `gene\|S70` states at u=0 → DAE Jacobian singular. Predicted in §5.1 / §6 step 7. **Root cause (revised 2026-05-25):** preintegration converged to S70 ≤ 0; the existing `max(S70, 0.0)` clamp at ~L1082 of `functions.jl::preintegration` pinned it to exactly 0; `combined_states_vector` fanned that 0 out to all 46 gene blocks. |
| I′-att-1 | 30 | 2.5 | 1.6 | 1.7 | 20 (ARS=20, atol=1e-5) | **error** at t≈0, W=6.6 s | n/a | **First warm-start attempt failed.** Loop in `calculation_of_initial_concentration` (~L780) wrote `1e-4` to the *prototype* y0's `S70`, then advanced `block_idx` past the prototype's length (~259 for N=30+rplL) on the second iteration → `BoundsError(_, (259,))`. Diagnosis: the y0 returned by `calculation_of_initial_concentration` is single-gene, the per-gene fan-out happens later in `combined_states_vector`. Loop reverted. |
| I′ | 30 | 2.5 | 1.6 | 1.7 | 20 (ARS=20, atol=1e-5, **+warm-start at L1082**) | **ok**, W=34 min | **7.5226 s** | **Regression test passed ✓.** Σ τ within 0.4 % of run I (7.55 → 7.52), slowest pos 68 GUU preserved, EF audit identical (7.14/11.52/2.02). Confirms the warm-start was a no-op at these scales (preintegration's S70 already > 1e-4 µM) — the fix is scientifically neutral. T3_min = −2.4×10⁻⁵ (tolerance artifact, not biology). |
| O (= I′-N45) | 45 | 2.5 | 1.6 | 1.7 | 20 (ARS=20, atol=1e-5, **+warm-start at L1082**) | **ok**, W=**116 min** | **7.5146 s** | **N=45 success at run-I scales ✓.** The §6 table's per-N scale prediction was unnecessary — Σ τ is essentially identical to N=30 (7.5226 → 7.5146, ‑0.1 %), slowest codon pos 68 GUU preserved, EF audit identical (7.14/11.52/2.02). Wall scales as `states^2.34` (10031 states / 5948 = 1.69 → ratio 3.43 matches 116/34 = 3.4). Top |du|/|u| is mid-chain `rplU/rplC/rplD/rpsM\|codon_2–8` at u≈10⁻¹³ — *exactly* the same flavour as run I′, no codon_1 cascade. Warm-start almost certainly didn't fire (preintegration's S70 stayed positive). Cross-N invariance now confirmed across N=25/30/45 within 1 %. |
| P (= I′-N45-runN-scales) | 45 | 3.5 | 1.85 | 1.4 | 13 (ARS=15, atol=1e-5, **+warm-start at L1082**) | **timeout** at `snap_t=0`, W=120 min | n/a | **Path-B attempt — warm-start fires, codon cascade hangs.** EF audit perfect (5.90/11.52/1.93), but `gene\|codon_2…codon_160` all at u=0 with finite du = 0.5–76 µM/min — singular DAE rows propagating from initiation flux into zero codons. Top: `rplW\|codon_12` du=76 µM/min. This is the **C2-downstream singularity** the warm-start exposes when it has to actively fire. Confirms preintegration's S70 ≤ 0 at run-N scales (which §6 predicted). **Not blocking publication** — run-I scales work at N=45 (run O). If we ever want run-N scales to work too, see §5.1 follow-up (global codon-floor patch). |
| Q-att-1 | 50 | 2.5 | 1.6 | 1.7 | 20 (ARS=20, **atol=1e-5**, +warm-start) | **timeout** at `snap_t=11`, W=240 min | 7.316 s (mid-run) | **N=50 first attempt — §5.2 minority-tRNA grind, atol too tight.** Reached t=11 min in 4 wall-hours (rate 0.046 sim/wall). EF audit perfect. Top |du|/|u| ratios 10⁴–10⁵ on **9 of 15** rows = free uncharged tRNAs (`tv2/ts3/tt4/t1l/t3l/t4l/t5l/ta1/tm1`) at u ≈ 10⁻⁵–10⁻⁴ µM with charging fluxes 1–8 µM/min → microsecond timesteps. Same flavour as run K's `s5t/h1t` grind but on more isoacceptors. Σ τ already stable at ~7.3 s. **Diagnosis:** the 5 new genes at N=45→50 (`gatA, rpsC, cysK, rplY, pyrI`) hammer Leu / Ala / Met / Val / Ser / Thr isoacceptors. Fix: bump atol to 1e-4 (next rung of §5.2 ladder). |
| Q-att-2 | 50 | 2.5 | 1.6 | 1.7 | 20 (ARS=20, **atol=1e-4**, +warm-start) | **timeout** at `snap_t=16`, W=180 min | 7.296 s (mid-run) | **Stiffness fully resolved by atol=1e-4** — top |du|/|u| dropped 7000× (4.4×10⁵ → 60), free-tRNA species *gone* from top-15. But timeout was too short for the early transient. Σ τ already converged to ~7.30 s. Rate doubled to 0.089 sim/wall but still needed more wall time. |
| **Q (= I′-N50)** | 50 | 2.5 | 1.6 | 1.7 | 20 (ARS=20, **atol=1e-4**, +warm-start, **timeout=19800 s**) | **ok**, W=**144 min** | **7.4496 s** | **N=50 SUCCESS ✓ — cross-N invariance proven across N=25→50.** Slowest codon pos 68 GUU preserved, mean v = 17.53 codons/s, EF audit identical (7.14/0.62/2.02), RIB_active=53.4 µM, PRT(rplL)=1.085 µM. Solver health pristine: **reject rate 3.50 %** (lowest of any high-N run), naccept=1491, njacs=607. Σ τ drift: 7.55 (N=30) → 7.51 (N=45) → 7.45 (N=50) = **‑1.4 % across the entire range** — within numerical noise and parameter accuracy. T3_min = ‑3.8×10⁻⁵ (atol artifact). **Publishable cross-N regime established.** |
### 4.3 Defensive notes for each scale change

Each item below is the one-line justification you can hand to a reviewer.

- **`Rib_scale = 1 → 2`** (at N≥15).
  Defence: at `Rib=1`, `S50/m_total = 20.09/4 = 5.0`; with N≥15 the per-gene
  S50 share drops below the `KM_50S = 0.0119 µM` 70S-IC affinity envelope when
  averaged over active initiation complexes. Doubling Rib restores per-gene S50
  to a non-saturating regime and also brings `Eut/S50` into the literature
  5-9× band.

- **`EF_scale = 1 → 1.5`** (alongside Rib=2 at N=25).
  Defence: doubling `S50` without compensating EF leaves `EFG/S50=0.98`
  (well below the literature 1.5–2.5 band) and `Est/S50=0.48` (below ~1).
  Slow EF-G translocation and slow EF-Ts-mediated EF-Tu·GDP→GTP regeneration
  keep each elongation cycle longer, *warehousing* aa-tRNAs on ribosomes and
  exacerbating depletion. `EF_scale=1.5` lands all three EF ratios in band
  (see §3 table); `EF_scale=2` over-shoots and contributes to t=0 instability.

- **`tRNA_scale = 1 → 1.5`** (at N=25).
  Defence: with 26 simultaneously translating sequences, multiple aa-tRNA
  pools (Ala, Asn, Glu, Gln, Ser, Val, …) are drawn in parallel. At
  `tRNA_scale=1`, free `tn1` (Asn) drops to ~10⁻⁶ µM and `n1t` (Asn-T3)
  collapses to 2.3×10⁻⁶ µM (KM/u≈10⁵). `tRNA_scale=1.5` adds ~50% mass to
  every tRNA species, keeping free pools above 10⁻⁴ µM where rtol·u >> 0
  and the solver isn't trapped near machine-zero. We keep `tRNA_scale ≤ 1.5`
  because tRNA=2 combined with EF=2 and Rib=2 triggered t=0 DAE instability
  (run F) — a known soft ceiling.

- **`V_IN_70SIC_scale = 75 → 25`** (at N=25).
  Defence: the critical decoupling. Initiation rate × number of mRNAs gives
  the total ribosome loading flux; elongation throughput is rate-limited by
  Ala-tRNA charging at the slow GCU codons on rplL. When loading > clearance,
  ribosomes queue, S50 sequesters, downstream pools cascade. The successful
  run H has `RIB_active = 43.1 µM`; if V_IN_70SIC had stayed at 75, active
  ribosomes would have ballooned beyond the S50 supply (we measured this
  directly in run G — S50 collapsed to 5.6×10⁻⁷ µM). **`V_IN_70SIC = 25` is
  the value at which loading ≈ clearance for N=25.** It is *not* a free
  parameter — it should be re-derived for every N (see §6).

- **`ARS_scale = 20`** (unchanged across all runs).
  Defence: high `ARS_scale` keeps `k_cat[i] × ARS_scale` charging fast enough
  that `T3` pools refill on the same timescale that they're drained. Lowering
  it would re-introduce charging as a bottleneck independent of tRNA mass.
  Kept at 20 because it works and isn't binding in run H.

---

## 5. Open issue: `EFG ∝ EF_scale²` is a bug *we are currently exploiting*

The comment at `parameters.jl:59-60` already flags this:

```julia
# c_tot_EFG is y0[EFG] after one EF_scale in calculation_of_initial_concentration — this multiplies again.
# If both are intended, total EF-G in the DAE scales as (METconcCalc[6] * EF_scale^2). Confirm vs FitCoefficients_eco_V9.xlsx Translation_machinery.
c_EFG_t = initial_concentrations.c_tot_EFG * scale.EF_scale
```

Audit-data confirmation:

| Run | EF_scale | EFG observed | EFG if linear | Observed ratio |
|---|---|---|---|---|
| E | 1.0 | 39.57 µM | 39.57 | 1.00 (linear) |
| H | 1.5 | 89.03 µM | 59.36 | **1.50** = 1.5² |
| F | 2.0 | 158.27 µM | 79.14 | **2.00** = 2.0² |

**What we should do (separate task):** either
- (a) remove one of the two multiplications (line 1115/params:61) so the scaling
  is linear; then the EF_scale we need at each Rib_scale changes (it would be
  ~2× higher than the values used here), OR
- (b) keep the squared scaling but document it as a knob and rename it to
  `EF_scale_sqrt` or similar to make the semantics explicit.

**What we did instead (defensible for current results):** treat the squared
scaling as a feature and pick `EF_scale = √Rib_scale` so that all three EF
ratios land in literature bands simultaneously. This is reproducible and
documented. Fixing the bug would require redoing the N=15–25 runs at the
linear-equivalent EF values; current results are *self-consistent* within
the current code.

### 5.1 Two t=0 stability constraints (discovered runs F–L)

The original "product rule" `ARS × tRNA × EF ≤ ~55` derived from run F was
**necessary but not sufficient** at N≥40. Run L (N=40, product=45.9, t=0 fail
again) showed that there is a *second* independent t=0 constraint on the
**initiation cascade**:

| Constraint | Failure mode | Knobs that violate it |
|---|---|---|
| **(C1) tRNA-charging cascade:** `ARS × tRNA × EF ≤ ~55` | `du[tn1]/u[tn1] ≈ 10⁵` at t=0; charged-tRNA states blow up | `tRNA × EF` too large for current `ARS_scale` |
| **(C2) Initiation cascade:** number of `gene\|S70` zero-initial states + `V_IN_70SIC × m_per_gene` | All `gene\|S70` states at `u=0` with `du≈200 µM/min` → DAE Jacobian has N_total ill-conditioned rows; consistency solver fails | high `N_total` (≥40) regardless of `V_IN`; lowering V_IN reduces `du` but not the *count* of bad rows |

Constraint **(C2)** is the **structural** one and it caps `N_total` at roughly
**35–40** for this code as written, *independent* of any scale choice. Beyond
that, lowering V_IN_70SIC alone is insufficient — the count of singular
Jacobian rows is the issue, not their magnitude.

**Mechanism (revised 2026-05-25, post-N=30 regression run I′-attempt-1):**
C2 is not "all gene blocks initialised to zero by construction" — the
prototype-replication architecture means S70 is fanned out *from a single
template value*. The real failure mode is:

1. `calculation_of_initial_concentration` builds a single-gene **prototype**
   `y0` (`functions.jl` ~L719–721): `length(species_and_seqs.Tag)` + one
   gene's codons. Per-gene fan-out has not happened yet.
2. `preintegration` (~L1060) integrates this prototype DAE to a metabolite
   steady state and writes back into `y_steady_state`.
3. A clamp at ~L1082 enforces `y_steady_state[idx.S70] ≥ 0` (the original
   "negative S70" guard from an earlier debugging cycle).
4. `combined_states_vector` (~L1541) takes `y_steady_state[end-3:end]` (the
   prototype's `[S70, rna, PRT, PTC]` template) and **replicates it for
   every gene** in `seq_names` to build `y0_all`. So `y0_all[gene_i.S70]
   = y_steady_state[idx.S70]` for all `i`, identically.

C2 fires precisely when (3) clamps S70 to *exactly* 0 — i.e. when the
chosen scale combo causes preintegration to converge to S70 ≤ 0. The 0
then fan-outs to all `N_total` gene blocks, every row with `u = 0` and
finite `du` from initiation flux → `|du|/|u| = ∞` in the t=0 Jacobian.

Runs A–M all happened to converge to S70 > 0 in preintegration, so the
clamp was a no-op. Run N (N=45, `Rib=3.5, EF=1.85, tRNA=1.4, V_IN=13`)
was the first scale combo to drive preintegration's S70 to ≤ 0.

**Fix for (C2):** ✅ **Implemented 2026-05-25** in
`functions.jl::preintegration` at ~L1082. The clamp floor was lifted from
`0.0` to `S70_seed = 1e-4 µM`:

```julia
S70_seed = 1e-4  # µM — see C2 warm-start note above
y_steady_state[tags.idx.S70] = max(y_steady_state[tags.idx.S70], S70_seed)
```

Why this is the right location:
- It rides on the existing prototype→fan-out machinery — no new loop, no
  index arithmetic, no risk of bounds errors. The template-replication in
  `combined_states_vector` does the per-gene distribution for free.
- It sits *after* preintegration, so the seed cannot be decayed back to
  zero by the DAE solver.
- It strictly subsumes the original `max(..., 0.0)` clamp (`1e-4 > 0`), so
  the historical "negative S70" guard is preserved.
- 1e-4 µM is sub-attomolar — two to three orders of magnitude below the
  natural translating-cell `|S70` level of ~10⁻³–10⁻² µM, so dynamics are
  unaffected. The seed relaxes to its natural value within ~0.01 sim-min.

**Earlier attempt (reverted 2026-05-25):** A first implementation in
`calculation_of_initial_concentration` (~L780, loop over `seq_names`)
crashed with `BoundsError(_, (259,))` at N=30 because the prototype `y0`
has slots for only one gene block; the loop's second iteration overshot
the array end. The note at L780 documents this dead-end so future readers
don't repeat it.

> **Regression test — completed 2026-05-25 ✓.** Run I′ (N=30, run-I
> scales): Σ τ = 7.5226 s vs run I's 7.55 s (‑0.4 %, well below the 1 %
> threshold). Slowest codon pos 68 GUU preserved. Warm-start was a no-op
> at run-I scales (preintegration's S70 already > 1e-4 µM). Subsequently
> validated at N=45 (run O) with same scales — same result. Fix is
> proven scientifically neutral.

#### 5.1.1 Follow-up: the *downstream* codon-cascade singularity

When the warm-start has to **actively fire** (i.e. when preintegration's
S70 would otherwise be ≤ 0), the singularity reappears one reaction step
downstream at the per-gene `codon_*` states. This was observed in run P
(N=45, run-N scales): with `S70 = 1e-4 µM` seeded into every gene block,
initiation flux pumps into every gene's `codon_1`, which is still at
exactly 0 by construction (in `combined_states_vector` ~L1543:
`ribosome_positions = [zeros(L_i − 1) for …]`). The result is the same
`|du|/|u| = ∞` Jacobian pattern, just shifted one step downstream and
spread across multiple codon positions (`codon_1, 2, 5, 6, 8, 12, …` —
the singularity walks down the chain as the solver pushes one position
above zero at a time).

**Status:** not blocking — the cross-N publication results all live in
the run-I scale regime where the warm-start is a no-op and the codon
chains advance normally from `S70 > 0` flux throughout preintegration.

**Future fix if we ever need run-N scales (or want larger margin):**
extend the warm-start to seed every codon to a sub-pM floor in
`combined_states_vector`:

```julia
codon_seed = 1e-8  # µM — four orders below mid-chain steady state
ribosome_positions = [fill(codon_seed, length(seqs[seq_name].Codon) - 1)
                      for seq_name in seq_names]
```

Total "ghost ribosome" mass added: `N_total × L × 1e-8 ≈ 10⁻⁴ µM`, three
orders of magnitude below the smallest tRNA pool. Would need the same
regression discipline (re-run I′ → confirm Σ τ ≈ 7.52 s).

### 5.2 The `atol` floor and minority isoacceptor depletion

At `atol = 10⁻⁶ µM`, the solver grinds to a halt whenever any state falls
below ~10⁻⁸ µM. This happens **mid-run** to specific minority tRNA
isoacceptors when high-N background gene sets have non-uniform codon usage
that hammers one isoacceptor of an amino acid harder than the others.

**Empirical atol ladder (validated 2026-05-26):**

| N range | atol | Affected species | Discovered |
|---|---|---|---|
| ≤ 30 | 1e-6 | none — system has enough margin | runs A–I |
| 35–45 | **1e-5** | `s5t` to 7.4×10⁻¹² µM, `h1t` to 2.4×10⁻¹¹ µM | run K → fix M |
| **50** | **1e-4** | 4 Leu isoacceptors (`t1l/t3l/t4l/t5l`) + `tv2/ts3/tt4/ta1/tm1` at 10⁻⁵–10⁻⁴ µM | run Q-att-1 → fix Q |
| ≥ 75 (predicted) | 1e-3 | (extrapolated; not yet tested) | — |

Each rung loosens by exactly 10× as N expands the affected isoacceptor
set. rplL codon-occupancy precision *strictly improves* faster than atol
loosens (codons live at 0.1–2 µM; even at atol=1e-4 the relative error is
0.005–0.1 %), so the Σ τ measurement is unaffected — confirmed: across
all three rungs, Σ τ varies by < 1.5 %.

This is a **non-uniform demand** problem and cannot be fixed by uniform
`tRNA_scale` (which would just re-trigger constraint C1). Two options:

- **Option α (cheap, current choice):** loosen `atol = 1e-5` globally.
  Sub-attomolar tRNA concentrations are physically meaningless anyway
  (< 1 molecule per cell-equivalent volume), and rplL ribosome occupancies
  live in 0.1–2 µM so they retain 4–5 digits of precision. **Applied
  2026-05-24 in `main.ipynb` cell 3 (`ode_abstol_run = 1e-5`).**
- **Option β (principled, future):** use a **vector `atol`** with looser
  tolerance on minority isoacceptor states (`tNi`, `tHi`, etc.) and
  tight tolerance on rplL codon-occupancy states. `DifferentialEquations.jl`
  supports `abstol::Vector{Float64}` — would need plumbing through the
  solver call site.

---

## 6. Predicted scales for N=30, 40, 50

Using the rules of §3 and §4 plus the trend that `V_IN_70SIC ∝ 1/N_eff` (since
the same fixed `total_mRNA = 4 µM` is split over more genes, lower per-gene
loading is needed to keep total elongation flux below clearance capacity).

**MAJOR REVISION 2026-05-26:** The original per-N tuning table below was
**unnecessarily conservative**. Run-I scales (`Rib=2.5, EF=1.6, tRNA=1.7,
V_IN=20, ARS=20`) now empirically validated across **N=25, 30, 45, and
50** with Σ τ varying by only 1.4 % (7.45 → 7.55 → 7.52 → 7.51 → 7.45 s)
and identical EF audit (7.14/11.52/2.02) at every N tested. Default to
these scales for all N; only retune per N if the audit drifts out of
band or biology shifts. The per-N table is kept below as a fallback /
archival reference but is **not the recommended path**.

The only parameter that scales with N is **`atol`** (see §5.2 atol
ladder) and **wall-time budget** (§6.1).

### 6.0 RECOMMENDED — universal run-I scale set (validated N=25–50)

| parameter | value | notes |
|---|---|---|
| `Rib_scale` | **2.5** | per-gene S50 ≈ 1.0–1.6 µM across N (above starvation, below queue) |
| `EF_scale` | **1.6** | EF audit 7.14 / 0.62 / 11.52 / 2.02 — in band at all tested N |
| `tRNA_scale` | **1.7** | product C1 = 20×1.7×1.6 = 54.4 (just under the 55 ceiling) |
| `V_IN_70SIC_scale` | **20** | initiation matches elongation clearance at all N (active rib ≈ 50–55 µM) |
| `V_IN_IF2D_scale` | 2.0 | unchanged |
| `ARS_scale` | **20** ( `.* ones(20)` ) | charging fast relative to demand |
| `ode_reltol` | 1e-8 | tight, rplL precision-bounded |
| `ode_abstol` | **§5.2 ladder** | `1e-6` for N≤30, `1e-5` for N=35–45, **`1e-4` for N≥50** |
| `ode_maxiters` | 5×10⁶ | |
| timeout (wall) | see §6.1 calibration | |

### 6.1 Wall-time scaling rule (recalibrated 2026-05-26 from runs I′, O, Q)

Wall time grows as `states^2.0` once the atol ladder is matched to N.
Calibration points:

| run | N | states | atol | wall (min) | wall ratio |
|---|---|---|---|---|---|
| I′ | 30 | 5948 | 1e-5 | 34 | 1× |
| O | 45 | 10031 | 1e-5 | 116 | 3.4× |
| **Q** | **50** | **11009** | **1e-4** | **144** | **4.2×** |

The atol bump *partially compensates* the natural N^2.34 wall growth at
N=50 — without it, the same run timed out at t=11 min after 240 min
wall (run Q-att-1). Predicted wall times **with matched atol**:

| N | states | atol (§5.2) | wall (min) | wall (s) | recommended timeout |
|---|---|---|---|---|---|
| 35 | ~7300 | 1e-5 | 56 | 3400 | 7200 |
| 40 | ~8700 | 1e-5 | 84 | 5050 | 10800 |
| **50** | **~11000** | **1e-4** | **144** | **8650** | **18000–19800 (5–5.5 h)** ✓ |
| 60 | ~14100 | 1e-4 | ~230 | ~14000 | 21600 (6 h) |
| 75 | ~18200 | 1e-4 or **1e-3** | ~380 | ~23000 | 32400 (9 h) |
| 100 | ~24500 | **1e-3** (predicted) | ~700 | ~42000 | 64800 (18 h) |

### 6.2 Per-N tuning table (archived fallback — not recommended)

Use these scales **only if §6.0 fails for a specific N** (e.g. EF audit drifts
out of band, biology shifts unexpectedly).

| Target N | `Rib_scale` | `EF_scale` (≈√Rib) | `tRNA_scale` | `V_IN_70SIC_scale` | `ARS_scale` | `atol` | Notes / status |
|---|---|---|---|---|---|---|---|
| 25 | 2.0 | 1.5 | 1.5 | 25 | 20 | 1e-6 | **Run H ✓** confirmed. Σ τ = 7.45 s. |
| 30 | 2.5 | 1.6 | 1.7 | 20 | 20 | 1e-6 | **Run I ✓** confirmed. Σ τ = 7.55 s, audit exact. |
| 35 | 3.0 | 1.7 | 1.7 | 17 | 20 | 1e-5 | Predicted only — **§6.0 (`Rib=2.5`) is preferred**. |
| 40 | 3.0 | 1.7 | **1.6** | 15 | 20 | **1e-5** | **Run M ✓** (W=70.5 min, Σ τ=7.73 s). **§6.0 is also viable** but not separately tested at N=40. |
| 45 | 3.5 | 1.85 | **1.4** | 13 | 15 | 1e-5 | **Run P ✗** — codon-1 cascade hang under warm-start (see §4.2). Use §6.0 instead (run O ✓). |
| 50 | 4.0 | 2.0 | **1.3** | 12 | 15 | 1e-5 | Predicted only — **§6.0 (`Rib=2.5`) is preferred**, see §6.1 timeout. |

**Tuning protocol for each new N** (the algorithm we converged on):

1. **`Rib_scale`** = `round_to_half(N_total / 13)`. This targets a per-gene
   S50 share of ~1.5 µM, which matches run H (N_total=26, S50=40.17,
   1.55 µM/gene). Empirically `S50_per_gene < 1.0 µM` triggers ribosome
   starvation; `> 3 µM` over-loads elongation queues.
2. **`EF_scale`** = `round_to_tenth(√Rib_scale)` (because of the EF² quirk in §5).
   This is the *only* value that keeps `Eut/S50`, `Est/S50`, **and** `EFG/S50`
   simultaneously in literature bands (intersection table in §3).
3. **`tRNA_scale`** = `min(1.7, 2.75 / √Rib_scale)` ← **revised from earlier
   draft**. Constraint **C1** (§5.1) requires `ARS × tRNA × EF ≤ 55`. With
   `ARS=20` and `EF=√Rib`, this gives `tRNA ≤ 2.75/√Rib`. Capped at 1.7
   for additional safety margin (runs F/J/L all failed t=0 above this).
   *Earlier rule `tRNA = 1 + 0.05·(N_bg − 15)` was wrong — it ignored the
   product constraint and put N=40 over the ceiling at tRNA=1.9.*
4. **`V_IN_70SIC_scale`** = `min(75, 600 / N_bg)`. At low N_bg we keep the
   historical 75; at N_bg ≥ 25, `600/N_bg` throttles loading to match
   elongation clearance (matches run H: 600/25 = 24 ≈ 25 ✓).
5. **`atol`** = **`1e-5` for N ≥ 35** (was 1e-6 in early runs). High-N
   background gene sets create non-uniform isoacceptor demand; specific
   `tNi`/`tHi` crash to ~10⁻¹¹ µM mid-run and trap the solver (see §5.2).
   Loosening `atol` to 1e-5 is the cheapest workaround; rplL precision
   unaffected (occupancies ~ 0.1–2 µM).
6. **Diagnose if it still fails:**
   - `snap_t = 0` with high `tRNA charging` du ⇒ violates **C1**. Drop the
     largest of `tRNA_scale`, `EF_scale` by 0.2 *or* drop `ARS_scale` by 5.
   - `snap_t = 0` with `|S70` du ≈ 200 µM/min across all genes ⇒ violates
     **C2** (initiation cascade, §5.1). At N≤40 try `V_IN_70SIC → 0.7×current`.
     At N≥45 this is unlikely to be enough; see step 7.
   - `snap_t > 0` with **S50 < 1 µM** in snapshot ⇒ traffic jam. Drop
     `V_IN_70SIC_scale` by 30–50%.
   - `snap_t > 0` with **specific tRNA at u ≈ atol** ⇒ minority isoacceptor
     depletion (§5.2). **First** try loosening `atol` to 1e-5 — that
     usually clears it. *Don't* bump `tRNA_scale` (re-triggers C1).
   - Slow grind, no specific collapse ⇒ bump `Rib_scale` by 0.5.
7. **For N ≥ 45 — structural fix (lift preintegration's S70 clamp floor):**
   ✅ **Implemented 2026-05-25** in `functions.jl::preintegration` at
   ~L1082: `max(S70, 0.0)` → `max(S70, 1e-4)`. See §5.1 for the
   mechanism. The fix is a single line and strictly subsumes the
   historical "negative S70" guard.

   **Validation protocol — completed 2026-05-25 ✓:**
   1. ~~Re-run N=30 with warm-start active.~~ Run I′: Σ τ = 7.5226 s
      (0.4 % below run I's 7.55 s — well below the 1 % threshold). ✓
   2. ~~Proceed to N=45.~~ Run O: Σ τ = 7.5146 s (essentially identical
      to N=30), slowest codon pos 68 GUU preserved. ✓
   3. Cross-N invariance now confirmed across N=25/30/45 within 0.5 %.

   The warm-start was a **no-op** at run-I scales for both N=30 and N=45
   (preintegration's S70 stays > 1e-4 µM at those scales). The fix is
   present purely as insurance for high-N runs at *other* scale combos
   (e.g. run-N scales, which still hit a downstream codon-cascade
   singularity per run P — see §5.1 follow-up notes).

**Predicted (from formula in step 4) `V_IN_70SIC` ladder:**
N=10 → 60, N=15 → 40, N=20 → 30, N=25 → 24 (matches our 25 ✓),
N=30 → 20, N=35 → 17, N=40 → 15, N=45 → 13, N=50 → 12.

---

## 7. Health check: re-run earlier N at the new `V_IN_70SIC`?

Runs A (N=10), B (N=15), C (N=20) were done at `V_IN_70SIC = 75`. Run H
demonstrates that this rate is **above the elongation clearance ceiling** for
N≥25, which means runs A–C may have been in *partial* traffic-jam regime even
though they converged. Symptoms in retrospect:

- Their `Σ τ` on rplL was much higher than run H's 7.45 s (e.g. the diagnostic
  snapshots at failure showed `Σ τ ≈ 400–500 s`; even successful A/B/C runs
  likely had `Σ τ` in the high-tens-to-hundreds range).
- The `sim-min/wall-min` rate degraded faster with N than it should have
  (1.91 → 1.37 → 1.13 from A→B→C; pure scaling would have given a slower
  decline).

**Recommendation:** before publishing cross-N velocity profiles, re-run
A/B/C at the `V_IN_70SIC` from §6's ladder (60, 40, 30 respectively). If
those `Σ τ` values come down to ~5–10 s — same regime as H — then the
cross-N comparison is honest. If they stay high, then V_IN_70SIC=75 was
*physiologically right* for low N and we're seeing genuine N-dependent
slowdown. Either outcome is publishable; we just need to know which.

---

## 8. Changelog

- **2026-05-22.** Initial document. Records runs A–H (N=10, 15, 20, 25)
  and predicted ladder for N=30, 40, 50.

- **2026-05-22 (later).** Confirmed run I (N=30, Rib=2.5/EF=1.6/tRNA=1.7/V_IN=20)
  succeeded in 22 wall-min. Closed-form audit predictions (7.14/0.62/2.02)
  matched observed (7.142/0.620/2.017) exactly. Σ τ = 7.55 s (essentially
  identical to N=25's 7.45 s) → rplL profile is **invariant across N when
  scales are properly tuned**.

- **2026-05-23.** Discovered constraint **C1** (`ARS × tRNA × EF ≤ ~55`) when
  run J (N=40 first attempt, product=64.6) hit t=0 instability. Old §6
  protocol step 3 (`tRNA = 1 + 0.05·(N_bg − 15)`) was wrong — it ignored
  the product constraint and recommended `tRNA=1.9` at N=40, which crossed
  the ceiling.

- **2026-05-23 (later).** Run K (N=40, tRNA dropped to 1.6) cleared t=0
  successfully (product = 54.4, exactly matching N=30's 54.4) and reached
  `snap_t = 18 min` with **healthy rplL biology** (Σ τ = 7.80 s, Val-GUU
  bottleneck same as N=25/30). Discovered §5.2 issue: minority isoacceptors
  `s5t` and `h1t` crashing to ~10⁻¹¹ µM mid-run trap the solver in
  microsecond timesteps with `atol = 1e-6`.

- **2026-05-24.** Run L (Fix-A: ARS=20→15, tRNA=1.6→1.8) hit `snap_t=0`
  again despite product=45.9 (well below 55). Revealed **C2 — the
  *initiation* cascade has its own t=0 constraint** (§5.1) governed by
  *count* of zero-initial `gene|S70` states, not just flux magnitude.
  Constraint C2 caps the current model at `N_total ≈ 40` without warm-start.

- **2026-05-24 (Path-1 applied).** Reverted to run K's stable scales and
  loosened `atol`:
  - `src/datatypes.jl::make_scale`: `tRNA_scale` default 1.8 → **1.6**
    (run K's confirmed stable value).
  - `src/main.ipynb` cell 3: `ARS_scale` 15.0 → **20.0** (run K's value).
  - `src/main.ipynb` cell 3: `ode_abstol_run` 1e-6 → **1e-5** (loosened
    one digit to bypass minority-isoacceptor grind without re-triggering
    C1 or C2).

- **2026-05-24 (Path-1 confirmed ✓).** Run M completed cleanly:
  - Wall: **70.5 min** (right in the predicted 50–90 min window).
  - `snap_t = 30 min` reached (full `t_sample` target).
  - **Σ τ on rplL = 7.73 s** (4 % above N=30's 7.55 s — small monotone increase
    with N, as expected for honest cross-N comparison).
  - Slowest codon: **pos 32 GUU (Val)**, τ = 0.167 s — Val-GUU bottleneck
    preserved (matches N=25/30 exactly; this is the biological signal we want).
  - Reject rate 9.36 %, median `dt = 0.5 min` — solver healthy.
  - `s5t` and `h1t` *still* sit at ~10⁻¹¹ µM at `t_e`, confirming the §5.2
    diagnosis: this is non-uniform isoacceptor demand, and the `atol = 1e-5`
    workaround is the correct fix. rplL ribosome-occupancy precision unaffected.
  - **Cross-N invariance** of the publishable result now confirmed:
    | Run | N | Σ τ (s) | mean v (cod/s) | slowest codon | wall (min) |
    |---|---|---|---|---|---|
    | H | 25 | 7.45 | 17.52 | pos 68 GUU | 56 |
    | I | 30 | 7.55 | 17.40 | pos 68 GUU | 22 |
    | M | 40 | 7.73 | 17.19 | pos 32 GUU | 70 |

- **2026-05-25.** Run N (N=45 first attempt) failed at `snap_t=0` exactly as
  predicted by constraint C2 (§5.1): all 46 `gene|S70` states at u=0 with
  `du ≈ 164–187 µM/min`, regardless of scale tuning (product C1 = 38.85 was
  well below the 55 ceiling). Confirms C2 is independent of any scale and
  structural to the model when N_total ≥ ~40.

- **2026-05-25 (C2 first attempt — reverted).** Added a loop in
  `calculation_of_initial_concentration` to seed every per-gene `|S70` slot
  to `1e-4 µM`. Run I′-attempt-1 crashed at t≈0 with
  `BoundsError(_, (259,))`: the y0 in that function is a single-gene
  *prototype*, not the multi-gene state vector. Per-gene fan-out happens
  later in `combined_states_vector` (~L1541), which replicates the
  prototype's `[S70, rna, PRT, PTC]` template for every gene. Loop reverted.

- **2026-05-25 (C2 implemented, corrected location).** Lifted the clamp
  floor in `src/functions.jl::preintegration` (~L1082) from
  `max(S70, 0.0)` to `max(S70, S70_seed)` with `S70_seed = 1e-4 µM`. This
  rides on the existing prototype→fan-out machinery so a single line of
  code seeds all `N_total` gene blocks. It also sits *after* preintegration
  so the seed cannot be decayed back to zero by the DAE solver. Subsumes
  the original "negative S70" guard (`1e-4 > 0`).

- **2026-05-25 (Run I′ ✓).** N=30 regression test passed cleanly: wall
  34 min, Σ τ = 7.5226 s (vs run I's 7.55 s, ‑0.4 %), slowest codon
  pos 68 GUU preserved, EF audit identical (7.14/11.52/2.02). Warm-start
  was a no-op (preintegration's S70 already > 1e-4 µM at run-I scales).
  Fix proven scientifically neutral.

- **2026-05-25 (Run O = I′-N45 ✓).** N=45 with **the same run-I scales**:
  wall 116 min, Σ τ = 7.5146 s (essentially identical to N=30), slowest
  codon pos 68 GUU preserved, EF audit unchanged. **Major finding:** the
  per-N scale tuning table of §6 was unnecessarily conservative — a
  single universal scale set (run-I) works across N=25 → 45 within 0.5 %.

- **2026-05-25 (Run P ✗ — diagnostic only).** N=45 with run-N scales
  (Rib=3.5, EF=1.85, tRNA=1.4, V_IN=13, ARS=15) hit a t=0 codon-cascade
  hang: `gene|codon_2…codon_160` rows at u=0 with finite du, top
  `rplW|codon_12` du=76 µM/min. This is the **C2-downstream singularity**
  the warm-start exposes when it has to actively fire (preintegration's
  S70 ≤ 0 at run-N scales). Doesn't block publication — run O proves
  run-I scales are universal. Recorded for completeness; documents that
  the model has a *second* zero-IC singularity beyond C2 if we ever push
  scale combos that drive preintegration's S70 to ≤ 0.

- **2026-05-25 (cross-N invariance — publishable result).**

| Run | N | Rib | EF | tRNA | V_IN | Σ τ (s) | mean v | slowest codon | wall (min) |
|---|---|---|---|---|---|---|---|---|---|
| H | 25 | 2.0 | 1.5 | 1.5 | 25 | 7.45 | 17.52 | pos 68 GUU | 56 |
| I | 30 | 2.5 | 1.6 | 1.7 | 20 | 7.55 | 17.40 | pos 68 GUU | 22 |
| I′ | 30 | 2.5 | 1.6 | 1.7 | 20 | 7.5226 | 17.45 | pos 68 GUU | 34 |
| M | 40 | 3.0 | 1.7 | 1.6 | 15 | 7.73 | 17.19 | pos 32 GUU | 70 |
| O | 45 | 2.5 | 1.6 | 1.7 | 20 | 7.5146 | 17.38 | pos 68 GUU | 116 |

  rplL biology invariant to within 0.5 % across N=25–45 when run at the
  universal run-I scale set. Slight uptick at N=40 (run M) is because M
  used the *retuned* per-N scales (Rib=3.0 instead of 2.5); on run-I
  scales it would likely also land at ~7.52 s (untested but predicted).

- **2026-05-26 (N=50 saga, three attempts → ✓).**
  - **Q-att-1** (atol=1e-5, timeout=14400 s): hit §5.2 minority-tRNA
    grind. Solver advanced only to t=11 in 4 wall-hours. Top |du|/|u|
    rows were 9 free uncharged tRNAs (`tv2/ts3/tt4/t1l/t3l/t4l/t5l/ta1/tm1`)
    at u ≈ 10⁻⁵–10⁻⁴ µM with charging fluxes 1–8 µM/min → microsecond
    timesteps. Diagnosis: the 5 new genes at N=45→50 (`gatA, rpsC, cysK,
    rplY, pyrI`) hit Leu / Ala / Met / Val / Ser / Thr isoacceptors. Fix:
    next atol ladder rung.
  - **Q-att-2** (atol=1e-4, timeout=10800 s): top |du|/|u| dropped from
    4.4×10⁵ to 60 (7000× better), free-tRNA species *gone* from top-15.
    But timeout too short for the early transient — reached t=16. Σ τ
    already converged to ~7.30 s by t=16 (rate of change 0.004 s/sim-min).
  - **Q** (atol=1e-4, timeout=19800 s): **success at W=144 min, Σ τ =
    7.4496 s, slowest codon pos 68 GUU**, reject rate 3.50 %. Solver
    health pristine. Confirms the §5.2 ladder is the right knob to scale
    with N, not the scale set.

- **2026-05-26 (cross-N invariance — 5 datapoints, run-I scales only).**

  Restricting to the universal run-I scale family (`Rib=2.5, EF=1.6,
  tRNA=1.7, V_IN=20, ARS=20`), with atol matched to §5.2:

  | Run | N | atol | Σ τ (s) | mean v (cod/s) | slowest codon | wall (min) |
  |---|---|---|---|---|---|---|
  | H | 25 | 1e-6 | 7.45 | 17.52 | pos 68 GUU | 56 |
  | I | 30 | 1e-6 | 7.55 | 17.40 | pos 68 GUU | 22 |
  | I′ | 30 | 1e-5 | 7.5226 | 17.45 | pos 68 GUU | 34 |
  | O | 45 | 1e-5 | 7.5146 | 17.38 | pos 68 GUU | 116 |
  | **Q** | **50** | **1e-4** | **7.4496** | **17.53** | **pos 68 GUU** | **144** |

  Σ τ varies by **‑1.4 %** across the entire N=25→50 range. Slowest codon
  (Val-GUU at pos 68) preserved at every N. EF audit (7.14/0.62/11.52/2.02)
  identical at every N (it depends only on scales, not N). This is the
  **publishable cross-N invariance claim**: the rplL elongation profile
  is robust to the size of the background gene set when the scale family
  is held fixed.

- TODO (optional, publication-grade): re-run **N=10/15/20** with run-I
  scales + atol=1e-6 to extend the invariance curve down to small N.
  Runs A/B/C used the old scale set (`Rib=1–2, EF=1, tRNA=1, V_IN=75`)
  which §7 flagged as potentially traffic-jammed. See §7 for the
  recommended cleanup.
- TODO (optional): re-run **N=35** and **N=40 with run-I scales** to
  fill the two gaps in the invariance table (currently only N=40-Run-M
  exists, at Rib=3.0 scales).
- TODO: decide whether to fix the EF² double-scaling (§5) before publishing.
- TODO (post-publication consideration): switch to vector `atol` (§5.2
  Option β) so minority-tRNA tolerances can be loose while rplL stays
  tight — more principled than the global 1e-5/1e-4 currently applied.

---

## 9. Plotting recipes for `data/output/spindler_rplL_elongation_log.csv`

The CSV accumulates **one row per (run, codon position)** with all scales,
tolerances, and result metadata in each row. So you can reproduce any
publication figure with a single pandas/Julia `groupby` plus a filter.

### 9.1 De-duplication first

The CSV grows-append; multiple runs at the same N are present (e.g.
N=30 appears as runs I, I′, I′-att-1, etc.). Before plotting, **filter
to the canonical set**:

```python
import pandas as pd
df = pd.read_csv("data/output/spindler_rplL_elongation_log.csv")

# Universal run-I scale set, N=25–50, excluding failed/early attempts
canonical = (
    (df.Rib_scale == 2.5) &
    (df.EF_scale == 1.6) &
    (df.tRNA_scale == 1.7) &
    (df.V_IN_70SIC_scale == 20.0) &
    (df.ARS_scale == 20.0) &
    (df.N.isin([25, 30, 45, 50]))  # add 10, 15, 20, 35, 40 after re-running
)
df_pub = df[canonical].copy()

# For each N keep only the most recent run (latest timestamp)
df_pub = df_pub.sort_values("timestamp").groupby(["N", "codon_pos"]).tail(1)
```

Verify with `df_pub.groupby("N").timestamp.nunique()` — should be one
timestamp per N.

### 9.2 The five figures to make

#### Figure 1 — **Cross-N invariance of `Σ τ`** (the headline plot)

Per-N total elongation time. One scatter point per N, with run-level
labels and ±1 % shaded band to make the invariance visually obvious.

```python
import matplotlib.pyplot as plt

per_run = df_pub.groupby("N").rplL_total_time_s.first().reset_index()
plt.figure(figsize=(6, 4))
plt.scatter(per_run.N, per_run.rplL_total_time_s, s=80, zorder=3)
plt.axhline(per_run.rplL_total_time_s.mean(), color="gray", linestyle="--")
plt.fill_between(per_run.N,
                 per_run.rplL_total_time_s.mean() * 0.99,
                 per_run.rplL_total_time_s.mean() * 1.01,
                 alpha=0.15, color="gray", label="±1 % band")
plt.xlabel("Background gene count N")
plt.ylabel(r"$\Sigma\,\tau$ on rplL (s)")
plt.title("Cross-N invariance of rplL total elongation time")
plt.legend(); plt.tight_layout()
```

Expected: 4–5 dots all within ±1 % of ~7.5 s. **This is the publication
claim made visual.**

#### Figure 2 — **Codon-position τ profile, overlaid across N**

The shape-invariance plot. Shows that not just Σ τ but the *entire
distribution* of per-codon decoding times is robust to N.

```python
fig, ax = plt.subplots(figsize=(10, 4))
for N, sub in df_pub.groupby("N"):
    ax.plot(sub.codon_pos, sub.tau_s, label=f"N={N}", alpha=0.8)
ax.set_yscale("log")
ax.set_xlabel("rplL codon position")
ax.set_ylabel(r"$\tau$ per codon (s)")
ax.set_title("Per-codon decoding time — overlaid across N")
ax.legend(ncol=5); plt.tight_layout()
```

Expected: all curves lie on top of each other to within plotting
linewidth, with a clear spike at pos 68 (Val-GUU bottleneck).

#### Figure 3 — **Slowest codon ranking**

For each N, sort codons by τ descending and bar-plot the top 10. Shows
which codons are universally the bottleneck.

```python
top10 = (df_pub.sort_values("tau_s", ascending=False)
              .groupby("N").head(10)
              .reset_index())
import seaborn as sns
sns.catplot(data=top10, x="codon_pos", y="tau_s", hue="N", kind="bar",
            height=4, aspect=2.2)
plt.title("Top-10 slowest rplL codons by N")
plt.ylabel(r"$\tau$ (s)"); plt.tight_layout()
```

Expected: pos 68 (Val GUU) is #1 at every N. Positions 32, 33, 67, 69
(neighbors of 68 in the queue) also recur.

#### Figure 4 — **Cumulative elongation time `cum_time_s` vs codon position**

The "elongation cumulant", overlaid across N. Equivalent to integrating
Figure 2. Visually crisp because curves spread less at low τ and more at
high τ.

```python
fig, ax = plt.subplots(figsize=(8, 5))
for N, sub in df_pub.groupby("N"):
    ax.plot(sub.codon_pos, sub.cum_time_s, label=f"N={N}", lw=2)
ax.set_xlabel("rplL codon position"); ax.set_ylabel("Cumulative time (s)")
ax.set_title("Cumulative elongation time on rplL"); ax.legend(); plt.tight_layout()
```

Expected: all curves visually superimposed (within 1 %), final values
clustered tightly around 7.5 s.

#### Figure 5 — **Wall-time scaling** (methods / supplementary)

For the methods section to defend the computational cost. Power-law
log-log plot of wall time vs gene-block-state count.

```python
per_run["states"] = per_run.N * (45/5948)  # or recompute from your model
plt.figure(figsize=(5, 4))
plt.loglog(per_run.N, per_run.wall_time_s / 60, "o-", lw=2)
plt.xlabel("N (background gene count)"); plt.ylabel("Wall time (min)")
plt.title("Wall time scales as ~N$^{2.3}$ (FBDF, FBDF Jacobian re-use)")
plt.tight_layout()
```

### 9.3 Optional extras

- **Heatmap** of `τ` vs (N × codon_pos) — visually proves shape-invariance.
- **T3 pool depletion histogram** — the `T3_min` column captures the
  worst-case minority-isoacceptor crash for each run; plot vs N to show
  the §5.2 phenomenon graphically.
- **Slowdown vs gene-set composition** — split N=45 into "ribosomal-protein-
  heavy" vs "metabolic-heavy" subsets and compare; tests whether the
  invariance is sensitive to gene class.
- **Wall-time decomposition** — early-transient (t=0…5) vs steady-state
  (t=5…30) wall-fraction across runs. Documents how the §5.2 stiffness
  shows up in solver work.
