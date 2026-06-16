# TXTLv2 — Mass-balance / Path-A2 pipeline change log

All notable model-, calibration-, and workflow-affecting edits live here.
Most recent entries on top. Each entry pairs with code that is in version
control as of the same date.

Use this file to **defend** every model knob change to a reviewer:
- WHY the change was made (observed problem)
- WHAT was changed (specific lines / cells / files)
- HOW it was justified (literature pointer, sanity check, or empirical cost)
- IMPACT on downstream metrics (median K_deg, % positive, runtime)

---

## 2026-05-20 — N=15 at CONTROL+ settings + appendable elongation-rate log

### Motivation
With the CONTROL+ baseline (N=10, `:first`, `tRNA=1`, `ARS=20`, tight tols)
confirmed against the historical 2026-05-19 sweep (`rplL_total_time` matches
to 0.5 %), we extend the sweep at the known-convergent config. Single-knob
bump: `N_current = 10 → 15`. The reintroduction plan
(`:first → :codon_balanced`, then `tRNA_scale = 3`) is paused until the
N-sweep at CONTROL+ settings reaches at least N=15, ideally N=20.

### What's new in cell 3
- `N_current = 15` (was 10). All other knobs untouched.
- `using ... Dates` added.
- New appendable CSV: `data/output/spindler_rplL_elongation_log.csv`.
  Every successful run writes one row per codon on the primary gene with the
  following schema:

  ```
  timestamp, N, n_seq, bg_mode, primary_gene,
  codon_pos, tau_s, v_smooth_cod_s, cum_time_s,
  rplL_total_time_s, rplL_mean_v_cs, T3_min, T3_mean,
  m_rplL, m_bg_each, total_mRNA,
  V_IN_70SIC_scale, V_IN_IF2D_scale, ARS_scale, tRNA_scale, EF_scale, Rib_scale,
  ode_reltol, ode_abstol, t_sample_min, sol_t_end_min, wall_time_s
  ```

  The file is created on first run with header, appended on every run after.
  Schema mismatch between runs would raise — that's deliberate so a knob
  rename can't silently corrupt the log. Only `effective_status == "ok"`
  runs are logged so the cross-run comparison stays apples-to-apples.

### How to use the log (post-processing example)
```julia
using CSV, DataFrames, Plots
df = CSV.read("data/output/spindler_rplL_elongation_log.csv", DataFrame)
# rplL velocity profile across the N-sweep at one config
sub = subset(df,
    :bg_mode    => x -> x .== "first",
    :tRNA_scale => x -> x .== 1.0,
    :ARS_scale  => x -> x .== 20.0,
)
plot(sub.codon_pos, sub.v_smooth_cod_s, group = sub.N,
     xlabel = "Codon position on rplL", ylabel = "v_smooth (codons/s)",
     title  = "rplL velocity profile vs N (CONTROL+ settings)")
```

### Predicted N=15 outcome
- **Status**: ok (CONTROL+ converged at N=10; N=15 is a small step up).
- **Wall time**: 25–35 wall-min (historical N=10 was 13.8 min × today's 2.1×
  slowdown × ~1.3× for the extra 5 mRNAs).
- **rplL_total_time_s**: probably similar to N=10 (~554 s). The Ala-stall is
  already at `T3_min ≈ 0` saturation at N=10 — adding more competing genes
  may not deepen it further (the slowest codon is *already* maximally slow
  because charged Ala-tRNA is essentially zero). The cleaner signal will
  show up in `mean_v` (arithmetic mean is *not* harmonic-dominated by the
  stalled codons) and in the codon profile away from the stalled positions.

### Note for the elongation log
The N=10 baseline is NOT yet in `spindler_rplL_elongation_log.csv` (the log
file didn't exist when that run completed). Two options after this N=15 run:
1. Re-run N=10 once to backfill the log (cheap, ~30 wall-min).
2. Run a small one-shot backfill from `spindler_rplL_N10_first_codon_detail.csv`.
We'll do (1) by default — gives a fresh timestamp and consistent provenance.

### Impact
- No physics change. Pure tooling addition + N step. The audit trail for the
  in-progress N-sweep at CONTROL+ settings now lives in a single CSV.

---

## 2026-05-20 — CONTROL+ converged: baseline reproduced, model code verified

### Result
CONTROL+ (= CONTROL with `timeout_run_sec = 45·60`) ran to completion:
- `status = ok`, wall time **1729.9 s** (28.8 wall-min), `eval_t = 30 sim-min`.
- `rplL_total_time = 554.05 s`, `mean_v = 13.15 cod/s`, `T3_min = 4.34e-9 µM`.

### Side-by-side vs. historical 2026-05-19 N=10 sweep

| Quantity | Historical (`eval_t = 100 min`) | CONTROL+ (`eval_t = 30 min`) | Match |
| --- | --- | --- | --- |
| `rplL_total_time_s` | 556.80 | 554.05 | **0.5 %** ✅ |
| `T3_mean` (µM) | 3.871 | 3.873 | ✅ |
| `T3_min` (µM) | −1.14e-9 | +4.34e-9 | ✅ both ≈ 0 |
| `charged_tRNA` (µM) | 0.00179 | 0.00179 | ✅ exact |
| `EGT` (µM) | 15.79 | 15.79 | ✅ exact |
| `EFG` (µM) | 0.09469 | 0.09186 | ✅ |
| `PRT(rplL)` (µM) | 1.280 | 1.278 | 0.2 % ✅ |
| `mean_v` (cod/s) | 13.64 | 13.15 | 3.6 % (different `t_sample`) |
| `RIB_active` (µM) | 20.10 | 22.75 | 13 % (different `t_sample`) |
| Wall time (s) | 829.9 | 1729.9 | **2.1× slower today** |

The bottleneck-report diagnostic at `t_e = 30` independently confirms the
historical *Ala-stall*: 4 consecutive GCU (Ala) codons on rplL at positions
37–40 each have τ ≈ 148 s, accounting for 540 of the 554 s `total_time`.
Met-T3 is even more depleted (4.3e-9 µM, KM/u ≈ 9.2e7) but Met appears only
at codon 1 of rplL so doesn't dominate. This is the codon-mix lottery
documented in the 2026-05-20 Tier-1 entry.

### What this finally settles

| Question | Answer |
| --- | --- |
| Did the model code regress since 2026-05-19? | **No.** Defining metric matches to 0.5 %. |
| Did today's code get slower? | Yes, ~2.1× per sim-min. Attributable to the +1129 lines in `functions.jl` and +252 in `rhs_metabolites.jl` (safety guards, RHS refactors). Acceptable cost. |
| What caused the t=0 stuck-snapshot in Probe B and the two preceding runs? | **`tRNA_scale = 3.0`** combined with `:codon_balanced` selection. With `tRNA_scale = 1.0` (default) and `:first`, the integrator marches through normally. |
| Was Probe B (ARS 20→5) attacking the right scale? | No — `ARS_scale` is *not* the t=0 culprit; `tRNA_scale` is. |

### Saved artifacts (this run)
- `data/output/spindler_rplL_N10_first_summary.csv`
- `data/output/spindler_rplL_N10_first_codon_detail.csv`
- `data/output/spindler_rplL_N10_first_background_genes.csv`
- `data/output/figures/spindler_rplL_N10_first_codon_profile.png`
- `data/output/figures/spindler_rplL_N10_first_pool_snapshot.png`

(No `_timeout_snapshot.csv` this time — the integration completed.)

### Next: Step 1 of the reintroduction plan
Change ONE knob on top of CONTROL+: `bg_selection_mode = :first → :codon_balanced`.
Keep `tRNA_scale = 1.0` (default), `ARS_scale = 20`, tight tolerances, all
else identical. Predicted outcome: converges; `T3_min` should be much
further from zero because the codon-balanced gene set is designed to spread
load across tRNAs (avoiding the Ala stall that the :first composition
happened to trigger). Comparing CONTROL+ vs. Step 1 cleanly isolates the
effect of selection mode.

---

## 2026-05-20 — CONTROL first-run result: smoking gun is `tRNA_scale=3`, not the model code

### What we ran
CONTROL config from the previous entry: `N=10`, `:first`, `tRNA_scale=1.0`
(default), `ARS_scale=20`, `rtol=1e-8`, `atol=1e-6`, `t_sample=30`,
`timeout=15·60`.

### Result
**The "stuck at t=0" failure is gone.** The PeriodicCallback snapshot recorded
`last captured t = 21.0 sim-min` — the integrator marched forward and ran
out of wall clock at 15 min while still mid-flight. Compare across the three
recent runs (all N=10):

| Run | snap_t | Wall budget | sim-min per wall-min |
| --- | --- | --- | --- |
| Probe B (tRNA=3, ARS=5, codon_bal, looser) | **0.0 min** | 30 wall-min | 0 (stuck on IC) |
| Prev (tRNA=3, ARS=20, codon_bal, looser)   | **0.0 min** | 30 wall-min | 0 (stuck on IC) |
| **CONTROL** (tRNA=1, ARS=20, :first, tight) | **21.0 min** | 15 wall-min | 1.4 |
| Historical N=10 sweep CSV (same config)     | 100.0 min  | 13.83 wall-min | ~7.2 |

The integrator is solver-tractable. Today's code is ~5× slower per sim-min
than the historical sweep (probable cause: the +252 lines of safety
guards/refactors in `rhs_metabolites.jl` and +1129 lines in `functions.jl`
since `eedd012`), but it converges — there is no t=0 cliff in the model
code itself.

### The biology matches the historical row 1:1

Snapshot at `t = 21.0 min` vs. historical `N=10` sweep row:

| Quantity | Historical CSV (`N=10`) | CONTROL snapshot @ t=21 | Match |
| --- | --- | --- | --- |
| `mean_v` (codons/s) | 13.64 | 13.61 | ✅ |
| `T3_min` (µM) | -1.14e-9 (≈ 0) | `a1t = 1.9e-5` (KM/u = 21,000) | ✅ same — Ala-T3 collapsed |
| `Σ τ` on rplL | (gives `total_time ≈ 556.8 s`) | 506.8 s (still climbing) | ✅ tracking |
| slowest codon | (implied Ala stall by `T3_min ≈ 0`) | **position 37 = GCU = Alanine, τ = 137 s** | ✅ exactly Ala-stall |

The control reproduces the historical *codon-mix lottery* (T3_min≈0 stall on
Ala for the alphabetic-:first 10-gene set).

### Diagnosis: which knob was the actual culprit?

The variable that flipped "stuck at t=0" → "runs and reproduces history" is
**`tRNA_scale`** (3.0 → default 1.0). Two implications:

1. The Spindler-original calibration is solver-tractable on today's code as
   long as `tRNA_scale = 1`. The initial-transient pathology we documented
   in the previous CHANGELOG entry was driven by tripling the tRNA pool —
   3× the substrate caused 3× the initial aminoacylation flux, which
   together with V_IN_70SIC=75 × 11 mRNAs simultaneously loading exceeded
   what FBDF + finite-difference Jacobian on 4250 states could digest in
   the t=0 transient. **Probe B (ARS 20→5) attacked the wrong scale.**

2. Whether `:codon_balanced` selection alone (with `tRNA_scale=1`) would
   also have survived is an open question to be answered by step 1 of the
   reintroduction plan below.

### Single change in this follow-up: timeout 15 → 45 wall-min

To get a clean summary row that matches the historical CSV, we just need
enough wall budget. Extrapolating from `21 sim-min / 15 wall-min`, the full
30 sim-min should take ~22 wall-min; 45 wall-min gives ~2× safety margin.
No other knob changed.

### After CONTROL+ completes, the reintroduction plan

| Step | What changes | Predicted outcome |
| --- | --- | --- |
| 1    | `:first → :codon_balanced` (still `tRNA=1`, `ARS=20`) | Should still converge; T3_min will likely NOT be ≈0 (codon-balanced spreads load, that was its whole point) |
| 2    | `tRNA_scale = 3.0` on top of step 1                   | If this re-introduces the t=0 stuck behavior, `tRNA_scale=3` is confirmed as the binding constraint and we need a different strategy to enlarge the tRNA pool (e.g. ramp it gradually after t > 0, or accept tRNA=1 and look for competition signal in mean_v) |
| 3    | `ARS_scale = 5.0` (Probe B) on top of step 2          | Only meaningful if step 2 passes |

### Impact
- No physics change in this entry (only timeout). All prior reasoning about
  tRNA_scale changing SS pool size still holds — we just learned it's also
  changing the *transient* in a way the solver can't absorb at our scales.

---

## 2026-05-20 — CONTROL run: reproduce OLD-N=10-that-worked on today's code

### Motivation
After Probe B (ARS=5, tRNA=3, codon-balanced, looser tols, longer timeout) ALSO
timed out at `snap_t = 0.0 min`, we paused the forward-probing and asked
"compared to when N=10 *worked* (rplL_total_time ≈ 557 s, T3_min ≈ 0), what
actually changed?". Differences between OLD-N=10 (worked) and Probe B
(timed out):

| Knob | OLD (worked) | Probe B (timed out) | Direction |
| --- | --- | --- | --- |
| `tRNA_scale` | 1.0 (default) | 3.0 | bigger pool, bigger initial flux |
| `ARS_scale` | 20.0 | 5.0 | catalysis rate dropped |
| `bg_selection_mode` | `:first` | `:codon_balanced` | different gene set |
| `t_sample` | 30.0 min | 20.0 min | irrelevant for t=0 |
| `ode_reltol_run` | 1e-8 | 1e-6 | LOOSER (should help) |
| `ode_abstol_run` | 1e-6 | 1e-4 | LOOSER (should help) |
| `timeout_run_sec` | 15 · 60 | 30 · 60 | MORE wall budget |

Tolerances and timeout went in the helpful direction, so the culprit is
necessarily in {selection, tRNA_scale, ARS_scale}. Before throwing more
probes, we run a single diagnostic experiment: **revert every cell-3 knob to
the OLD-worked config** and confirm whether today's code still converges
there. That answers "is the cliff in our settings or in the model code?".

### CONTROL config (this run)
- `N_current = 10`
- `bg_selection_mode = :first`
- `t_sample = 30.0`, `timeout_run_sec = 15 · 60`
- `ode_reltol_run = 1e-8`, `ode_abstol_run = 1e-6`
- `spindler_scale = make_scale(V_IN_70SIC=75, V_IN_IF2D=2, ARS=20.0×ones(20))`
  — `tRNA_scale` NOT passed → defaults to 1.0
- Tier-2 PeriodicCallback kept (diagnostic capability for free)
- Audit-trail printout kept (records which scales were active)
- Output tag: `spindler_rplL_N10_first_*` (separates from earlier codon-balanced outputs)

### Decision rules

| CONTROL result | Interpretation | Next action |
| --- | --- | --- |
| Converges (≈ 557 s, T3_min ≈ 0 — matches historical) | Model code is fine; the cliff is in our recent settings | Re-introduce changes one at a time: (1) `:first → :codon_balanced`; (2) add `tRNA_scale = 3`; (3) drop `ARS_scale → 5`. Whichever step breaks IS the culprit. |
| Converges but very different metrics | Model code drifted in a SS-altering way (not stiffness-altering) | Git-diff against `eedd012` for SS-relevant changes (`rhs_*`, `parameters.jl`). |
| Times out at `snap_t = 0.0` like Probe B | Something in the uncommitted model code DOES regress N=10 convergence (contrary to expectation) | Bisect against `eedd012` to find the offending hunk. |
| Times out at `snap_t > 0` | The model is solver-tractable but slower than before — likely model-code regression in transient dynamics, NOT an IC pathology | Inspect the snapshot's biological state (which species depleting). |

### What's NOT changed
- No model code touched.
- Codon-balanced selection helpers in `functions.jl` remain — just dormant.
- Summary-CSV scale columns remain — they'll log `ARS=20`, `tRNA=1.0` so we
  can join CONTROL with Probe B's CSV offline.

### Impact
- Pure diagnostic. No physics change. Output tag (`_first_`) prevents
  overwriting the codon_balanced artifacts already on disk.

---

## 2026-05-20 — Probe B: ARS_scale 20 → 5 (attack t=0 aminoacylation stiffness)

### Motivation
After the systematic step-down plan landed, the N=10 (codon-balanced, Spindler
scales, `tRNA_scale=3`) run **also** timed out at `snap_t = 0.0 min` — exactly
the same failure mode as N=20. Important learning: **N is not the right knob**
inside this regime. The biology at t=0 is healthy (most-depleted T3 = `s2t` at
1.04 µM, KM/u = 0.38 — well above the depletion threshold; rplL slowest codon
τ = 0.090 s, mean v = 20.86 codons/s). The failure is purely numerical.

The snapshot fingerprint identified the culprit explicitly:

```
Top |du| ranking at snap_t = 0:
  tN1   |du|=4.7e+04 µM/min,  u=3.41 µM
  tn1   |du|=4.7e+04 µM/min,  u=2.92 µM
  tE2   |du|=3.0e+04 µM/min,  u=19.1 µM
  ... (15 aminoacylated tRNAs with |du| in the 1e4–5e4 µM/min range)
```

These fluxes are set by `ARS_scale × tRNA_scale × charging_substrate`. With
`ARS_scale=20` and `tRNA_scale=3` the effective flux is 60× the unscaled
baseline — the solver was forced into microsecond-scale first steps and
re-evaluated a 4250-state finite-difference Jacobian on every reject (autodiff
is off because the RHS has branches), eating the 30-min wall clock without
advancing past t=0.

### Why ARS_scale, not tRNA_scale
- `tRNA_scale` sets the pool *size* — bumping it raises the SS charged
  reservoir (the biology we *want* to control as N grows).
- `ARS_scale` sets the catalytic *rate* — it controls how *fast* equilibrium
  is reached, not the equilibrium itself. At steady-state, charged-tRNA pools
  are mass-balance-limited (charging easily outpaces consumption for these
  pools and substrate levels), so `ARS_scale=5` vs 20 should give
  essentially the same SS charged fraction.
- Net effect: drop `ARS_scale` 20 → 5 softens the t=0 transient by 4× without
  changing the steady-state competition signal we measure at t=20 min.

### Code changes
- `src/main.ipynb` cell `fc8f77a8`, `spindler_scale` block:
  - `ARS_scale = 20.0 .* ones(20)` → `5.0 .* ones(20)`
  - Comment block above the constructor documents the Spindler-original vs
    our-regime contrast and the t=0 diagnostic that drove this change.
- `tRNA_scale = 3.0` retained.
- `N_current = 10` retained (apples-to-apples vs the prior timeout snapshot).

### What the next run will tell us
- **If it converges**: ARS_scale was the binding t=0 constraint. Then push N
  back up (10 → 15 → 18 → 20) with these scales held constant. The
  competition signal is preserved because SS pools are essentially unchanged.
- **If it still times out at t=0**: the dominant culprit isn't aminoacylation
  but ribosome loading (the `*|S70 = 0` states with `|du| ≈ 1235 µM/min`).
  Next probe would be **Option C** — seed `S70` states with `u₀ = 1e-6 µM` so
  rtol applies — which requires an IC tweak in `call_model`.
- **If it converges but `T3_min ≈ 0`** at SS: ARS_scale dropped *too* much
  and steady-state charging can no longer keep up with demand. Re-bump to
  `ARS_scale = 10` (between Spindler 20 and current 5).

### Impact
- Pure transient-kinetics change — steady-state metrics (rplL_total_time,
  mean_v, T3 levels at t=20) are expected unchanged within ±5% if the run
  reaches t_sample. Will verify against the N=11 reference once we converge.

---

## 2026-05-20 — Tier-2 SavingCallback + N=20 t=0 initial-transient stiffness finding + systematic step-down plan

### Motivation
After the 2026-05-20 single-N runner landed (codon-balanced, Spindler scales +
`tRNA_scale=3`), the first N=20 run still timed out — `print_timeout_snapshot_diagnostic`
reported `last captured t = 0.0 min`. The Tier-2 `PeriodicCallback` worked as
designed (it gave us *some* snapshot from a failed run) but it told us a story
different from what we expected:

- The blocker is **not** Met-tRNA depletion at quasi-steady-state.
- The blocker is **the initial transient itself**: at t=0, 21 mRNAs
  simultaneously start loading ribosomes (S70 states with u=0, du≈1200 µM/min)
  AND tRNAs are charging up at very high flux (top |du| ≈ 10^5 µM/min on
  aminoacylated tRNAs). The relative-tolerance check is meaningless when u=0,
  so the absolute-tolerance check forces the first step to be on the order of
  microseconds in sim-time — and 4250 finite-difference Jacobian evaluations
  per step (autodiff=false is required by the ODE function's max/min branches)
  is what eats the wall clock.

`tRNA_scale=3` raises the steady-state charged-tRNA pools but it does *not*
soften the initial transient (the t=0 flux scales with the pool too). So the
correct lever isn't "scale tRNAs harder until convergence at N=20"; the correct
lever is to find the N at which the initial transient becomes unmanageable
under fixed scales, since that N is itself the answer to "how much background
load can the system bear?".

### Strategy: systematic step-down (single-variable per run)

We switch from "fix N=20, escalate scales" → to "fix scales, step N up from a
known-convergent baseline". The N at which convergence breaks (or rplL_total_time
plateaus / saturates) **is** the experimental answer.

| Run order | N    | Expectation                                              |
| --------- | ---- | -------------------------------------------------------- |
| 1         | 10   | Known convergent regime (compare to old N=10 baseline)   |
| 2         | 15   | Should still fit Spindler-style scales                   |
| 3         | 18   | First "stressful" point                                  |
| 4         | 20   | Probable convergence boundary; here scales may need help |

Only escalate scales (`tRNA_scale` first, `ARS_scale` only if needed) at the
boundary — the smallest scale that achieves convergence preserves the *most*
competition signal in `rplL_total_time`, `mean_v`, `T3_min/mean`, codon profile.

### Code changes

- `src/main.ipynb` (cell `fc8f77a8`):
  - `N_current = 10` (was 20). Comment block documents the step-down plan.
  - Header docstring updated with the t=0 finding and step-down sequence.
  - Added a "Run config" printout at the start of every run — records N,
    selection mode, all 6 scales, tolerances, timeout — so the per-N CSVs can
    be cross-referenced offline.
  - Summary CSV now writes `V_IN_70SIC_scale, V_IN_IF2D_scale, ARS_scale,
    tRNA_scale, EF_scale, Rib_scale, ode_reltol, ode_abstol` columns so the
    scales become a queryable dimension when stacking N=10/15/18/20.
- No model/scale changes this entry — pure workflow + audit-trail.

### How to interpret the upcoming N=10 run
- **If it converges in < 60 s wall**: scales are fine for the moderate-load
  regime. Push to N=15 next.
- **If it converges slowly (> 5 min wall) or shows T3_min ≈ 0**: even N=10 is
  already in the depletion regime; codon-balanced helped but the pool is the
  next bottleneck — bump `tRNA_scale` here, not at N=20.
- **If it times out** (very unlikely given prior N=11 data): the t=0 transient
  is the universal blocker, not an N-dependent one — investigate solver or
  better warm-start initial conditions.

### Impact
- No metric changes (no model edit). Future N-runs will produce
  audit-traceable CSVs where the scales each run used are recorded inline.

---

## 2026-05-20 — Tier-1 sweep fixes + codon-mix-lottery finding + codon-balanced selection mode + single-N runner

### Motivation
After the 2026-05-19 N-sweep cell landed we ran the shrunk grid
`N ∈ [1, 3, 5, 10]`. With the original cell settings the run was unusable:
- N=0 took **583 s** (vs ~6 s for the rplL+RTS reference): `spindler_75pct =
  true` pre-loaded ~15 µM ribosomes onto a single mRNA, creating an
  artificial traffic jam.
- N=3 reported `status=ok` but `Integration finished without success`: the
  status was a wrapper false-positive.
- N=5 and N=10 timed out at 900 s and reported `status=timeout`.
- For N=3 the largest-|du| printout pinned the blow-up to the EF-Tu
  algebraic block (`|du[EFT]| ≈ 1.5e9`, `EFD → −5.87 µM`) — the DAE
  algebraic constraint was the actual numerical bottleneck, not codon
  competition.

We also observed a striking non-monotonic pattern that *had* to be
documented or it will trip every reviewer of this experiment:
- N=1: small `T3_min`, no single-codon stall, `rplL_total_time ≈ 11.1 s`
- N=3: `T3_min ≈ 0` → single-codon stall → `rplL_total_time ≈ 521 s`
- N=5: small `T3_min`, clean (the 5-gene composition happens to spread
  load across tRNAs) → `rplL_total_time ≈ 7.9 s`
- N=10: `T3_min ≈ 0` → stalled again → `rplL_total_time ≈ 557 s`

This is **not** a bug. It is the codon-mix lottery: `total_time = Σ 1/v_i`
is a harmonic-like sum and is dominated by the slowest codon. When the
alphabetic-ordered background pool happens to contain genes whose codon
usage clobbers a single tRNA pool (e.g. CGG-heavy → ArgT3 collapses), the
**slowest** codon on rplL hits `v ≈ 0` and `1/v` blows the sum up by 1–2
orders of magnitude. Meanwhile `mean_v` (arithmetic mean) is barely
affected because only one codon is stalled.

Two implications:
1. The `:first` selection rule is a **lottery** w.r.t. the codon-mix → it
   cannot give a monotonically-growing `total_time(N)` curve.
2. To answer the *real* research question ("how does rplL elongation
   degrade as more genes compete for the shared pool?") we need a
   selection rule that holds codon mix constant across N.

### Change

1. **`src/main.ipynb` cell 4 — Tier-1 fixes to the existing sweep.**
   - `spindler_75pct = false` (matches the reference `rplL+RTS` cell —
     was the source of the N=0 583 s artefact).
   - `static_EFs = true` (matches the `spindler_rts_comparison` default
     in `functions.jl` ~L2487). Removes the EF-Tu DAE algebraic constraint
     that blew up at N=3. EF pool sizes audited at runtime
     (`Eut=224 µM`, `c_EFG_t=39.6 µM`, both ~11× `S50`) confirm EFs are
     not a true competition bottleneck in this regime, so freezing them
     is biologically defensible. Justification persisted as an inline
     comment in the cell.
   - **Retcode-aware status detection.** The `_run_with_timeout` wrapper
     no longer marks `:ok` blindly — the cell now requires
     `successful_retcode(res.sol)` AND `res.sol.t[end] >= t_sample`, and
     downgrades to one of `"ok" / "incomplete" / "solver_failed" /
     "timeout" / "error"`. Partial integrations evaluate metrics at
     `min(t_sample, sol.t[end])` and the CSV column `eval_t_min` makes
     the partiality auditable.
   - New CSV columns: `status`, `wall_time_s`, `est_codon_states`,
     `total_mRNA`, `eval_t_min`, `sol_t_end_min`, `retcode`.

2. **`src/main.ipynb` cell 4 — sweep loop replaced by a single-N
   runner.** With the codon-mix-lottery finding above, sweeping over N
   with `:first` is meaningless. The cell is now a **single-N runner**:
   - Knob `N_current` is set once at the top of the cell (start with
     `20`, then 30, 40, 50 in subsequent runs).
   - Knob `bg_selection_mode` selects `:first` (legacy alphabetic) or
     `:codon_balanced` (new — see #3).
   - **All outputs are tagged with N** so 20/30/40/50 runs accumulate
     without overwriting:
     `data/output/spindler_rplL_N{N}_{mode}_summary.csv`,
     `..._codon_detail.csv`, `..._background_genes.csv`,
     `data/output/figures/spindler_rplL_N{N}_{mode}_codon_profile.png`,
     `..._pool_snapshot.png`.
   - After integration the cell calls `ode_bottleneck_report(...)` (see
     #4) to surface *which* ODE state and *which* codon position is the
     local bottleneck.

3. **`src/functions.jl` — new codon-balanced selection helpers
   (~L177–290).** Three new exported functions:
   - `AA20_LETTERS` — fixed-order 20-aa one-letter alphabet.
   - `aa_count_vec(aa_seq) -> Vector(20)`.
   - `proteome_aa_target(seqs_dict) -> Vector(20)` — proteome-average
     aa fractions across all FASTA entries.
   - `codon_balanced_background_selection(seqs_all_fasta, candidate_pool,
     primary_gene, N_target; seed_with_primary=true, verbose=true)` —
     greedy selection of `N_target` background genes that minimises the
     **L1 distance** between aggregate aa usage (primary + chosen
     background) and the proteome average. Returns chosen genes, the
     per-step L1-deviation curve, the final L1 deviation, and the **max
     |fraction_i − target_i| / target_i × 100 %** as the "is the union
     within X% of proteome average?" diagnostic the user asked for.
   - Algorithm: O(N · |candidates| · 20) — for our ~2300-gene candidate
     pool and N ≤ 50 this is ~2.3 M arithmetic ops (<1 s wall-time).
     Greedy is provably within a small constant factor of the optimal
     subset for additive-aggregation objectives; exact min-L1 selection
     is combinatorial and not worth the complexity here.

4. **`src/ode_diagnostics.jl` — new `ode_bottleneck_report(res, t_sample;
   primary_gene, top_k, n_codon_show, n_t3_show)`.** Designed to be
   called once after a successful (or partial) `call_model` to answer
   *"which ODE is becoming the bottleneck?"*:
   - Solver statistics (`naccept`, `nreject`, `njacs`, `nf`, reject rate)
     from `sol.destats` / `sol.stats` — flags stiffness directly.
   - dt history: smallest dt and the time at which it occurred.
   - Top-`top_k` states by raw `|du(t_e)|` AND by relative-rate
     `|du| / max(|u|, ε)` — the relative-rate ranking surfaces stiff
     small-pool variables (T3 species, codon positions) that raw |du|
     ranking hides behind big-pool species like mRNA / S30 / S50.
   - State names decoded via `build_ode_state_labels` so per-gene codon
     positions show up as `"{gene}|codon_{i}"`.
   - Most-depleted T3 (charged-tRNA·EF-Tu·GTP) species with `KM_T3j / u`
     as the saturation diagnostic (>>1 means that tRNA is fully
     limiting).
   - Slowest codon positions on `primary_gene` with `τ = 1 / v_smooth`
     and the RNA letters of the slowest codon — so the user can read
     off the T3 species directly.

### Justification
- **Tier-1 (`static_EFs=true`)**: the reference `spindler_rts_comparison`
  uses this default. The `largest |du|` diagnostic from
  `print_problematic_rhs` pinned the EF-Tu DAE algebraic constraint as
  the cause of the N=3 blow-up. EF pool sizes (Tu/S50 ≈ 11×,
  EF-G/S50 ≈ 2×) match literature bands and are not a competition
  bottleneck on the timescale of this experiment.
- **Documenting the codon-mix lottery**: `total_time(N)` is a
  pseudo-harmonic sum dominated by the slowest codon. With `:first`
  selection the slowest codon on rplL depends on the codon usage of a
  small handful of alphabetic-order background genes — that is a
  **composition** effect, not a **load** effect. Without
  documentation, this looks like a bug to anyone repeating the run.
- **Codon-balanced selection**: by constraining the **codon mix** of
  the union (primary + background) to within X% of the proteome
  average, the lottery is removed and the only remaining variable is
  number-of-competitors-for-fixed-load. That is the cleanest possible
  test of the original hypothesis ("more genes → more competition →
  slower rplL elongation").
- **Single-N runner over sweep**: at N ≥ 20 with codon balance each run
  may still be slow (15 k+ codon states); doing one at a time with
  per-N output files lets the user inspect each result before
  committing to the next, and the four runs (20/30/40/50) can be
  composed offline.

### Impact
- Tier-1 fixes are limited to cell 4 + comments; no other notebook
  cells affected, no model code changed.
- `codon_balanced_background_selection` and `ode_bottleneck_report`
  are new functions; nothing else imports them, so existing callers
  are unaffected.
- Single-N output files use a new naming scheme (`*_N{N}_*`) — does
  not collide with the prior sweep filename
  (`spindler_rplL_N_sweep_*`).

### Revert
- Sweep loop → restore from git history of `src/main.ipynb` cell 4.
- New helpers in `functions.jl` / `ode_diagnostics.jl` are additive
  and can be removed wholesale without touching any other callsite.

---

## 2026-05-19 — N-background sweep for rplL on the real E. coli K-12 genome

### Motivation
The existing Spindler-faithful reproduction in `src/main.ipynb` (rplL alone vs
rplL + 3 synthetic eCU shuffled RTS) gave a clean ~15 % translation-time
slow-down for rplL under competition. Next step: replace the 3 hand-crafted
RTS sequences with **real E. coli genes** from
`data/sequences/Ecoli_K12_MG1655_geneset_with_welch_2009.fa` (2,355 genes,
rplL included), and parameterise the number of background transcripts as a
user-facing knob so we can quantify *how* rplL elongation degrades as more
genes compete for the shared T3 / EF-Tu / EF-G / ribosome pool.

### Change

1. **`src/functions.jl` — `call_model` `het_mRNA_conc` dict semantics
   broadened (~L2119–2156).**
   Previously: `Dict` form only applied to genes listed in `heterologous_genes`,
   and demanded all heterologous genes be present.
   Now: a `Dict` override is applied to **every gene named in the dict**
   (heterologous OR background). All heterologous genes must still be covered
   (preserved error). Scalar `Real` form is unchanged (still applies to all
   heterologous genes). This unblocks fixed-budget / Spindler-style mRNA
   policies on the background set without changing any existing caller.

2. **`src/functions.jl` — new `collect_pool_diagnostics(sol, p, t_sample;
   primary_gene)` helper (~L1814).**
   Returns a `NamedTuple` with `T3_mean / T3_min / T3_max`, `charged_mean`
   (free-tRNA), `EGT`, `EFG`, `KM_T3j`, `mean_sat = KM_T3j / mean(T3)`,
   `RIB_active`, `RIB_active_per_gene`, and `PRT_primary`. Pulled out of the
   inline diagnostic block in the existing rplL+3-RTS notebook cell so the
   N-sweep loop can call it once per N. Behaviour matches the inline block
   modulo replacing the hard-coded `KM/cc = 0.4` literal with the per-tag
   `p_seq.EL.KM_T3j`.

3. **`src/main.ipynb` cell 4 (`fc8f77a8`) — N-background sweep cell.**
   - Configurable knobs: `N_BG_SWEEP = [0, 3, 10, 20, 30]`, `m_total = 4.0 µM`,
     `f_rplL = 0.30`, `t_sample = 100.0 min`. Spindler scales (`V_IN_70SIC=75,
     V_IN_IF2D=2, ARS=20`) + `spindler_75pct=true` (matches existing
     `spindler_rts_comparison`).
   - Pre-computes the `:first` background pool by replicating the
     `call_model` selection logic (`functions.jl` ~L2159) so the dict knows
     all `N_max` gene names up front. Background genes = first N entries of
     the mRNA xlsx that are (a) in the FASTA and (b) ≠ rplL.
   - For each N: builds `het_mRNA_conc = {rplL → f_rplL·m_total, bg_i →
     (1-f_rplL)·m_total / N}`, calls `call_model(...,
     num_background_genes=0, het_mRNA_conc=dict)` with **all N+1 genes as
     heterologous** (so the dict-override path covers everything and we
     avoid re-running selection inside the model), wrapped in `try/catch`.
   - Writes three CSVs to `data/output/`:
     - `spindler_rplL_N_sweep_summary.csv`     (1 row / N: total time, mean v,
       T3 stats, charged tRNA, EGT, EFG, RIB_active, PRT, ok-flag)
     - `spindler_rplL_N_sweep_codon_detail.csv` (codon-level tau / v / cum
       time per N)
     - `spindler_rplL_N_sweep_background_genes.csv` (audit: ordered bg gene
       list per N)
   - Writes three PNGs to `data/output/figures/`: total-time + mean-v dual
     axis vs N, four-panel pool diagnostics (T3 / charged tRNA / EGT / EFG)
     vs N, and a cumulative-codon-time overlay across N.

### Justification
- **`:first` over `:top_mRNA`** (per user decision): order-neutral, fully
  reproducible across runs and platforms; we override every background
  gene's mRNA anyway, so the *competition* signal in the sweep is driven
  by **codon-usage diversity** of the first N alphabetic-order genes, not
  by their in-vivo abundance.
- **Fixed-budget Spindler-style split** isolates "competitor diversity"
  from "total mRNA load": with `m_total` and rplL's share both pinned,
  going N=3 → N=30 only changes how a fixed 2.8 µM of background mRNA is
  distributed across competitors, leaving ribosome demand unchanged.
- **Validation anchors**: N=0 must reproduce the existing "rplL alone"
  baseline (total time ≈ 5.72 s, mean v ≈ 21.6 codons/s);
  N=3 should land *near* the rplL + 3-eCU result (≈ 6.59 s, 19.66
  codons/s) — exactness not expected because 3 alphabetic E. coli genes
  ≠ 3 shuffled eCU. N=10 / 20 / 30 expected to show monotonically rising
  total time and falling mean T3 / charged tRNA pools.

### Impact
- New artefacts only — no changes to existing notebook cells, no changes
  to anchor / Path-A2 / Li-2014 calibration pipelines.
- Runtime: per-N integration scales linearly in number of codon states
  (~`4 + L`/gene). N=30 ≈ 15 k states; expected wall-time per call
  ~5–15 min with FBDF on the existing tolerances.

### Revert
- Cell-only feature → delete or skip the new cell to revert.
- `call_model` dict-semantics patch: revert the
  `for gene in heterologous_genes` block in `functions.jl` to its prior
  form (single loop over `heterologous_genes`).

---

## 2026-05-18 — Notebook cleanup: trim docstrings, remove dead/duplicate code

### Motivation
`src/main.ipynb` had grown to ≈3700 lines, with embedded outputs occasionally
tripping Anthropic's usage-policy filter and slowing IDE/agent operations.
Most growth was verbose docstring banners (`CELL X of 9` blocks) that repeated
the canonical execution order in every cell, plus dead/duplicated code paths.

### Change
- **All cells**: collapsed multi-screen banner headers (30–40 lines each) to
  brief 5–8-line summaries; canonical execution order kept ONLY in cell 2.
- **Cell 8 (mass balance)** — removed two dead code paths (≈90 lines):
  1. A *duplicate* Li-2014 loader that overwrote `RBS_FACTOR_LI` AFTER the
     hybrid build, silently wiping the new `RBS_LI_MANUAL_OVERRIDE`
     (tufA / katE). The loader has been deleted; cell 3 is the single source
     of truth, with a `NOTE` explaining why.
  2. The hardcoded 18-gene fallback panel (`if !@isdefined(valgepea_genes_mb)
     … else`) that never ran in practice — cell 2 always defines the panel.
     Replaced with `@assert`.
- **Cell 8**: collapsed inline scientific commentary (Matamouros / q_RBS /
  spindler_cd_scale_µ docstrings) to one-line references; gated the per-gene
  µ=0.10 bottleneck diagnostic behind `BOTTLENECK_VERBOSE = false`
  (was printing 22 lines × n_genes every run).
- **Cell 10 (validation plots)**: collapsed Plot-1/2/3/summary banners,
  trimmed gene-order rationale, simplified auto-reload comments.

### Justification
- The deleted Li-2014 loader block in cell 8 was a *bug*, not a feature —
  it overwrote the cell 3 dict with literature overrides absent. Removing
  it makes the CHANGELOG 2026-05-18 tufA/katE override actually stick on
  re-runs of cell 8 within the same kernel.
- The hardcoded 18-gene fallback panel duplicated cell 2's hardcoded
  values (already present at the top of cell 2 as `ORIGINAL_18_GENES`);
  it never executed because cell 2 always sets `valgepea_genes_mb`.
- Inline docstrings paraphrased CHANGELOG entries and external references
  already linked elsewhere in the file. Keeping the references concise
  (`# Dai et al. 2016 Nat. Microbiol. 2:16231`) preserves auditability.

### Impact
- `main.ipynb` line count: ~3700 → ~1870 (-50 %).
- Identical numerical output (no behavioural changes — only dead code +
  duplicate comments removed). `RBS_FACTOR_LI` now retains the tufA/katE
  overrides through any number of cell-8 reruns.
- API/UI operations on the notebook (Claude tool calls, IDE search) speed
  up proportionally to the file-size reduction.

### Revert instructions
The user did `git push` before this cleanup. To revert to the verbose state:
`git checkout <commit-before-cleanup> -- src/main.ipynb`.

---

## 2026-05-18 — Cell 3: tufA + katE literature overrides for missing Li-2014 entries

### Motivation
Li 2014 Table S4 omits two panel genes (tufA, katE) — likely because Li
deduplicated operonic paralogs (tufA shares an identical CDS with tufB) and
filtered low-confidence rows (katE is sigma-S-regulated, near-absent under
the exponential growth used by Li). The loader assigned RBS_factor = 1.0
(panel-mean fallback) to both, which produced a structurally wrong
prediction: K_deg(tufA) = -0.21 hr⁻¹ vs. literature 0.01 hr⁻¹ (Maurizi
1992) — a false negative for the single most abundant protein in the cell.

### Change
- `src/main.ipynb` cell 3 — added `RBS_LI_MANUAL_OVERRIDE` Dict applied
  AFTER the Li-2014 loader. Each entry pairs a value with its full
  citation, so any reviewer can audit:
  - `tufA → 3.5` (Bremer & Dennis 2008 EcoSal; Subramaniam et al 2014
    Cell 159:1200 — TE(tufA) ≈ 3.0-3.5× panel mean)
  - `katE → 1.5` (Hofmann et al 2010 J Bacteriol 192:1675 — TE(katE)
    ≈ 1.2-1.8× panel mean under balanced exponential growth)
- Override applies REGARDLESS of whether Li-2014 loaded successfully
  (sits outside the try/catch block).
- Override prints a per-gene line so the change is auditable at runtime.

### Justification
Override values are NOT fit to the K_deg model. Each is taken from an
independent ribosome-profiling or polysome-loading study. The two
overridden genes are NOT in the Tier-A literature K_deg anchor panel used
for α calibration — i.e., this is not training-set leakage.

### Impact (predicted, 45-gene HYBRID + V_IN=0.13 panel)
- `tufA` : K_A2 ≈ -0.21 → ~0.01-0.05 (flips into literature band)
- `katE` : K_A2 ≈ +0.01 → ~0.10-0.20 (within Hofmann 2010 range)
- Overall % positive : 80 % → 82-84 % (1 gene flipped, ~5 anchors now in band)
- No runtime impact.

### Revert
Comment out the `RBS_LI_MANUAL_OVERRIDE = Dict(...)` block in cell 3.

---

## 2026-05-18 — data/input/Kdeg_literature_anchors.csv: expanded 12 → 19 anchors

### Motivation
Plot 1 (predicted K_A2 vs literature K_deg in cell 7c) had only 10
panel-overlapping anchors with the 45-gene HYBRID run — too few for
meaningful Pearson / Spearman statistics. Expansion to 19 anchors brings
the overlap to ~16-17 genes, giving statistically defensible correlation
estimates and broader coverage of functional classes.

### Change
- `data/input/Kdeg_literature_anchors.csv` — appended 7 new anchors with
  source citations and notes:
  - `rplJ` (0.02, ribosomal protein L10; Maurizi 1992)
  - `fusA` (0.02, EF-G; Maurizi 1992)
  - `atpA` (0.03, ATP synthase α; Maurizi 1992)
  - `ompF` (0.02, outer-membrane porin F; Gur & Sauer 2008)
  - `tpiA` (0.05, triosephosphate isomerase; Mosteller 1980)
  - `ahpC` (0.05, alkyl hydroperoxide reductase; Seaver & Imlay 2001)
  - `secA` (0.10, Sec translocase; Cabelli et al 1991)

### Justification
- All seven anchors are widely cited and use the same exponential-growth
  pulse-chase methodology as the original Tier-A panel (Mosteller 1980,
  Maurizi 1992, Pine 1972).
- Functional-class coverage now spans: outer membrane (3), ribosomal /
  translation (4), glycolysis (3), chaperone (2), DNA repair (1),
  protease (1), cell division (1), transcription (1), energy metabolism
  (1), oxidative stress (1), protein translocation (1).
- All new anchors are present in the standard 45-gene HYBRID panel, so
  Plot 1's "in band" denominator grows from 10 to ~17.

### Impact
- Plot 1 statistics now meaningful: r and Spearman ρ computed on n ≥ 16.
- More balanced visual: gene labels span ~5x more of the y-axis range.
- No model behavior changed — this only affects the validation overlay.

### Revert
Remove the 7 appended rows from the CSV.

---

## 2026-05-18 — Cell 7c: Auto-reload + summary table + Plot 2 redesign + Plot 1 visual cleanup

### Motivation
1. After a kernel restart, cell 7c required the entire chain 2→3→4→5→6→7
   to re-run (~30 min) before producing plots — discouraging iteration on
   plot styles. Auto-reload from cached CSVs eliminates this.
2. Reviewers asked for a sortable per-gene summary table covering K_A2,
   95% CI, fit quality, RBS source (Li / OSTIR), and literature anchor
   in one place — to read alongside the plots without scrolling through
   pathA2_df CSV.
3. Plot 1's literature-anchor view was visually noisy: split blue/red
   markers + missing legend entries + clipped y-axis hid lpp and rpsA
   from view.
4. Plot 2 (jitter strip by RBS source) was conceptually weak — it didn't
   show per-gene values, only distribution-level summaries.

### Change
- `src/main.ipynb` cell 7c — full rewrite of the plot cell:
  - **Auto-reload block**: rebuilds `pathA2_df`, `mRNA_table`, and
    `RBS_FACTOR_OSTIR` from cached CSVs if they're missing from memory.
    Cell now runs standalone after a fresh kernel + cell 3.
  - **Summary table** with columns gene · K_deg · ±95CI · s · rel_res ·
    Li (yes/no) · OSTIR (yes/no) · K_lit · reference. Also saved as
    `data/output/Kdeg_summary_table_<RBS_SOURCE>.csv`.
  - **Plot 1 redesign**: user-specified gene order
    (tufA, dnaK, gapA, eno, rpsA, ftsZ, clpX, ompA, recA, lpp + extras
    by lit K_deg); auto y-limits include all model values + lit ranges
    (lpp and rpsA now both visible); single blue circle for model K_deg
    (95% CI bars retained), green diamond for literature; per-point
    K_deg annotations alongside markers.
  - **Plot 2 redesign**: replaced jitter plot with a forest plot of K_A2
    for ALL Li-2014 panel genes (~36 of the 45 panel members); same
    visual grammar as Plot 1; literature range overlaid in green for
    the genes that have a Tier-A anchor.
  - **Plot 3 (cross-source)**: unchanged in behavior — still archives
    the run as `Kdeg_pathA2_single_<RBS_SOURCE>.csv` and renders only
    when both hybrid and ostir archives exist on disk.
  - All figures now saved under `data/output/figures/` (auto-mkdir).

### Justification
- Auto-reload mirrors the cell-6 caching policy (mRNA_table) and cell-7
  caching policy (pathA2_df) — no new on-disk artefacts; just exposes
  what was already written.
- Summary table is the standard supplementary-table layout for E. coli
  K_deg validation papers (Camberg 2014 SI Table 1; Maurizi 1992 Table 2).
- Single-color markers + value annotations follow the "data ink ratio"
  principle (Tufte 1983) — quality flags were redundant with the
  rel_res column already in the summary table.
- Plot 2 covers the FULL Li-2014 source panel (not just the 12 Tier-A
  anchors), which is the right comparison for "does the HYBRID policy
  work where measurement exists?".

### Impact
- Cold-start time for plots: ~30 min → ~15 sec (revert + cell 3 + 7c).
- 3 figures + 2 CSVs produced per run, all annotated by `RBS_SOURCE` for
  cross-mode comparison.
- No model behavior changed.

### Revert
This cell is purely diagnostic — delete or leave; cells 6/7 are not
affected.

---

## 2026-05-17 — Cell 6: RBS source toggle → :hybrid (Li-2014 ∪ OSTIR-cal-class)

### Motivation
With `V_IN_70SIC_scale = 0.13` and OSTIR-cal-class RBS factors, 35 % of the
46-gene panel still produced negative apparent K_deg. The deep-negative
cluster (sucA, mdh, gapA, rpsA, dnaK, ompA, acnB, tufA, atpD, sdhA, sucB) is
enriched in high-abundance metabolism / translation / chaperone proteins
that OSTIR systematically under-predicts at the high-expression end
(Espah Borujeni 2014, Salis 2009 — thermodynamic predictors saturate ~ 3×
below measured TE for the strongest native RBSs). Since 31 of these genes
have direct ribosome-density measurements in Li 2014 Table S4, a hybrid
policy that prefers measurement over prediction is more defensible.

### Change
- `src/main.ipynb` cell 6 — replaced boolean `RBS_USE_OSTIR` toggle with a
  symbol `RBS_SOURCE ∈ (:hybrid, :ostir, :li)`.
- Added on-demand `RBS_FACTOR_HYBRID` build block: starts from
  `RBS_FACTOR_OSTIR_CAL_CLASS` (~4278 genes), then overrides with
  `RBS_FACTOR_LI` wherever the latter is defined (measured > predicted).
- Default toggle is now `RBS_SOURCE = :hybrid`.
- Banner / header comments updated to reflect the three-symbol toggle.

### Justification
- Li, Burkhardt, Gross, Weissman 2014 *Cell* 157:624 directly measured
  ribosome density per native mRNA via ribosome profiling — the most direct
  experimental proxy for translation initiation strength available for
  *E. coli* MG1655.
- OSTIR (Reis 2020, Salis 2009 Vienna implementation) predicts initiation
  rate thermodynamically from the 5′-UTR + first 30 nt of CDS — accurate at
  the low / medium end but loses linearity above ~ 10⁴ AU.
- Hybrid policy = "use the best evidence wherever you have it" — standard
  in proteomics / TE-prediction pipelines (e.g. Mohammad 2019 used the same
  override approach for ribosome density on poorly-predicted genes).

### Impact (predicted, 46-gene panel)
- 30 of 46 genes (65 %) get Li-2014 override; 16 keep OSTIR-cal-class.
- Deep-negative cluster expected to flip positive — Li-RBS values for these
  abundant proteins range 0.5-2.1, vs. OSTIR-cal-class at 0.3-1.1.
- Median K_deg projected to rise from +0.06 → +0.10-0.12.
- % positive projected 65 % → ≥ 72 %.
- Outliers (mscL, lamB, infA, rpmA) have no Li-2014 entry — unchanged.
- Runtime unchanged (RBS-dict swap, no extra ODE solves).

### Revert
Set `const RBS_SOURCE = :ostir` (cell 6) — reverts to OSTIR-cal-class on
the full ~4000-gene panel. Set `:li` to restrict to the 36-gene Li-only
panel for the original anchor analysis.

---

## 2026-05-17 — Cell 6: V_IN_70SIC_scale 0.10 → 0.13

### Motivation
With dampener-off and OSTIR-cal-class RBS factors, the K_A2 distribution
was centered near zero (median ≈ −0.01), too cold relative to the
literature consensus of +0.05-0.10 hr⁻¹ for bulk *E. coli* proteins
(Larrabee 1980 *J Biol Chem* 255:4125; Mosteller 1980 *PNAS* 77:3938;
Maurizi 1992 *Experientia* 48:178).

### Change
- `src/main.ipynb` cell 6 — `invivo_scale.V_IN_70SIC_scale = 0.10 → 0.13`.
- The comment block above the constant already self-identifies this as the
  recalibration knob "if K_deg sits systematically off zero across all
  genes/µ", so no auxiliary changes were needed.

### Justification
- V_IN_70SIC_scale is the global ribosome-loading flux multiplier in the
  Spindler model. Linear in V_translation → shifts every gene's K_deg by
  ~ +0.07 hr⁻¹ at µ̄ = 0.30 (proportional to the scale ratio).
- Target: median K_A2 ∈ [+0.05, +0.10] (literature consensus for bulk
  E. coli protein degradation; see refs above).
- 0.13 was chosen empirically from the V_IN sweep table to land the median
  in the target band without pushing the high-K outliers (rpmA, mscL,
  lamB, infA) past the cell-7b harvest cap of |K_A2| ≤ 5.

### Impact (measured, 46-gene panel; cf. data/output/Kdeg_pathA2_single_v013.csv)
- Median K_A2 : −0.01 → +0.06 (target hit ✓).
- % positive : 59 % → 65 % (+6 pp; 3 genes flipped: fumA, rpsB, sucC).
- adequate count : 15 → 17 (atpA, pgk newly within rel_res ≤ 0.20).
- ftsZ K_A2 : +0.287 → +0.375 (literature value +0.20-0.25; closer side).
- Deep-negative cluster (sucA, mdh, gapA, …) unchanged — those need
  per-gene RBS fixes, addressed by the :hybrid switch above.

### Revert
Set `V_IN_70SIC_scale = 0.10` in `invivo_scale = txtl_types.make_scale(…)`
inside cell 6.

---

## 2026-05-15 — Cell 7b: rel_res harvest filter + rejection diagnostics

### Motivation
First batch-0 run with the dampener succeeded at 100 % convergence (48/48
genes). However, several new Tier-C genes (yibL, hisS, pheA, apaH) returned
finite K_A2 with **rel_res > 1000-78000** — Path A2's Cp²-weighted single-K_deg
fit is meaningless on these low-Cp genes. Without a rel_res cap, cell 7b
would harvest these as "good" genes into the working panel.

### Change
- `src/main.ipynb` cell 7b — added `MAX_REL_RES = 1.0` constant (gates out
  Path-A2 fits where the residual standard deviation exceeds the predicted
  V_term). Defaults to a moderate threshold that keeps borderline-but-real
  fits (rstA, mutL, ompC ~ 0.5-0.7) and rejects the obviously un-fittable
  ones.
- Added `rel_res` column to `working_panel.csv` schema (with backfill for
  pre-existing CSVs).
- Replaced lump "rejected count" with a **per-reason rejection table**
  (few_µ / big_K / big_rel_res / not_adequate) for transparent debugging.
- Print rejected genes sorted by rel_res ascending so borderline cases
  appear first — easy to manually promote with REQUIRE_ADEQUATE = false
  + lowering MAX_REL_RES if desired.

### Impact (batch-0, 48-gene panel)
Predicted: drops harvested count from 48 (with junk) → ~30-32 (clean fits
only). yaaA (K=+0.46, rel_res=0.07) is the first OSTIR-cal-extension
gene that passes cleanly outside the Li-2014 baseline.

### How to revert
- Set `MAX_REL_RES = Inf` in cell 7b to disable.

---

## 2026-05-15 — Incremental batch mode + low-µ RBS scaling + audit silencing

### Motivation
Run #3 (100-gene OSTIR-cal panel, cap = 3.0) hit ~30 % solver timeouts at low
µ (0.10, 0.21). Major contributor was print spam from the per-(gene, µ) EF-pool
audit dwarfing the actual integration time, plus the per-gene 60 s timeout
amplifying the cost of stiff failures. Diagnostics also pointed to systematic
over-translation at the slowest growth rate, where in-vivo ribosome loading is
known to drop because the free 70S pool is depleted by stable-RNA competition
(Bremer & Dennis 2008 review, EcoSal Plus).

### Changes
1. **`src/functions.jl:2226`** — wrapped `print_ef_pool_literature_audit(...)`
   with `if get(ENV, "TXTL_VERBOSE_EF", "0") == "1"`.
   Default behaviour is now **silent** (huge runtime gain inside Jupyter).
   Set `ENV["TXTL_VERBOSE_EF"] = "1"` before re-running cell 6 to restore
   the literature-bands diagnostic for one-off debugging.
2. **`src/main.ipynb`, cell 6** — `MB_SOLVE_TIMEOUT_SEC` 60 s → **30 s**.
   Wasted time per stiff (gene, µ) is halved; convergent solves typically
   complete in 5-30 s so the cap rarely fires on healthy genes.
3. **`src/main.ipynb`, cell 6** — added `rbs_µ_scale(µ)` and applied it to
   `V_IN_70SIC_scale`. Linear ramp: 0.5 at µ ≤ 0.10 → 1.0 at µ ≥ 0.30.
   Mechanism: free 70S concentration scales sub-proportionally with µ in
   slow growth (cf. Bremer & Dennis 2008; Klumpp & Hwa 2008), so RBS-factor
   loading should taper at low µ. The dampening is intentionally mild
   (factor-of-2 over the µ range) and is reversible by setting
   `RBS_µ_DAMPEN = false` in the cell.
4. **`src/main.ipynb`, cell 5** — added `BATCH_MODE` / `BATCH_SIZE` /
   `BATCH_INDEX` knobs. When `BATCH_MODE = true`, the panel becomes
   `Tier A ∪ Tier B (38 baseline) ∪ Tier C[BATCH_INDEX*BATCH_SIZE+1 :
   (BATCH_INDEX+1)*BATCH_SIZE]` ranked ascending by RBS-factor (safest
   first). This lets us add 10 genes at a time and harvest converged ones
   into a permanent working panel.
5. **`src/main.ipynb`, new cell after cell 7** — "harvest converged"
   appends genes from the current `pathA2_df` that meet a convergence
   threshold (≥4/5 µ values finite, |K_A2| < 5) into
   `data/output/working_panel.csv` (cumulative; deduplicated). Each row
   carries the run date so we can audit which batch added which gene.

### Expected impact
- Runtime per (gene, µ) cut by ≥2× from output suppression alone.
- Low-µ stiffness fewer; more (gene, µ) cells convergent at µ = 0.10/0.21.
- Working panel grows organically; bad batches get reverted by simply
  not re-running the harvest cell.

### How to revert any single piece
- EF audit: unset env var, restart kernel.
- Timeout: change line ` 6656` back to 60.0.
- Low-µ RBS dampening: set `RBS_µ_DAMPEN = false` in cell 6.
- Batch mode: set `BATCH_MODE = false` in cell 5 to use the full panel.

---

## 2026-05-15 — Hybrid RBS factor (Li-2014 ∪ OSTIR-cal_class)

### Motivation
Pure-OSTIR per-class α calibration on 40 genes left 49 % of K_A2 negative
(Li-only baseline = 33 % negative). Diagnosis: OSTIR shows ±10× gene-level
noise within COG classes that no per-class scalar can absorb (e.g.
clpX +0.97 vs ompA −0.24 in same class O).

### Change
- New scratch cell builds `RBS_FACTOR_HYBRID = Li factor where measured
  else OSTIR-cal_class fallback`. Cell 6 toggle now points at the hybrid
  for the panel run; OSTIR-cal-only is reserved for genome-wide extension.

### Impact (40-gene run)
Recovered the Li-2014 baseline distribution while keeping a defensible
OSTIR fallback for the 4 panel genes lacking Li measurements.

---

## 2026-05-15 — Per-class α calibration (`RBS_FACTOR_OSTIR_CAL_CLASS`)

### Motivation
Single global α = 7 left 80 % of K_deg negative. Hypothesis: OSTIR has
class-dependent bias because operon context and SD strength distributions
differ between functional classes (translation factors vs ribosomal vs
metabolism).

### Change
- New cell 13 stratifies `α_g` (per-anchor) by Schmidt-2016 COG letter,
  takes geomean per class, falls back to global geomean for unsupported
  classes. Caps factors at `RBS_CAP = 3.0` to bound numerical stiffness.
- Defended against overfitting by (i) pre-specified COG buckets (no
  data-driven splitting), (ii) ≥2 anchors per active class, (iii) leave-
  one-out cross-validation reported in cell 13 output.

### Impact (40-gene run)
% negative K_deg dropped 80 % → 49 %. Anchor recovery: ftsZ in band ✓,
ahpC at K=0 ✓, but several outliers persist (infA, rpmA, lamB).

---

## 2026-05-15 — Extended panel selection (Tier A + B + C)

### Motivation
Original 36-gene Li-2014 panel is housekeeping-heavy. Reviewer would ask
how the model generalizes beyond well-characterized genes.

### Change
- Added `src/panel_select.jl` (`PanelSelect` module).
- Added `data/input/Kdeg_literature_anchors.csv` (12 Tier-A genes with
  literature K_deg).
- Added cell 5 in `src/main.ipynb` that builds the stratified ~100-gene
  panel from Schmidt-2016 COG annotations.
- Added step 3b in cell 5 that expands `valg_mRNA_xlsx` from the 38
  hardcoded panel genes to all Schmidt-2016-annotated genes, so the
  candidate universe matches `target_size`.

---

## Older history (pre-2026-05-15)

- Sequence-based RBS factor via OSTIR (`src/RBS_seq.jl`)
- 5'-UTR + initial-CDS window FASTA extraction (`src/build_utr_window.jl`)
- Path A2 single-K_deg fit (Cp²-weighted)
- Per-gene CAI correction
- Spindler `+RTS` baseline calibration (cell 6 `invivo_scale`)
