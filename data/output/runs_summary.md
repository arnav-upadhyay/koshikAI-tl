# Successful run index — `spindler_rplL_elongation_log.csv`

This file is the **canonical row identifier** for the elongation log. The log has accumulated
rows from many tuning attempts; this document lists which `(timestamp, N)` rows are the
**physically valid** ones and which to ignore.

> Each row in the elongation log corresponds to **one codon position** of rplL.
> A single simulation contributes **120 rows** (the length of the rplL ORF in codons).
> Disambiguate runs by the `(timestamp, N)` pair — same `(timestamp, N)` ⇒ same simulation.

## TL;DR — canonical rows for cross-N analysis

| N  | timestamp                  | Σ τ (s)   | mean v (cod/s) | T3 mean (µM) | atol  | wall (s) | scale regime |
|---:|----------------------------|----------:|---------------:|-------------:|------:|---------:|:--|
| 10 | `2026-05-27T10:12:49.800`  | **7.749** | 16.84          | 6.29         | 1e-4  | 1265     | run-I (low-N rerun) |
| 20 | `2026-05-27T15:37:32.860`  | **7.538** | 17.44          | 6.33         | 1e-5  | 1453     | run-I (low-N rerun) |
| 25 | `2026-05-22T23:19:21.121`  | **7.450** | 17.52          | 5.64         | 1e-6  | 3376     | run-F era (V_IN=25, tRNA=1.5) |
| 30 | `2026-05-25T15:19:51.616`  | **7.523** | 17.45          | 6.24         | 1e-5  | 2035     | run-I′ (post-warm-start) |
| 35 | `2026-05-23T16:43:56.419`  | **7.727** | 17.21          | 6.40         | 1e-6  | 1900     | run-K (V_IN=17, tRNA=1.8) |
| 40 | `2026-05-24T11:35:35.939`  | **7.732** | 17.19          | 5.66         | 1e-5  | 4228     | run-M (V_IN=15, tRNA=1.6, ARS=20) |
| 45 | `2026-05-25T22:19:22.412`  | **7.515** | 17.38          | 6.20         | 1e-5  | 6971     | run-O |
| 50 | `2026-05-26T14:32:38.477`  | **7.450** | 17.53          | 6.26         | 1e-4  | 8659     | run-Q |

**8 datapoints, N ∈ [10, 50]. Σ τ ∈ [7.45, 7.75] s ⇒ spread = 0.30 s ≈ 4.0% of mean.**

The 4% spread is dominated by the three runs at non-run-I scales (N=25, 35, 40).
Restricted to **strict run-I scale set** (`V_IN=20, ARS=20, tRNA=1.7, EF=√2≈1.41, Rib=2`),
we have **N ∈ {10, 20, 30, 45, 50} → Σ τ ∈ [7.45, 7.75] s, spread 4.0%** (still dominated by
N=10 at 7.75; the other four cluster tightly at 7.45–7.54 with spread of only 0.09 s = 1.2%).

## Universal "run-I" scale set (validated for N ≥ 20)

```
Rib_scale        = 2.0
EF_scale         = sqrt(Rib_scale) ≈ 1.4142    # works around EFG-double-scaling bug, see scale_tuning_log.md §3.4
tRNA_scale       = 1.7
V_IN_70SIC_scale = 20.0
V_IN_IF2D_scale  = 2.0
ARS_scale        = 20.0  .* ones(20)            # all 20 amino-acid ARSs
ode_reltol       = 1e-8
ode_abstol       = U-shaped by N (see below)
```

**atol "U-shaped" ladder** (post-warm-start code only — `scale_tuning_log.md` §5.2):

| N range  | atol | rationale |
|:--|:--|:--|
| 1 – 10   | **1e-4** | Severe Ala-T3 transient at low N (high-mRNA Ala-rich genes dominate the budget); microsec dt at tighter atol |
| 15 – 45  | **1e-5** | Balanced mRNA distribution, no transient spikes |
| 50+      | **1e-4** | Minority isoacceptor depletion (Leu/Val/Ser/Thr) grinds tight-atol solver |

## How each row was generated (canonical-row provenance)

### N = 10  (atol = 1e-4)
- **First 10 genes**: `lpp, rplE, rplU, rpmJ, tufA, fusA, ompF, tufB, acpP, rpsK`
- **mRNA**: `m_rplL = 1.2 µM, m_bg_each = 0.28 µM, total = 4.0 µM`
- **Result**: Σ τ = 7.75 s, mean v = 16.84 cod/s, slowest codon pos 68 (GUU). RIB_active = 54 µM, PRT_rplL = 2.82 µM.
- **Note**: Σ τ is the highest of the cluster (~3% above the N≥20 mean). This is real biology: with only 10 background genes, the rplL gets a larger share of ribosomes, slightly slower per-codon dwell from local tRNA depletion at high-flux codons. Not a numerical artifact (RIB_active and T3_mean are healthy).

### N = 20  (atol = 1e-5)
- **First 20 genes** (10 above + ): `gapA, rpsL, rpsM, rplX, rplO, cspC, rpsS, ompA, rplN, rplK`
- **mRNA**: `m_bg_each = 0.14 µM, total = 4.0 µM`
- **Result**: Σ τ = 7.54 s, mean v = 17.44, slowest codon pos 68 (GUU). RIB_active = 55.6 µM, PRT_rplL = 1.07 µM. T3_min = -5.7e-5 (negligible negative-overshoot artifact, ubiquitous across all clean runs).
- **Note**: An earlier N=20 attempt at atol=1e-4 (ts `2026-05-27T10:37:47`) drifted to RIB_active = 406 µM and T3_mean = 15.5 µM — physically impossible. **Do NOT use that row.** The 1e-5 retry is the canonical row.

### N = 25  (atol = 1e-6 — tuning-era F)
- **scales**: `V_IN=25, ARS=20, tRNA=1.5, EF=1, Rib=1` — these are **not** run-I scales.
- **mRNA**: `m_bg_each = 0.112 µM, total = 4.0 µM`
- **Result**: Σ τ = 7.45 s, mean v = 17.52. The result is healthy because at this gene count the model is more forgiving — `scale_tuning_log.md` §4.2 Run F documents this.
- **Status**: kept in the cross-N plot for continuity, but flagged: different scale regime. A rerun at strict run-I scales would be valuable but is not strictly necessary because Σ τ already sits in the cluster.

### N = 30  (atol = 1e-5) — Run I′
- **First 30 genes**: continues the mRNA-xlsx order through `rpoA, rplM, rpsJ, rpsH, icdA, rplV, rpmB, rplJ, b1824, rplW` (last 10).
- **mRNA**: `m_bg_each = 0.0933 µM`
- **Result**: Σ τ = 7.52 s, mean v = 17.45. Anchors the cluster.
- **Note**: This row is the **post-warm-start** rerun. An earlier N=30 row (`2026-05-23T00:48:04.722`, Σ τ = 7.55) is the **pre-warm-start** Run I and is biologically nearly identical (Σ τ difference: 0.025 s = 0.3%). Keeping only the I′ row.

### N = 35  (atol = 1e-6 — Run K, non-run-I scales)
- **scales**: `V_IN=17, ARS=20, tRNA=1.8, EF=?, Rib=?` — tuning-era K.
- **Result**: Σ τ = 7.73 s, mean v = 17.21. Slightly elevated like N=10. Plausibly the same low-N-like behavior re-emerging because V_IN was lowered (less initiation → fewer rplL ribosomes per cycle → slightly slower per-codon dwell).
- **Status**: kept for continuity; would be cleaner to rerun at strict run-I (V_IN=20, tRNA=1.7) but Σ τ sits within the cluster band.

### N = 40  (atol = 1e-5 — Run M, non-run-I scales)
- **scales**: `V_IN=15, ARS=20, tRNA=1.6, EF=?, Rib=?`.
- **Result**: Σ τ = 7.73 s, mean v = 17.19. Same comment as N=35.

### N = 45  (atol = 1e-5) — Run O
- **First 45 genes** — continues through `rpsD, gatZ, infC, prlA, pyrB, rplA, rplF, rpsG, rplD, yceD, rpsA, rplC, nlpD, rpmG, gatY`.
- **mRNA**: `m_bg_each = 0.0622 µM`
- **Result**: Σ τ = 7.51 s, mean v = 17.38. Right in the cluster.

### N = 50  (atol = 1e-4) — Run Q
- **First 50 genes** — adds `gatA, rpsC, cysK, rplY, ?`.
- **mRNA**: `m_bg_each = 0.056 µM`
- **Result**: Σ τ = 7.45 s, mean v = 17.53. RIB_active = 53.4 µM. Cleanest run of the set.
- **Note**: atol loosened to 1e-4 to handle minority isoacceptor depletion (Leu/Val/Ser/Thr fell to 10⁻⁴ µM at atol=1e-5 mid-run — see `scale_tuning_log.md` §4.2 Run Q-att-1). Wall time: 8659 s = 2.4 h.

---

## Rows to IGNORE in the elongation log

These rows are in the CSV but represent failed / non-converged / numerically-broken
simulations. Filter them out by **timestamp** in any cross-N analysis.

| N  | timestamp                  | Σ τ (s)   | reason for exclusion |
|---:|----------------------------|----------:|:--|
| 0  | `2026-05-26T21:42:36.526`  | 2.0e6     | rplL elongation frozen (v ≈ 6e-5 cod/s). No background genes ⇒ pool dynamics get stuck. Not a usable N=0 datapoint. |
| 1  | `2026-05-21T11:59:28.097`  | 11.55     | OLD scales (V_IN=75, tRNA=1, ARS=20, EF=1, Rib=1). Pre-tuning baseline. |
| 1  | `2026-05-26T21:46:19.071`  | 2449.6    | Stalled mid-run with run-I scales (mean v = 10.4 cod/s ≪ expected 17). N=1 may have similar issues to N=0 — too few background genes to stabilize pools. |
| 10 | `2026-05-21T11:38:42.467`  | 550.0     | OLD scales failure (V_IN=75, tRNA=1). |
| 15 | `2026-05-20T23:41:49.137`  | 1903.2    | OLD scales failure. |
| 20 | `2026-05-21T01:04:49.071`  | 677.5     | OLD scales failure. |
| 20 | `2026-05-27T10:37:47.492`  | 6.99      | atol=1e-4 drift: T3_mean = 15.5 µM (2.5× normal), RIB_active = 406 µM (impossibly above S50 pool). Σ τ looks plausible but underlying state is unphysical. |
| 30 | `2026-05-23T00:48:04.722`  | 7.55      | Pre-warm-start Run I. Functionally identical to Run I′ but the warm-start changed the t=0 code path; the I′ row is the canonical post-fix one. (You can use either; just not both.) |

### Status of N = 0 and N = 1

Neither N=0 nor N=1 has produced a usable datapoint yet. Hypothesis: at very low N the
**ribosome / tRNA / 70SIC steady-state landscape qualitatively changes** — there isn't enough
background protein synthesis to "drain" the pools to the working level seen at N ≥ 10. The
shared metabolite pools end up in a different basin and rplL elongation effectively halts
(N=0) or grinds (N=1).

If you want a usable low-N reference, **N = 10 is the practical floor** with the current
model and scale set. Trying to push lower would require a separate tuning study and is
probably not scientifically informative.

---

## Quick lookup: which raw artifact files belong to which canonical row?

The per-N tagged files (`spindler_rplL_N{N}_first_*.csv`) are **overwritten** each time you
run that N, so they always reflect the **latest** run for that N. As of 2026-05-27:

| N | latest summary file (= canonical row) | matches log timestamp |
|---:|:--|:--|
| 10 | `spindler_rplL_N10_first_summary.csv` | **STALE** — overwritten by failed atol=1e-5 retry (status=timeout, all pool fields NaN). The canonical N=10 row in the elongation log still has codon-level τ and T3 data, but if you want a fresh pool snapshot (RIB_active, PRT, EGT) at N=10, rerun N=10 at atol=1e-4 to regenerate this file. |
| 20 | `spindler_rplL_N20_first_summary.csv` | 2026-05-27T15:37:32.860 ✓ |
| 25 | `spindler_rplL_N25_first_summary.csv` | 2026-05-22T23:19:21.121 ✓ (never re-run since) |
| 30 | `spindler_rplL_N30_first_summary.csv` | 2026-05-25T15:19:51.616 ✓ |
| 35 | `spindler_rplL_N35_first_summary.csv` | 2026-05-23T16:43:56.419 ✓ |
| 40 | `spindler_rplL_N40_first_summary.csv` | 2026-05-24T11:35:35.939 ✓ |
| 45 | `spindler_rplL_N45_first_summary.csv` | 2026-05-25T22:19:22.412 ✓ |
| 50 | `spindler_rplL_N50_first_summary.csv` | 2026-05-26T14:32:38.477 ✓ |

The summary CSVs have additional fields not in the elongation log: `RIB_active`,
`PRT_primary`, `EGT`, `EFG`, `aa_balance_L1`, etc. — useful for the pool-stability plots.

---

## How `plot_cross_N.py` uses this index

The plotting script (`src/plot_cross_N.py`) hard-codes the canonical-timestamp dict from
the table above. Whenever you add a new canonical row (e.g., a fresh N=20 rerun at strict
run-I scales), update **both** this MD file and the `CANONICAL_TIMESTAMPS` dict at the
top of `plot_cross_N.py`. The script silently filters the log to those rows only.
