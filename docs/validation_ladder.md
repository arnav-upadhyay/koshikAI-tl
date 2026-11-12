# Validation ladder — translation first, then promoter, then RBS

Goal: show that **translation flux** (elongation + shared tRNA/ribosome pools) can explain protein differences **when promoter and RBS are the same**. Only after that do we add promoter, then RBS. Native proteomes (Schmidt, Balakrishnan) mix all three and come last.

PRT in this model is chains made by harvest (`du[PRT] = V_TL_TE`, no µ dilution). Compare **folds** and **rates** to Subramaniam. Do not expect raw µM = copies/cell.

Freeze stays **off**. Production equal-`m_bg` / `f_rplL = 0.46` cell is a different experiment and is not this ladder.

---

## Step 1 — Subramaniam 2013 (promoter and RBS fixed)

**Paper:** Subramaniam, Pan, Cluzel. *PNAS* (2013) 110:2419. Same pZS*11 promoter and RBS, YFP CDS synonymous Leu variants, rich vs Leu limitation. Paper: rich PRT nearly equal; starved CTA ≪ CTG (up to ~100-fold). Mechanism: isoacceptor competition for charging, not codon usage tables.

**Model map:** same `V_IN`, same `het_mRNA_conc`; only YFP CDS (CTG / C1A / CTA); Leu-ARS scale = charging proxy. Sites: `[42, 44, 46, 53, 60, 64, 68]`; C1A = site 42 only.

### 1A. Reporter alone (`n_seq = 1`) — done

Sweep: `src/main.ipynb` YFP cell. Table: `data/output/pnas_yfp_cds_only.csv`. Freeze off, `V_IN = 1`, mRNA = 1 µM, t = 100 min.

| Construct | Leu-ARS | PRT (µM) | vs CTG | starved / rich | ctaτ | Notes |
|-----------|---------|----------|--------|----------------|------|--------|
| CTG | 1 | 1.746 | 1.000 | — | 0 | |
| C1A | 1 | 1.742 | 0.998 | — | 0.073 | |
| CTA | 1 | 1.751 | 1.003 | — | 1.23 | CTA pool empty even when rich |
| CTG | 0.1 | 0.331 | 1.000 | 0.19 | 0 | RIB 40 |
| C1A | 0.1 | 0.504 | 1.52 | 0.29 | 0.11 | **invalid:** RIB 140, CTApool 99 µM |
| CTA | 0.1 | 0.315 | 0.95 | 0.18 | 1.56 | CTA pool empty |

**What holds:** rich-medium degeneracy. Same initiation, synonymous CDS → PRT within ~0.3%. That is the paper’s grey bars (their spread was < 1.4-fold).

**What does not:** starved ranking. CTA/CTG at Leu-ARS = 0.1 is **0.95**, not ≪ 1. All three drop ~5×. CTA dwell and empty CTA pool are real; harvest protein does not split. Cause: the reporter **is** the proteome, so cutting Leu-ARS hits every Leu codon on that one mRNA. In the paper the reporter is a small extra on a CTG-heavy chromosome.

**Do not claim:** C1A “wins” under starvation (numerical blow-up). Do not claim freeze. Do not map 1.75 µM to copies/cell.

### 1B. Same reporters + reduced transcriptome (5 head + 5 tail, `n_seq = 11`) — **done; fail at Leu-ARS = 0.1**

Notebook 1B cell; table `data/output/reduced_yfp_leu_ars.csv`. Reporter 0.143 µM (5% of 2.86). Background: lpp, rplE, rplU, rpmJ, tufA + alaS, oppB, ppc, tpx, acrB. Freeze off.

| Construct | Leu-ARS | PRT (µM) | vs CTG | starved/rich | RIB | ctaτ | CTApool | wall |
|-----------|---------|----------|--------|--------------|-----|------|---------|------|
| CTG | 1 | 0.006585 | 1.000 | — | 20.14 | 0 | 0.584 | 5038 s |
| C1A | 1 | 0.005096 | 0.774 | — | 20.10 | 0.058 | 0.658 | 3276 s |
| CTA | 1 | 0.004275 | 0.649 | — | 20.11 | 0.684 | 0.668 | 3466 s |
| CTG | 0.1 | 0.009397 | 1.000 | **1.43** | 20.10 | 0 | 1.10 | 2313 s |
| C1A | 0.1 | 0.009020 | 0.960 | **1.77** | 20.09 | 0.056 | 1.08 | 4260 s |
| CTA | 0.1 | 0.007270 | **0.774** | **1.70** | **48.71** | 0.631 | 1.09 | 5623 s |

**Verdict:** adding a CTG-heavy background did **not** produce Subramaniam’s starved ranking. CTA/CTG at 0.1 is **0.77**, *closer* to 1 than rich 0.65 (paper: ≪ 1, up to ~100×). All three harvests **rose** when Leu-ARS was cut (reallocation onto the 5% reporter). CTA pool stayed ~1 µM (not emptied). C1A’s one rare codon did nothing under starve (0.77 → 0.96 vs CTG).

CTA 0.1: RIB = 48.7 (S50 ≈ 20) and T3min = −0.004 — strained, not 1A-invalid (RIB 140 / CTApool 99). Do not treat it as a 100-fold stall.

**Do not:** protein multiplier, freeze, 21-gene FASTA, production `f_rplL=0.46` cell.

**If continuing 1B:** independent MM failed even at Leu-ARS = 0.01 (CTG PRT still 0.00940). Next is **competitive Leu-RS** in `rhs_metabolites.jl` (five isoacceptors share one V_max). Re-run packages, then 1B into `reduced_yfp_leu_ars_competitive.csv` (keep the old CSV). Pass remains: starved CTA/CTG ≪ 1 and ≪ rich 0.65, RIB not exploding.

Full session state: `docs/HANDOFF.md`.

---

## Step 2 — still translation, still fixed P+RBS

| After 1B | Paper | Why |
|----------|--------|-----|
| 2A | Subramaniam 2014 *Cell* ribo-seq (Leu/Ser starve) | Occupancy at CTA vs CTG on many genes; mRNA mix should not be the story |
| 2B | Goodman et al. 2013 *Science* **slice**: one promoter, one RBS, N-terminal codon variants | Initiation-region codon/structure with P+RBS locked. If protein moves and 5′ fold is the cause, that is **not** elongator starvation (your `q_RBS` is occupancy, not SD) |

---

## Step 3 — add promoter (RBS + CDS fixed)

Kosuri et al. 2013 *PNAS* 1301301110: **slice** one RBS, same sfGFP, vary promoter. mRNA should follow promoter; protein/mRNA stay similar. In the model: `V_TC` / mRNA scale, not a new CDS. Pinning mRNA from a table skips this test.

---

## Step 4 — promoter + RBS + CDS

Kosuri 2013 full promoter × RBS grid (same GFP). Goodman 2013 full promoter × RBS × N-terminus. Mutalik et al. 2013 *Nat Methods* (BIOFAB). Needs an initiation module that is not TASEP `q_RBS`. Do not use Paracoccus RBS models on this *E. coli* core.

---

## Step 5 — native proteome (all knobs mixed)

Balakrishnan 2022 *Science*; Schmidt 2016. Compare \(v/\mu\) (time in min ⇒ divide by µ/60 if µ is h⁻¹), then volume if copies/cell. Last, not first.

---

## What we will not do on this ladder

- Protein multiplier on PRT.
- Equal mRNA on 20 genes or `f_rplL = 0.46` as a Subramaniam cell.
- Freeze as a substitute for drop-off or for moving past a shared empty tRNA.
- Calling Kudla 2009 or Goodman N-termini “elongation only” without checking 5′ structure.
