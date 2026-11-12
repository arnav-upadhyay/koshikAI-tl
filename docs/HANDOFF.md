# Handoff — TXTLv2 Subramaniam validation (read this in a new chat)

**Repo:** `e:\TxTl_v2_main\TXTLv2`  
**Date of this file:** 2026-11-11  
**User:** one notebook solve at a time; do not start the 6-case script in parallel with the notebook.

Paste into a new Agent chat: *Read `docs/HANDOFF.md` and `docs/validation_ladder.md`. Independent-MM 1B failed. Competitive Leu-RS is in `rhs_metabolites.jl`. Re-run packages, then 1B into the competitive CSV. Do not turn freeze on. Do not edit the stall-freeze plan. Commits only if asked.*

---

## Goal

Mechanistic translation model (Julia ODE). Validate **translation flux** with **promoter and RBS locked**, then later add promoter, then RBS. Subramaniam 2013 is step 1: same initiation, YFP CDS CTG vs C1A vs CTA, Leu-ARS = charging proxy.

PRT = integral of termination (`du[PRT]=V_TL_TE`), **no µ dilution**. Compare **folds**, not µM vs copies/cell. **No protein multiplier.**

---

## Do not

- `freeze_stalled = true` (shared tRNA is not fixed by gene freeze)
- Delete repo files or edit a stall-freeze plan file
- Equal mRNA on background genes or production `f_rplL = 0.46` as a Subramaniam cell
- `include("packages.jl")` inside the 1B cell (`replacing module txtl_types`)
- Point `load_csv` at a FASTA (`column name :gene not found`)
- Paracoccus RBS models on this *E. coli* core
- Treat Kudla/Goodman N-termini as elongation-only
- `q_RBS` is TASEP occupancy, **not** Shine–Dalgarno
- Rebuild the 21-gene (10+10) FASTA — it hung >8 h

---

## Already done (do not redo)

### Model

- Track A: `t3_efg_saturation` — MM rate 0 at empty T3. N=5 freeze-off: PRT_rplL≈2.03, wall ~455 s.
- Track B: delay vs protein; freeze does not refill shared T3.
- Production notebook cell: N=20 equal `m_bg`, `f_rplL=0.46`, V_IN=75, ARS=20 — **different experiment**. Leave it.

### 1A — YFP alone (`n_seq=1`) — done

`data/output/pnas_yfp_cds_only.csv`. V_IN=1, mRNA=1 µM.

- Rich matches paper: PRT 1.746 / 1.742 / 1.751.
- Starved fails: CTA/CTG=0.95; all drop ~5×. Reporter *is* the proteome.
- C1A at 0.1 **invalid** (RIB=140, CTApool=99).

### Background (current files)

`n_head=5`, `n_tail=5`, `f_reporter=0.05`, `m_total=2.86`.  
Genes: **rplL**, lpp, rplE, rplU, rpmJ, tufA | alaS, oppB, ppc, tpx, acrB.  
Proteome CTG/Leu≈0.62; reduced≈0.65. Tail µM **can** exceed head (bulk stand-ins).

Files: `data/output/reduced_cell_load.csv`, `data/sequences/reduced_yfp/bg_{CTG,C1A,CTA}.fa`, `data/sequences/reduced_cell_h5_t5.fa`.

### 1B six-case — **done; fail at Leu-ARS = 0.1** (`n_seq=11`)

Table: `data/output/reduced_yfp_leu_ars.csv` (keep these six rows).  
Notebook: last code cell in `src/main.ipynb`. Packages cell first. Freeze off. `V_IN=1`, `t_end=100`.

| lab | Leu-ARS | PRT | vs CTG | starved/rich | RIB | ctaτ | CTApool | wall |
|-----|---------|-----|--------|--------------|-----|------|---------|------|
| CTG | 1 | 0.006585 | 1.00 | — | 20.14 | 0 | 0.584 | 5038 s |
| C1A | 1 | 0.005096 | 0.77 | — | 20.10 | 0.058 | 0.658 | 3276 s |
| CTA | 1 | 0.004275 | 0.65 | — | 20.11 | 0.684 | 0.668 | 3466 s |
| CTG | 0.1 | 0.009397 | 1.00 | **1.43** | 20.10 | 0 | 1.10 | 2313 s |
| C1A | 0.1 | 0.009020 | 0.96 | **1.77** | 20.09 | 0.056 | 1.08 | 4260 s |
| CTA | 0.1 | 0.007270 | **0.77** | **1.70** | **48.7** | 0.631 | 1.09 | 5623 s |

**Paper test failed.** Starved CTA/CTG = **0.77**, not ≪ 1, and **worse** than rich 0.65 (gap closed). All three PRT **increased** when Leu-ARS dropped (reallocation onto 5% reporter, not a charging collapse). CTA pool ~1 µM, not emptied. C1A did not become an intermediate stall (0.77 → 0.96 vs CTG).

CTA 0.1 RIB=48.7 vs S50≈20 and T3min=−0.004: strained; not 1A-invalid. CTG/C1A at 0.1 were healthy (RIB≈20).

Rich: mild CTG/CTA split 1.54 vs paper &lt;1.4. PRT ≪ 1.75 because reporter is 5% of the cell — expected; **do not** scale PRT.

**1B cell last independent-MM knobs:** `lab="CTG"`, `leu_ars=0.01`. CTG PRT at 0.01 = **0.009396** (same as 0.1). Stopped the 0.01 trio.

---

## NEXT — competitive Leu-RS (in the model now)

Independent MM failed because each of t1l…t5l drew a **full** `V_ARS[10]`. Charging is now

`v_i = V · (t_i/KM_i) / (D0 + Σ_j t_j/KM_j)` so `Σ v_i ≤ V`.

Other amino acids unchanged. Freeze still off. Old table `reduced_yfp_leu_ars.csv` is the independent-MM record — **do not append to it**.

**CTG/1 competitive vs independent is bit-identical** (PRT 0.006585, pools 0.5836 / 12.1245). That is expected if unloaded Leu tRNA stays tiny: both rate laws collapse to `v ∝ t_unloaded`. Competition only splits V when unloaded tRNA accumulates (starve). Do not re-run CTG/1. Next: C1A/1 then CTA/1 then 0.1. 1B cell now prints `LEU_RS_COMPETITIVE=true` — that line must appear. The starve trio is the real test.

**Run:** restart/re-run the **packages** cell (rhs changed). Then 1B cell:

```julia
lab     = "CTG"
leu_ars = 1.0
```

`out_csv` = `../data/output/reduced_yfp_leu_ars_competitive.csv`

Order: CTG/1 → C1A/1 → CTA/1 → CTG/0.1 → C1A/0.1 → CTA/0.1.

Rich may drop vs independent 0.0066 because total Leu charging is ≤1×V not 5×V. Still compare **folds**. Pass: starved CTA/CTG ≪ 1 and ≪ that run’s rich fold; RIB not exploding.

Do not freeze. Do not multiply PRT. Do not rebuild 21 genes.

---

## How to run 1B

1. Re-run the packages cell (`include("packages.jl")`) so `rhs_metabolites.jl` reloads.  
2. Change **only** `lab` and `leu_ars`.  
3. `load_csv = "../data/output/reduced_cell_load.csv"` — never a `.fa`.  
4. Write to `reduced_yfp_leu_ars_competitive.csv`. Leave `reduced_yfp_leu_ars.csv` as the independent-MM archive.

CLI: from `src/`  
`julia --project=.. --threads=12 run_reduced_leu_ars.jl`  
User prefers **one notebook solve at a time**.

---

## Key code

| Path | Role |
|------|------|
| `src/packages.jl` | loads model + `reduced_cell_load.jl` |
| `src/functions.jl` | `call_model`, `t3_efg_saturation` |
| `src/rhs_metabolites.jl` | competitive Leu-RS (`V_leu1`…`V_leu5`) |
| `src/rhs_sequence.jl` | `du[PRT]=V_TL_TE` |
| `src/reduced_cell_load.jl` | mini-transcriptome + YFP FASTA writer |
| `src/make_reduced_background.jl` | builder knobs |
| `src/run_reduced_leu_ars.jl` | 6-case CLI (hardcoded 1.0 and 0.1 — edit if 0.01) |
| `src/data/PNAS_2013/sd02.txt` | yfp0 / yfp1 CDS |
| `docs/validation_ladder.md` | paper ladder |

Sites: `[42,44,46,53,60,64,68]`; C1A = 42. Reporter FASTA name: `rplL`. Leu-ARS index: `AA20_LETTERS` L → 10.

---

## Later ladder (not started)

2A Subramaniam 2014 Cell ribo-seq → 2B Goodman 2013 **one** promoter/RBS slice → 3 Kosuri promoter slice → 4 full P×RBS → 5 Balakrishnan/Schmidt via `v/µ`.

---

## Solver notes

- 11-gene 5+5 converges (~40–90 min). 21-gene hung >8 h.  
- `timeout_s` in the CLI script only **warns**.
