# TXTLv2 — Cell-Free / In Vivo Translation Model

**Branch:** `arnav-tripwire-test4feb26`

Written by: Arnav Kumar Upadhyay, Cedric Meyer — adapted from MATLAB code of Jan Spindler
([DOI: 10.1021/acssynbio.4c00578](https://doi.org/10.1021/acssynbio.4c00578))

---

## Summary of Changes (relative to last commit on this branch)

### Overview

| File | Lines changed | Category |
|------|--------------|----------|
| `src/rhs_metabolites.jl` | +76 | Nonnegativity guards, DAE fix, PIC fix |
| `src/rhs_sequence.jl` | +34 | q_R_elo clamping, mRNA degradation, C_vect guard |
| `src/functions.jl` | +376 | call_model rewrite, mRNA degradation, diagnostics, IC scaling |
| `src/datatypes.jl` | +26 | `make_scale()`, `delta` field |
| `src/parameters.jl` | +2 | EF_scale documentation |
| `src/packages.jl` | +2 | `ode_diagnostics.jl`, `successful_retcode` |
| `src/ode_diagnostics.jl` | new file | Diagnostic print functions |
| `src/main.ipynb` | +6000 | Taniguchi & Valgepea validation, sensitivity sweep, CAI correction |

---

## Detailed Change Log

### 1. mRNA Degradation Model

**What:** Added first-order mRNA degradation following Calabrese et al. *Nature Physics* (2024).
Previously mRNA was a constant (`tags.RNA`); now `rna` is a dynamic state variable with
`d[rna]/dt = production - delta * rna`, where `delta = ln(2) / t_half` (default t_half = 3 min).

**Why:** To allow mRNA levels to respond dynamically to translation load and enable future
growth-rate-dependent mRNA decay.

**Where:**
- `src/functions.jl`: Added `degradation_rate()`, `degradation_flux()`, `effective_mrna()` functions (~lines 1770-1830). `generate_local_tags()` now computes `delta` per gene (~line 1121).
- `src/datatypes.jl`: Added `delta::Float64` field to `data_struct_integration` (~line 496).
- `src/rhs_metabolites.jl`: Changed `du[tags.idx.rna]` from `V_IN_70SIC` to `0.0` (reset; degradation added in rhs_sequence) (~line 168).
- `src/rhs_sequence.jl`: Added `production_rate - deg_flux` to `du[tags.idx.rna]` (~lines 207-220).

---

### 2. Nonnegativity Guards (Numerical Stability)

**What:** Clamped all physically nonnegative species (tRNAs, aa-tRNAs, T3 complexes, EFT, EGD, S30, S50, S70, PIC) to `>= 0` using local copies before computing kinetic rates. Solver state vector is untouched.

**Why:** The implicit DAE solver can produce small negative concentrations during Newton iterations. A negative tRNA makes the Michaelis denominator (`KM + tRNA`) shrink, causing rates to explode and `|du|` to blow up, leading to `dt < epsilon` solver failure. This was the root cause of convergence failures with >33 background genes.

**Where:**
- `src/rhs_metabolites.jl` lines 32-75: All 39 unloaded tRNAs clamped to `>= 1e-12`; EFT, EGD to `>= 0`; S30, S50, S70, PIC to `>= 0`; all 39 aa-tRNAs and 39 T3 complexes to `>= 0`.
- `src/functions.jl` (preintegration copy) lines 735-780: Identical guards in `rhs_metabolites_preintegration!`.
- `src/rhs_sequence.jl` lines 96-100: Clamped `u[1:num_species]` to `>= 0` before computing `C_vect` (effective T3 pools for codon-specific elongation).

---

### 3. Ribosome Queuing Factor (`q_R_elo`) Clamping

**What:** Added `clamp(q_R_elo, 0.0, 1.0)` for the TASEP mean-field ribosome occupancy factor.

**Why:** At high local ribosome density the denominator of `1 - c_RIB/denom` becomes tiny, making `q_R_elo` hugely negative. This produces unphysical *reverse* elongation flux, causing runaway instability. Clamping to [0, 1] reflects the physical constraint that this is a probability.

**Where:**
- `src/rhs_sequence.jl` lines 61-65 and line 83: Both inner-codon and last-codon `q_R_elo` calculations.

---

### 4. DAE Algebraic Constraint Fix (EF-G / EF-Tu Conservation)

**What:** Changed the algebraic constraint equations for EF-G and EF-Tu mass conservation to use raw state values from `u[]` instead of the clamped local variables.

**Why:** The solver adjusts the algebraic variables to satisfy `residual = 0`. If the constraint uses clamped values (e.g., `max(EGD, 0)` instead of `u[idx.EGD]`), the solver sees a flat zero gradient when the variable is negative and cannot correct it. This caused EFT to diverge to -178,000 µM with 60+ genes.

**Where:**
- `src/rhs_metabolites.jl` lines 278-286: `du[idx.EFG] = raw_EFG + raw_EGT + raw_EGD - c_EFG_t` (was `EFG + EGT + EGD`).
- `src/functions.jl` (preintegration copy) lines 838-846: Same fix.

---

### 5. PIC Mass Balance Fix

**What:** Changed `du[tags.idx.PIC]` to use the pre-computed `V_IN_70SIC` variable instead of recomputing the initiation rate inline without the scale factor.

**Why:** The inline formula hardcoded `k_IN_70SIC * PIC / (1 + KM_mRNA/RNA + KM_50S/S50)` without `tags.scale.V_IN_70SIC_scale`, so scaling the initiation rate had no effect on PIC consumption. This made the model completely insensitive to `V_IN_70SIC_scale`.

**Where:**
- `src/rhs_metabolites.jl` line 134: `du[tags.idx.PIC] = ... - V_IN_70SIC` (was inline formula).
- `src/functions.jl` (preintegration copy) line 809: Same fix.

---

### 6. Initial Condition Scaling (`Rib_scale`, `IF_scale`, `tRNA_scale`)

**What:** Applied `Rib_scale`, `IF_scale`, and `tRNA_scale` to the corresponding initial concentrations in `calculation_of_initial_concentration`.

**Why:** These scale parameters existed in the `scale` struct but were never applied to ICs. Changing `Rib_scale` from 0.1 to 1.5 had zero effect on the simulation because the actual ribosome pools started at the same value regardless.

**Where:**
- `src/functions.jl` lines 613-621:
  ```
  y0[idx.S30] *= scale.Rib_scale
  y0[idx.S50] *= scale.Rib_scale
  y0[idx.IF1] *= scale.IF_scale
  y0[idx.IF2] *= scale.IF_scale
  y0[idx.IF3] *= scale.IF_scale
  y0[idx.b1t:idx.ty1] .*= scale.tRNA_scale
  ```

---

### 7. `call_model()` Rewrite

**What:** Rewrote `call_model()` with keyword arguments for sensitivity analysis: `custom_scale`, `growth_rate`, `mRNA_multiplier`, `het_mRNA_conc`, `background_selection` (`:first` or `:top_mRNA`), adaptive solver tolerances, `return_details` flag, and try/catch with diagnostic output on failure.

**Why:** The original function was hardcoded to `µ=1`, default scale, and first-N background selection. Running parameter sweeps and validation required manual code changes each time. The rewrite enables:
- Single-line sensitivity sweeps across any parameter
- Background gene selection by highest mRNA expression (`:top_mRNA`)
- Automatic tolerance scaling for large systems (>40 genes)
- Graceful failure reporting with RHS diagnostics

**Where:**
- `src/functions.jl` lines 1877-2130: Complete rewrite of `call_model()`.
- `src/datatypes.jl` lines 35-55: Added `make_scale(; kwargs...)` constructor.

---

### 8. `N_one` / `N_two` Clamping in `matrix_free_KT3!`

**What:** Clamped `C_vect[index] / c_codon` to `[0.0, 1.0]` and protected against division by zero with `max(c_codon, 1e-12)`.

**Why:** Negative T3 complex concentrations (numerical drift) made `c_codon ≈ 0`, causing `N_one`, `N_two` to blow up or go negative, which propagated into unphysical elongation rates.

**Where:**
- `src/functions.jl` lines 1483-1500: Both `N_one` and `N_two` loops.

---

### 9. Diagnostic Tooling

**What:** Added `src/ode_diagnostics.jl` with three print functions:
- `print_problematic_rhs(u, p, t)`: Evaluates the full RHS at a given state and prints the top-K largest `|du|` components with human-readable labels.
- `print_convergence_root_hint(u, p)`: Checks EF-Tu/EF-G specific terms (Michaelis denominator, mass conservation residual, sign of pools) and prints root-cause hints.
- `print_ef_pool_literature_audit(y, ic, p_met, idx)`: Compares EF-Tu/EF-Ts/EF-G/ribosome ratios against literature ranges (Dennis & Bremer, Neidhardt).

**Why:** Convergence failures gave cryptic solver errors (`dt < epsilon`). These functions identify *which* ODE/state is blowing up and *why*, enabling targeted fixes instead of guesswork.

**Where:**
- `src/ode_diagnostics.jl`: New file (168 lines).
- `src/packages.jl` line 35: `include("ode_diagnostics.jl")`.
- `src/main.jl` line 58: `print_ef_pool_literature_audit(...)` call after `EF_balance!`.

---

### 10. Sequence DataFrame Indexing Fix

**What:** Rewrote `sequence_dataframes()` to correctly index per-sequence output variables (S70, rna, PRT, PTC) and ribosome codon positions using cumulative offset tracking.

**Why:** The original indexing assumed a fixed layout that broke when multiple sequences had different lengths. The variable `idx` was advanced incorrectly, causing columns to contain wrong data.

**Where:**
- `src/functions.jl` lines 1545-1590: Rewritten `sequence_dataframes()`.

---

### 11. Gene Filtering (FASTA Safety)

**What:** `species_and_sequences()` now filters requested genes against what's actually in the FASTA file and warns about missing genes instead of throwing a `KeyError`.

**Why:** Removing a gene from the FASTA (e.g., `ompT`) while it was still selected from `mRNA_concentrations.xlsx` caused a `KeyError` crash.

**Where:**
- `src/functions.jl` lines 35-40: Filtered `seqs` dict and warned about skipped genes.
- `src/functions.jl` lines 1940-1960 (in `call_model`): Explicit FASTA membership check for heterologous genes; background gene pool filtered against `fasta_genes`.

---

### 12. CAI-Based Post-Hoc Correction (Validation Notebooks)

**What:** Implemented a Codon Adaptation Index (CAI) based correction for predicted protein concentrations in the Taniguchi and Valgepea validation cells of `main.ipynb`.

**Why:** The model predicts the same translation rate for any gene with the same length and mRNA level, regardless of codon composition. In reality, genes with codons matching abundant tRNAs (high CAI) are translated faster. The model systematically under-predicted protein by ~7x and lacked gene-specific variation.

**Math (log-space linear regression):**

The correction models the ratio of measured to predicted protein as a power law in CAI:

```
measured / predicted = beta * CAI^alpha
```

Taking logarithms of both sides:

```
log(measured / predicted) = alpha * log(CAI) + log(beta)
```

This is a standard linear regression `y = mx + b` where:
- `y = log(measured_prot / predicted_PRT)` for each data point
- `x = log(CAI)` for each gene
- `m = alpha` (slope) — the CAI exponent
- `b = log(beta)` (intercept) — log of the global scale factor

Fitted by ordinary least squares on the Valgepea dataset (18 genes x 5 growth rates = 90 points):

```
alpha_hat = sum((x - x_mean)(y - y_mean)) / sum((x - x_mean)^2)
log_beta  = y_mean - alpha_hat * x_mean
beta_hat  = exp(log_beta)
R^2       = 1 - SS_res / SS_tot
```

**Fitted values:** `alpha ~ 1.04`, `beta ~ 6.94`, `R^2 ~ 0.13`

**Corrected prediction:**

```
PRT_corrected = beta * CAI^alpha * PRT_predicted
```

**Interpretation:**
- `beta ~ 6.9`: The model under-predicts by ~7x on average (global calibration factor).
- `alpha ~ 1.0`: Approximately linear relationship — doubling CAI roughly doubles translational output.
- `R^2 = 0.13`: CAI explains ~13% of the residual variance (codon adaptation is one of many factors).

**CAI computation** (Sharp & Li 1987):
```
CAI = exp( (1/L) * sum_{i=1}^{L} log(w_i) )
```
where `w_i` is the relative adaptiveness of codon `i` (from E. coli K-12 codon usage table), `L` is the number of sense codons excluding Met/Trp/stop. High CAI (~0.7-0.8) = well-adapted codons; low CAI (~0.2-0.3) = rare codons.

**Results:**
- Valgepea: within-3x accuracy improved from 30% (27/90) to 91% (82/90)
- Taniguchi (tufA, CAI=0.82): raw ratio 0.26 → corrected ratio 1.49

**Where:**
- `src/main.ipynb` cell 12: Taniguchi validation with CAI functions, `ECOLI_W` table, `compute_cai()`, `read_fasta_dict()`, Valgepea-fitted `VALG_ALPHA=1.04`, `VALG_BETA=6.94`.
- `src/main.ipynb` cell 13: Taniguchi plots (parity + grouped bars) reading from CSV.
- `src/main.ipynb` cell 15: Valgepea validation with CAI computation and two-parameter log-space regression fit.
- `src/main.ipynb` cell 16: Valgepea plots (parity + grouped bars) reading from CSV.

---

### 13. Sensitivity Analysis (Notebook)

**What:** Added a focused parameter sweep in `main.ipynb` testing `k_TLEmaxj_scale`, `V_IN_70SIC_scale`, `Rib_scale`, and combined `num_bg + mRNA_multiplier` across multiple values for tufA.

**Why:** To identify which parameters the model is actually sensitive to and diagnose why background gene count had no effect on predictions.

**Where:**
- `src/main.ipynb` cell 14: Sensitivity sweep code generating `sensitivity_focused_sweep.csv` and tornado/sweep plots.

---

## Notebook Cell Map (`src/main.ipynb`)

| Cell | Type | Purpose |
|------|------|---------|
| 0 | markdown | Title / authors |
| 1 | code | Main simulation (includes, sequences, solve) |
| 2-11 | code | Plotting, results analysis, existing workflows |
| 12 | code | **Taniguchi validation** (10 genes, CAI correction) |
| 13 | code | **Taniguchi plots** (parity + grouped bars from CSV) |
| 14 | code | **Sensitivity sweep** (4 parameters) |
| 15 | code | **Valgepea validation** (18 genes x 5 growth rates, CAI fit) |
| 16 | code | **Valgepea plots** (parity + grouped bars from CSV) |

---

## Input Requirements for Simulating a New Gene

| Input | Source | Example |
|-------|--------|---------|
| Gene coding sequence | FASTA file in `data/sequences/` | Must be in-frame, ATG to stop |
| mRNA concentration | `data/input/mRNA_concentrations.xlsx` or `het_mRNA_conc` kwarg | 0.00664 µM (= 4 copies/cell) |
| Growth rate µ | Literature / experiment | 0.1 - 1.8 h⁻¹ |
| Number of background genes | User choice (≤50 recommended) | 25 |
| CAI (for post-hoc correction) | Computed from sequence via `compute_cai()` | 0.2 - 0.85 |

---

## Key References

- Taniguchi Y et al. *Science* 329:533 (2010) — mRNA & protein copies/cell for 1000 E. coli genes
- Valgepea K et al. *Mol Syst Biol* 9:689 (2013) — Growth-rate-dependent proteome/transcriptome
- Calabrese L et al. *Nature Physics* (2024) — mRNA degradation framework
- Sharp PM & Li WH. *Nucleic Acids Res* 15:1281 (1987) — Codon Adaptation Index
- Dennis PP & Bremer H. *J Bacteriol* (1974) — EF-Tu/ribosome ratios
