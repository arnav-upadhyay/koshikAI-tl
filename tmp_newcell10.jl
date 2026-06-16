# ═══════════════════════════════════════════════════════════════════════════════
#  Goodman 2013 rplS-variant batch — plots, ranking & RBS / CDS-type analysis
#  Reads the batch CSVs produced by the simulation cell and joins them with the
#  Goodman 2013 experimental table (Trans, RBS.Display, CDS.type, Promoter).
#
#  ── KEY MODEL FACT (read before interpreting the RBS panels) ──────────────────
#  The spliced mRNA = variant CDS head (11 codons) + WT rplS body. The RBS and
#  promoter are NOT part of the modeled sequence, so two variants that share a CDS
#  head but differ only in RBS produce IDENTICAL model output. The model is
#  therefore structurally blind to RBS strength. The RBS panels below do NOT test
#  "does the model reproduce the RBS effect" (it cannot) — they test whether the
#  experimental Trans is better *explained by elongation* within a given RBS class.
# ═══════════════════════════════════════════════════════════════════════════════
include("packages.jl")
include("goodman_2013.jl")          # spearman_rho_g, coerce_float!
using DataFrames, CSV, Plots, Printf, Statistics

# ─── KNOBS ───────────────────────────────────────────────────────────────────────
N_PLOT      = 1                      # which batch to plot (1 / 5 / 10 …) → picks CSV pair
BG_MODE     = "first"
RIB_MAX_UM  = 65.0                   # physical-validity cutoff: drop runs whose RIB_active
                                     #   exceeds this (healthy cluster ≈ 50 µM; > 65 ⇒ mass-
                                     #   conservation breakdown / unphysical basin).
out_fig_dir = "../data/output/figures"
mkpath(out_fig_dir)

base         = "../data/output/goodman_rplS_N$(N_PLOT)_$(BG_MODE)_batch"
codon_path   = "$(base)_codon_profile.csv"
summary_path = "$(base)_summary.csv"
tag_plot     = "goodman_rplS_N$(N_PLOT)_$(BG_MODE)"

# ─── 1. Load batch artifacts ─────────────────────────────────────────────────────
batch_codon   = CSV.read(codon_path,   DataFrame)
batch_summary = CSV.read(summary_path, DataFrame)
batch_summary.variant = String.(batch_summary.variant)
batch_codon.variant   = String.(batch_codon.variant)

# ─── 2. Join experimental Trans + RBS + CDS.type from Goodman 2013 ───────────────
goodman_df = CSV.read(joinpath("data", "Goodman_2013_Science.csv"), DataFrame;
                      missingstring = ["", "NA", "NaN", "#N/A"])
coerce_float!(goodman_df, [:Trans])
exp_trans = Dict{String,Float64}(); exp_rbs = Dict{String,String}()
exp_cds   = Dict{String,String}();  exp_prom = Dict{String,String}()
for r in eachrow(goodman_df)
    nm = String(r.Name)
    exp_trans[nm] = ismissing(r.Trans) ? NaN : Float64(r.Trans)
    exp_rbs[nm]   = ismissing(r."RBS.Display") ? "?" : String(r."RBS.Display")
    exp_cds[nm]   = ismissing(r."CDS.type")    ? "?" : String(r."CDS.type")
    exp_prom[nm]  = ismissing(r.Promoter)      ? "?" : String(r.Promoter)
end
batch_summary.Trans_exp = [get(exp_trans, v, NaN) for v in batch_summary.variant]
batch_summary.RBS       = [get(exp_rbs,  v, "?") for v in batch_summary.variant]
batch_summary.CDS_type  = [get(exp_cds,  v, "?") for v in batch_summary.variant]
batch_summary.Promoter  = [get(exp_prom, v, "?") for v in batch_summary.variant]

# ─── 3. Filter to physically valid converged runs ───────────────────────────────
n_total    = nrow(batch_summary)
is_ok      = batch_summary.status .== "ok"
has_v      = .!isnan.(batch_summary.rplL_mean_v_cs)
phys       = (batch_summary.RIB_active .<= RIB_MAX_UM) .& .!isnan.(batch_summary.RIB_active)
valid_mask = is_ok .& has_v .& phys
valid      = batch_summary[valid_mask, :]
sort!(valid, :Trans_exp)
n_timeout  = count(.!is_ok)
n_unphys   = count(is_ok .& has_v .& .!phys)
n_unique_v = length(unique(round.(valid.rplL_mean_v_cs, digits=4)))

println("─"^80)
@printf("Variants: %d total | %d ok | %d timed-out | %d ok-but-unphysical (RIB>%.0f µM) | %d VALID\n",
        n_total, count(is_ok), n_timeout, n_unphys, RIB_MAX_UM, nrow(valid))
@printf("Distinct model mean_v values among %d valid variants: %d  (duplicates ⇒ RBS-blind CDS sharing)\n",
        nrow(valid), n_unique_v)
println("─"^80)
if nrow(valid) < 3
    error("Fewer than 3 valid variants — nothing meaningful to plot. Loosen RIB_MAX_UM or rerun more variants.")
end

# Codon rows restricted to valid variants.
valid_set = Set(valid.variant)
cod = batch_codon[in.(batch_codon.variant, Ref(valid_set)), :]

# ─── Colour helpers ──────────────────────────────────────────────────────────────
trans_of = Dict(r.variant => r.Trans_exp for r in eachrow(valid))
lt   = [log10(trans_of[v]) for v in keys(trans_of)]
loL, hiL = minimum(lt), maximum(lt)
cmap = cgrad(:viridis)
vcolor(v) = cmap[(log10(trans_of[v]) - loL) / max(hiL - loL, eps())]

rbs_types   = sort(unique(valid.RBS))
cds_types   = sort(unique(valid.CDS_type))
cat_palette = [:steelblue, :firebrick, :seagreen, :goldenrod, :purple, :teal, :brown, :magenta]
rbs_color   = Dict(t => cat_palette[mod1(i, length(cat_palette))] for (i, t) in enumerate(rbs_types))
cds_color   = Dict(t => cat_palette[mod1(i, length(cat_palette))] for (i, t) in enumerate(cds_types))

# ═════════════════════════════════════════════════════════════════════════════════
#  FIGURE 1 — smoothed elongation rate v(codon) and dwell-time τ, coloured by Trans
# ═════════════════════════════════════════════════════════════════════════════════
pA1 = plot(title  = "Smoothed elongation rate v(codon) — $(nrow(valid)) valid rplS variants (N=$(N_PLOT))",
           xlabel = "Codon position", ylabel = "v_smooth (codons/s)",
           legend = false, size = (1250, 560),
           left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, right_margin = 14Plots.mm)
for v in unique(cod.variant)
    sub = cod[cod.variant .== v, :]
    plot!(pA1, sub.codon_pos, sub.v_smooth_cod_s; color = vcolor(v), lw = 1.2, alpha = 0.75)
end
# dummy series to force a Trans colourbar
scatter!(pA1, [NaN], [NaN]; zcolor = [loL], clims = (loL, hiL), c = :viridis,
         markersize = 0, label = "", colorbar_title = "log10(Trans exp)")

pA2 = plot(title  = "Codon dwell time τ = 1/v  (log y)",
           xlabel = "Codon position", ylabel = "τ (s)", yscale = :log10,
           legend = false, size = (1250, 560),
           left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, right_margin = 14Plots.mm)
for v in unique(cod.variant)
    sub = cod[cod.variant .== v, :]
    plot!(pA2, sub.codon_pos, sub.tau_s; color = vcolor(v), lw = 1.2, alpha = 0.75)
end
scatter!(pA2, [NaN], [NaN]; zcolor = [loL], clims = (loL, hiL), c = :viridis,
         markersize = 0, label = "", colorbar_title = "log10(Trans exp)")

figA = plot(pA1, pA2; layout = (2, 1), size = (1250, 1050))
savefig(figA, "$(out_fig_dir)/$(tag_plot)_elongation_overlay.png")
display(figA)

# ═════════════════════════════════════════════════════════════════════════════════
#  FIGURE 2 — predicted (sim) vs experimental Trans, coloured by RBS type
# ═════════════════════════════════════════════════════════════════════════════════
function scatter_vs_trans(yvec, ylab, ttl; logy = false)
    x = Float64.(valid.Trans_exp); y = Float64.(yvec)
    ρ = spearman_rho_g(x, y)
    r = cor(x, y)
    p = plot(title = @sprintf("%s\nSpearman ρ = %+.2f   Pearson r = %+.2f   (n=%d)", ttl, ρ, r, length(x)),
             xlabel = "Experimental Trans (Prot/RNA)", ylabel = ylab,
             xscale = :log10, legend = :outertopright, size = (640, 520),
             left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
    logy && plot!(p; yscale = :log10)
    for t in rbs_types
        idx = valid.RBS .== t
        any(idx) || continue
        scatter!(p, x[idx], y[idx]; label = "$t (n=$(count(idx)))",
                 color = rbs_color[t], markersize = 7, markerstrokewidth = 0.6)
    end
    return p, ρ, r
end

p_v,   ρv,   rv   = scatter_vs_trans(valid.rplL_mean_v_cs,    "mean v (codons/s)",  "Mean elongation rate")
p_tau, ρtau, rtau = scatter_vs_trans(valid.rplL_total_time_s, "Σ τ (s)",            "Total translation time"; logy = true)
p_prt, ρprt, rprt = scatter_vs_trans(valid.PRT_rplL,          "PRT_rplL (µM)",      "Protein output")

figB = plot(p_v, p_tau, p_prt; layout = (1, 3), size = (1950, 540))
savefig(figB, "$(out_fig_dir)/$(tag_plot)_pred_vs_exp_byRBS.png")
display(figB)

# ═════════════════════════════════════════════════════════════════════════════════
#  Per-group rank-correlation helper (used for both RBS and CDS.type)
# ═════════════════════════════════════════════════════════════════════════════════
function group_rho_table(groupcol::Symbol, group_label::String)
    groups = sort(unique(valid[!, groupcol]))
    rows = NamedTuple[]
    for g in groups
        sub = valid[valid[!, groupcol] .== g, :]
        n   = nrow(sub)
        ρv  = n >= 3 ? spearman_rho_g(sub.Trans_exp, sub.rplL_mean_v_cs)    : NaN
        ρt  = n >= 3 ? spearman_rho_g(sub.Trans_exp, sub.rplL_total_time_s) : NaN
        push!(rows, (group = string(g), n = n, rho_meanv = ρv, rho_tau = ρt))
    end
    df = DataFrame(rows)
    println("\n  Per-$(group_label) Spearman ρ (experimental Trans vs model output):")
    @printf("    %-12s %4s %14s %14s\n", group_label, "n", "ρ(Trans,meanv)", "ρ(Trans,Στ)")
    for r in eachrow(df)
        @printf("    %-12s %4d %14s %14s\n", r.group, r.n,
                isnan(r.rho_meanv) ? "  n/a (<3)" : @sprintf("%+.3f", r.rho_meanv),
                isnan(r.rho_tau)   ? "  n/a (<3)" : @sprintf("%+.3f", r.rho_tau))
    end
    return df
end

# ═════════════════════════════════════════════════════════════════════════════════
#  FIGURE 3 — RBS-type analysis
#    (a) mean_v distribution by RBS type (strip plot)
#    (b) per-RBS Spearman ρ bars: which RBS class is best "explained" by elongation
# ═════════════════════════════════════════════════════════════════════════════════
rbs_rho = group_rho_table(:RBS, "RBS")

# (a) strip plot of mean_v by RBS
pR1 = plot(title = "mean v by RBS type  (model is RBS-blind ⇒ spread = CDS content only)",
           xlabel = "RBS type", ylabel = "mean v (codons/s)", legend = false,
           xticks = (1:length(rbs_types), rbs_types), size = (760, 520),
           left_margin = 8Plots.mm, bottom_margin = 10Plots.mm)
for (i, t) in enumerate(rbs_types)
    sub = valid[valid.RBS .== t, :]
    xs  = fill(i, nrow(sub)) .+ 0.12 .* randn(nrow(sub))
    scatter!(pR1, xs, sub.rplL_mean_v_cs; color = rbs_color[t], markersize = 6,
             markerstrokewidth = 0.5)
end

# (b) per-RBS ρ bars
rbs_ok = rbs_rho[.!isnan.(rbs_rho.rho_meanv), :]
pR2 = bar(rbs_ok.group, rbs_ok.rho_meanv;
          title = "Per-RBS ρ(Trans, mean v)\n(closer to +1 ⇒ model tracks experiment better)",
          ylabel = "Spearman ρ", legend = false, xrotation = 20,
          color = [rbs_color[g] for g in rbs_ok.group], ylims = (-1, 1),
          size = (760, 520), left_margin = 8Plots.mm, bottom_margin = 12Plots.mm)
hline!(pR2, [0.0]; color = :black, lw = 1, ls = :dash)

figR = plot(pR1, pR2; layout = (1, 2), size = (1550, 540))
savefig(figR, "$(out_fig_dir)/$(tag_plot)_RBS_analysis.png")
display(figR)

# ═════════════════════════════════════════════════════════════════════════════════
#  FIGURE 4 — CDS.type analysis (this is what actually varies in the model)
# ═════════════════════════════════════════════════════════════════════════════════
cds_rho = group_rho_table(:CDS_type, "CDS.type")

pC1 = plot(title = "mean v by CDS.type  (the recoding the model DOES see)",
           xlabel = "CDS.type", ylabel = "mean v (codons/s)", legend = false,
           xticks = (1:length(cds_types), cds_types), size = (820, 520),
           xrotation = 20, left_margin = 8Plots.mm, bottom_margin = 12Plots.mm)
for (i, t) in enumerate(cds_types)
    sub = valid[valid.CDS_type .== t, :]
    xs  = fill(i, nrow(sub)) .+ 0.12 .* randn(nrow(sub))
    scatter!(pC1, xs, sub.rplL_mean_v_cs; color = cds_color[t], markersize = 6,
             markerstrokewidth = 0.5)
end

cds_ok = cds_rho[.!isnan.(cds_rho.rho_meanv), :]
pC2 = bar(cds_ok.group, cds_ok.rho_meanv;
          title = "Per-CDS.type ρ(Trans, mean v)",
          ylabel = "Spearman ρ", legend = false, xrotation = 20,
          color = [cds_color[g] for g in cds_ok.group], ylims = (-1, 1),
          size = (820, 520), left_margin = 8Plots.mm, bottom_margin = 12Plots.mm)
hline!(pC2, [0.0]; color = :black, lw = 1, ls = :dash)

figC = plot(pC1, pC2; layout = (1, 2), size = (1650, 540))
savefig(figC, "$(out_fig_dir)/$(tag_plot)_CDStype_analysis.png")
display(figC)

# ═════════════════════════════════════════════════════════════════════════════════
#  PRINTOUT — overall correlations, ranking, verdict
# ═════════════════════════════════════════════════════════════════════════════════
println("\n" * "█"^80)
println("  GOODMAN rplS  —  N=$(N_PLOT), :$(BG_MODE)  —  model vs experiment  ($(nrow(valid)) valid variants)")
println("█"^80)

@printf("\nOverall Spearman ρ (experimental Trans vs model output):\n")
@printf("    mean v (cod/s)   ρ = %+.3f    (POSITIVE expected if model matches)\n", ρv)
@printf("    Σ τ (s)          ρ = %+.3f    (NEGATIVE expected if model matches)\n", ρtau)
@printf("    PRT_rplL (µM)    ρ = %+.3f    (POSITIVE expected if model matches)\n", ρprt)

valid_sortv = sort(valid, :rplL_mean_v_cs; rev = true)
valid_sortt = sort(valid, :Trans_exp;      rev = true)
println("\n>>> Top 5 by SIMULATION (fastest mean v):")
@printf("    %-22s %12s %10s %10s %-8s %-10s\n", "variant", "Trans(exp)", "mean v", "Στ(s)", "RBS", "CDS.type")
for r in eachrow(first(valid_sortv, 5))
    @printf("    %-22s %12.0f %10.3f %10.3f %-8s %-10s\n",
            r.variant, r.Trans_exp, r.rplL_mean_v_cs, r.rplL_total_time_s, r.RBS, r.CDS_type)
end
println("\n>>> Top 5 by EXPERIMENT (highest Trans):")
for r in eachrow(first(valid_sortt, 5))
    @printf("    %-22s %12.0f %10.3f %10.3f %-8s %-10s\n",
            r.variant, r.Trans_exp, r.rplL_mean_v_cs, r.rplL_total_time_s, r.RBS, r.CDS_type)
end

# Which RBS / CDS class does the model track best?
function best_group(df)
    sub = df[.!isnan.(df.rho_meanv), :]
    isempty(sub) && return ("n/a", NaN)
    i = argmax(sub.rho_meanv)
    return (sub.group[i], sub.rho_meanv[i])
end
bestRBS, bestRBSρ = best_group(rbs_rho)
bestCDS, bestCDSρ = best_group(cds_rho)
println("\nBest-tracked subgroup (highest ρ between Trans and model mean v):")
@printf("    RBS class      : %-10s  ρ = %+.3f\n", bestRBS, bestRBSρ)
@printf("    CDS.type class : %-10s  ρ = %+.3f\n", bestCDS, bestCDSρ)

println("\nFigures saved:")
println("    $(out_fig_dir)/$(tag_plot)_elongation_overlay.png")
println("    $(out_fig_dir)/$(tag_plot)_pred_vs_exp_byRBS.png")
println("    $(out_fig_dir)/$(tag_plot)_RBS_analysis.png")
println("    $(out_fig_dir)/$(tag_plot)_CDStype_analysis.png")
println("\nNOTE: the model is RBS-blind (RBS/promoter are not in the spliced mRNA). Identical")
println("model mean_v for different variants ⇒ those variants share a CDS head but differ only")
println("in RBS. Per-RBS ρ tells you which RBS class's Trans is best explained by elongation.")
println("█"^80)
