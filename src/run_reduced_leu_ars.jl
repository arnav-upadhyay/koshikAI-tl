# ============================================================================
# run_reduced_leu_ars.jl
# ----------------------------------------------------------------------------
# Subramaniam 2013 on the reduced cell: CTG / C1A / CTA × Leu-ARS {1, 0.1}.
# Same promoter/RBS and the same mRNA dict; only the rplL CDS changes.
# Leu-RS is competitive (rhs_metabolites.jl: five isoacceptors share one V_max).
# Freeze off. Does not rebuild the background and does not touch stall_freeze.jl.
#
# Run from src/ (after make_reduced_background.jl):
#   julia --project=.. --threads=12 run_reduced_leu_ars.jl
#
# Knobs:
# ============================================================================

leu_levels   = [1.0, 0.1]          # Leu aaRS scale; the other 19 stay at 1
reporters    = ["CTG", "C1A", "CTA"]
timeout_s    = 36000               # warn if a finished solve took longer (does not kill FBDF)
ode_reltol   = 1e-8
ode_abstol   = 1e-4                # same as the YFP-only CDS sweep
ode_maxiters = 5_000_000
t_end        = 100.0               # min, harvest
resume       = true                # skip rows already marked ok in the CSV
load_csv     = "../data/output/reduced_cell_load.csv"
yfp_subdir   = "reduced_yfp"
out_csv      = "../data/output/reduced_yfp_leu_ars_competitive.csv"

# Leave false. Stall freezing is a separate plan and is not part of this test.
freeze_stalled = false

# ============================================================================

println("run_reduced_leu_ars: loading model…")
flush(stdout)

include("packages.jl")
using CSV, DataFrames, Printf

freeze_stalled == false ||
    error("run_reduced_leu_ars: freeze_stalled must stay false for this sweep")

isfile(load_csv) ||
    error("run_reduced_leu_ars: missing $load_csv — run make_reduced_background.jl first")

load = CSV.read(load_csv, DataFrame)
ncol(load) > 0 || error("run_reduced_leu_ars: $load_csv is empty")
genes = String.(load.gene)
concs = Float64.(load.m_assigned)
any(isnan, concs) && error("run_reduced_leu_ars: NaN in m_assigned")
reporter = genes[1]
het_conc = Dict{String,Float64}(zip(genes, concs))

function _headers(path)
    hs = String[]
    for line in eachline(path)
        startswith(line, ">") || continue
        push!(hs, String(split(chop(line; head = 1, tail = 0))[1]))
    end
    return hs
end

fastas = Dict{String,String}()
for lab in reporters
    rel = joinpath(yfp_subdir, "bg_$(lab).fa")
    path = joinpath(sequence_data_dir(), rel)
    isfile(path) || error("run_reduced_leu_ars: missing $path")
    hs = _headers(path)
    hs == genes ||
        error("run_reduced_leu_ars: $(lab) FASTA gene order does not match $load_csv")
    fastas[lab] = rel
end

function _yfp_cta_sites()
    pnas = _pnas_sd02_path()
    recs = _read_fa_ordered(pnas)
    getrec(p) = last(recs[findfirst(r -> startswith(first(r), p), recs)])
    yfp0, yfp1 = getrec("yfp0"), getrec("yfp1")
    cods0 = [yfp0[i:i+2] for i in 1:3:length(yfp0)]
    cods1 = [yfp1[i:i+2] for i in 1:3:length(yfp1)]
    length(cods0) == length(cods1) ||
        error("run_reduced_leu_ars: yfp0 and yfp1 differ in length")
    sens = findall(k -> cods0[k] != cods1[k], eachindex(cods0))
    isempty(sens) && error("run_reduced_leu_ars: yfp0 and yfp1 have no codon differences")
    return (
        CTG = Int[],
        C1A = Int[sens[1]],
        CTA = collect(sens),
    )
end

const CTA_SITES = _yfp_cta_sites()
const LEU_IDX = findfirst(==('L'), AA20_LETTERS)
LEU_IDX === nothing && error("run_reduced_leu_ars: leucine missing from AA20_LETTERS")

# Same scale as the CDS-only notebook sweep: initiation and every aaRS at 1,
# then leucine moved to `la`. Not the make_scale default (ARS = 20, V_IN = 75).
_make_scale = txtl_types.make_scale
function mk_scale(la)
    ars = ones(20)
    ars[LEU_IDX] = la
    return _make_scale(
        V_TC_scale = 1.0, k_TLEmaxj_scale = 1.0,
        V_IN_70SIC_scale = 1.0, V_IN_IF2D_scale = 1.0,
        Rib_scale = 1.0, IF_scale = 1.0, EF_scale = 1.0, tRNA_scale = 1.0,
        ARS_scale = ars, RibRecycling = 0.0, V_mDeg_scale = 0.0)
end

function decode_leu_pools(payload, gene, te)
    p = payload.p
    si = _find_seq_index(p, gene)
    tg = p.local_tags_list[si]
    us = extract_u_sequence(payload.sol(te), p, si)
    ns = tg.num_species
    cC = tg.t3_Codon_struc * reshape(max.(us[1:ns], 0.0), 1, ns)'
    ci = tg.species_and_sequences.cod_idx
    poolof(target) = begin
        for (k, v) in ci
            string(k) == target && return cC[v]
        end
        return NaN
    end
    return (cta = poolof("CUA"), ctg = poolof("CUG"))
end

mkpath(dirname(out_csv))
rows = NamedTuple[]
if resume && isfile(out_csv)
    prev = CSV.read(out_csv, DataFrame)
    for r in eachrow(prev)
        push!(rows, (
            label = String(r.label),
            leu_ars = Float64(r.leu_ars),
            PRT = Float64(r.PRT),
            fold_vs_CTG = Float64(r.fold_vs_CTG),
            RIB = Float64(r.RIB),
            cta_tau = Float64(r.cta_tau),
            cta_pool = Float64(r.cta_pool),
            ctg_pool = Float64(r.ctg_pool),
            T3_min = Float64(r.T3_min),
            wall_s = Float64(r.wall_s),
            n_seq = Int(r.n_seq),
            status = String(r.status),
        ))
    end
    println("Resume: ", count(r -> r.status == "ok", rows), " ok rows in ", out_csv)
end

function _done(lab, la)
    any(r -> r.status == "ok" && r.label == lab && isapprox(r.leu_ars, la; atol = 1e-9), rows)
end

function _save()
    CSV.write(out_csv, DataFrame(rows))
    flush(stdout)
end

@printf("Reduced Leu-ARS sweep   reporter=%s   n_seq=%d   m_reporter=%.4f µM   m_total=%.4f µM\n",
        reporter, length(genes), het_conc[reporter], sum(values(het_conc)))
println("  freeze_stalled = false")
println("  FASTAs: ", join((fastas[lab] for lab in reporters), ", "))
println("  CTA codon sites: C1A=$(CTA_SITES.C1A)  CTA=$(CTA_SITES.CTA)")
flush(stdout)

for la in leu_levels, lab in reporters
    if resume && _done(lab, la)
        @printf("\n%-12s (Leu-ARS=%.3g) skip (already ok)\n", lab, la)
        flush(stdout)
        continue
    end
    @printf("\n%-12s (Leu-ARS=%.3g) %s\n", lab, la, Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    flush(stdout)
    t0 = time()
    status = "error"
    payload = nothing
    try
        payload = call_model(fastas[lab], genes, 0;
            custom_scale = mk_scale(la),
            growth_rate = 1.0,
            het_mRNA_conc = het_conc,
            return_details = true,
            return_sol = true,
            static_EFs = true,
            spindler_75pct = false,
            freeze_stalled = false,
            t_end = t_end,
            ode_reltol = ode_reltol,
            ode_abstol = ode_abstol,
            ode_maxiters = ode_maxiters)
        status = "ok"
    catch e
        status = "error"
        println("FAILED")
        showerror(stdout, e)
        println()
        flush(stdout)
    end
    wall = time() - t0
    if wall > timeout_s
        @warn "solve exceeded timeout_s" label = lab leu_ars = la wall_s = wall timeout_s
    end

    ok = status == "ok" && payload !== nothing &&
         hasproperty(payload, :sol) && successful_retcode(payload.sol)
    if !ok
        st = status == "ok" ? "retcode_$(payload.sol.retcode)" : status
        push!(rows, (
            label = lab, leu_ars = la, PRT = NaN, fold_vs_CTG = NaN, RIB = NaN,
            cta_tau = NaN, cta_pool = NaN, ctg_pool = NaN, T3_min = NaN,
            wall_s = wall, n_seq = length(genes), status = st,
        ))
        @printf("  %s  wall=%.0fs\n", st, wall)
        _save()
        continue
    end

    te = float(payload.sol.t[end])
    dg = collect_pool_diagnostics(payload.sol, payload.p, te; primary_gene = reporter)
    sp = spindler_codon_translation_time(payload.sol, payload.p, reporter; t_sample = te)
    dp = decode_leu_pools(payload, reporter, te)
    cs = filter(p -> 1 <= p <= sp.L, CTA_SITES[Symbol(lab)])
    cta_tau = isempty(cs) ? 0.0 : sum(sp.tau[cs])
    push!(rows, (
        label = lab, leu_ars = la, PRT = dg.PRT_primary, fold_vs_CTG = NaN,
        RIB = dg.RIB_active, cta_tau = cta_tau, cta_pool = dp.cta, ctg_pool = dp.ctg,
        T3_min = dg.T3_min, wall_s = wall, n_seq = length(genes), status = "ok",
    ))
    @printf("  ok  PRT=%.4g  RIB=%.4g  ctaτ=%.3g  CTApool=%.3g  CTGpool=%.3g  T3min=%.3g  wall=%.0fs\n",
            dg.PRT_primary, dg.RIB_active, cta_tau, dp.cta, dp.ctg, dg.T3_min, wall)
    _save()
end

res = DataFrame(rows)
for la in leu_levels
    hit = findfirst(i -> res.label[i] == "CTG" && isapprox(res.leu_ars[i], la; atol = 1e-9) &&
                         res.status[i] == "ok", 1:nrow(res))
    hit === nothing && continue
    p_ctg = res.PRT[hit]
    p_ctg == 0 && continue
    for i in 1:nrow(res)
        isapprox(res.leu_ars[i], la; atol = 1e-9) || continue
        res.status[i] == "ok" || continue
        res.fold_vs_CTG[i] = res.PRT[i] / p_ctg
    end
end
CSV.write(out_csv, res)
println("\nSaved → ", out_csv)
show(res, allrows = true, allcols = true)
println()

println("\nPaper check (Subramaniam 2013 Fig. 1 — same promoter/RBS, CDS only, cell-like Leu demand):")
println("  rich (Leu-ARS=1):   PRT(CTG) ≈ PRT(C1A) ≈ PRT(CTA)")
println("  starved (0.1):      PRT(CTA)/PRT(CTG) ≪ 1; C1A in between")
for la in leu_levels
    println("\n  Leu-ARS = ", la)
    for r in eachrow(filter(row -> isapprox(row.leu_ars, la; atol = 1e-9), res))
        @printf("    %-8s  PRT=%8.4g  vs CTG = %6.3f  ctaτ=%.3g  status=%s\n",
                r.label, r.PRT, r.fold_vs_CTG, r.cta_tau, r.status)
    end
end
println("\nStarved / rich PRT fold (same CDS):")
for lab in reporters
    rich = filter(row -> row.label == lab && isapprox(row.leu_ars, 1.0; atol = 1e-9) && row.status == "ok", res)
    star = filter(row -> row.label == lab && isapprox(row.leu_ars, 0.1; atol = 1e-9) && row.status == "ok", res)
    (nrow(rich) == 0 || nrow(star) == 0) && continue
    @printf("  %-3s  PRT(0.1)/PRT(1) = %.3g\n", lab, star.PRT[1] / rich.PRT[1])
end
flush(stdout)
