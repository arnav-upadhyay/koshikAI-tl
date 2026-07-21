# ============================================================================
# select_background_genes.jl  —  standalone background-gene selector
# ----------------------------------------------------------------------------
# Two modes:
#   (A) Automatic, isoacceptor-aware selection (default):
#       picks the N_BG genes whose codons rest on ABUNDANT tRNA isoacceptors
#       (ranked by the usage-weighted harmonic mean of per-codon cognate
#       supply) and steers away from the minority pools (Leu1-5, Ser5, His)
#       that crash the ODE at high N. Pure pre-screen — same source data the
#       model uses (AS_and_Codon.dat + FitCoefficients_eco_V9.xlsx::reltRNA);
#       no scale changes.
#   (B) User-supplied list:
#       pass `custom_bg = ["geneA", "geneB", ...]` to use your own background
#       genes. They are validated against the sequence source; the isoacceptor
#       scoring is skipped.
#
# Requires helpers already provided by functions.jl (readmatrix, readfasta,
# seqs_import_transcripts, to_DataFrame). Call after `include("functions.jl")`.
#
# Returns a NamedTuple: (selected_bg, df, out_fasta, seqs_all).
# ============================================================================
using DataFrames, CSV, Printf, Statistics

function select_background_genes(;
        src_fasta::AbstractString              = "Ecoli_K12_MG1655_geneset_with_welch_2009.fa",
        n_bg::Int                              = 30,
        primary_gene::AbstractString           = "rplL",
        growth_rate::Float64                   = 1.0,          # µ (h⁻¹) — sets isoacceptor abundances
        custom_bg::Union{Nothing,AbstractVector{<:AbstractString}} = nothing,  # Mode B if given
        out_fasta::Union{Nothing,AbstractString} = nothing,
        require_mrna::Bool                     = true,         # keep only genes present in mRNA xlsx
        include_primary_in_out::Bool           = true,
        slow_pctl::Float64                     = 0.25,
        select_safest::Bool                    = true,         # false = pick the WORST genes (stress test)
        fragile_names::Vector{String}          = ["Leu1","Leu2","Leu3","Leu4","Leu5","Ser5","His"],
        write_fasta::Bool                      = true,
        verbose::Bool                          = true,
    )

    # default output filename reflects the chosen mode
    if out_fasta === nothing
        out_fasta = custom_bg === nothing ? "../data/sequences/bg_safe_$(n_bg).fa" :
                                            "../data/sequences/bg_custom_$(length(custom_bg)).fa"
    end

    # ---------- load candidate sequences (needed in both modes) ----------
    nrec = length(readfasta(joinpath("..", "data", "sequences", src_fasta)))
    seqs_all, names_all = seqs_import_transcripts(src_fasta, (1, nrec))

    local selected_bg::Vector{String}
    local df::DataFrame

    if custom_bg !== nothing
        # ================= MODE B: user-supplied background genes =================
        req     = String.(custom_bg)
        missing = [g for g in req if !haskey(seqs_all, g)]
        isempty(missing) ||
            error("select_background_genes: custom gene(s) not found in $src_fasta: " *
                  join(missing, ", "))
        selected_bg = String[g for g in req if g != primary_gene]
        df = DataFrame(gene = selected_bg,
                       L    = [length(seqs_all[g].Codon) for g in selected_bg],
                       source = fill("user-supplied", length(selected_bg)))
        verbose && println("Using $(length(selected_bg)) user-supplied background genes: ",
                           join(selected_bg, ", "))
    else
        # ============ MODE A: automatic isoacceptor-aware selection ============
        # 1. isoacceptor abundance at growth rate (mirrors calculation_of_initial_concentration)
        fit  = "../data/input/FitCoefficients_eco_V9.xlsx"
        tm   = Float64.(readmatrix(fit, "Translation_machinery", "B2:D14"))
        S50  = tm[2,3] + tm[2,1] * growth_rate / (tm[2,2] + growth_rate)      # 50S µM at µ
        rel  = Float64.(readmatrix(fit, "reltRNA", "B2:D40"))                 # 39 rows: [initiator; 38 elongators]
        tRNAconc    = (rel[:,1] .* exp.(-rel[:,2] .* growth_rate) .+ rel[:,3]) .* S50
        elong_abund = tRNAconc[2:39]                                         # row 1 = initiator (b1t) -> drop

        # 2. codon -> cognate isoacceptor rows & total supply (AS_and_Codon.dat)
        ascod    = DataFrame(CSV.File("../data/input/AS_and_Codon.dat", delim = ';', missingstring = ""))
        iso_name = String.(strip.(ascod.tRNAas))
        codon_supply = Dict{String,Float64}()
        codon_isos   = Dict{String,Vector{Int}}()
        for (t, codstr) in enumerate(ascod.codAS)
            for c in split(String(codstr), ",")
                cc = strip(c); isempty(cc) && continue
                codon_supply[cc] = get(codon_supply, cc, 0.0) + elong_abund[t]
                push!(get!(codon_isos, cc, Int[]), t)
            end
        end
        FRAGILE  = Set(filter(!isnothing, [findfirst(==(n), iso_name) for n in fragile_names]))
        slow_cut = quantile(collect(values(codon_supply)), slow_pctl)

        if verbose
            println(@sprintf("S50(µ=%.2f) = %.2f µM   |   slow-codon cutoff (p%d) = %.3f µM",
                    growth_rate, S50, round(Int, 100 * slow_pctl), slow_cut))
            println("10 scarcest codons (codon => total cognate supply µM):")
            for (c, s) in sort(collect(codon_supply), by = x -> x[2])[1:10]
                println(@sprintf("   %-4s %.3f   [%s]", c, s, join(iso_name[codon_isos[c]], ",")))
            end
        end

        # 3. candidate genes
        candidates = String[g for g in names_all if g != primary_gene]
        if require_mrna
            mrna_df    = to_DataFrame(readmatrix("../data/input/mRNA_concentrations.xlsx",
                                                 "mRNA_concentration", "A1:B2324"))
            mrna_genes = Set(String.(strip.(mrna_df.gene)))
            candidates = filter(g -> g in mrna_genes, candidates)
        end

        # 4. score each candidate
        rows = NamedTuple[]
        for g in candidates
            haskey(seqs_all, g) || continue
            codons = String.(seqs_all[g].Codon)
            length(codons) < 2 && continue
            sup = Float64[]; nfrag = 0
            for c in codons[1:end-1]                              # drop terminal (stop) codon
                haskey(codon_supply, c) || continue              # skip stops / unknowns
                push!(sup, codon_supply[c])
                isempty(intersect(Set(codon_isos[c]), FRAGILE)) || (nfrag += 1)
            end
            isempty(sup) && continue
            push!(rows, (gene = g, L = length(sup),
                         harm_supply  = length(sup) / sum(1.0 ./ sup),   # usage-weighted; dominated by scarce codons
                         min_supply   = minimum(sup),
                         frac_slow    = count(<=(slow_cut), sup) / length(sup),
                         frac_fragile = nfrag / length(sup)))
        end
        df = DataFrame(rows)
        sort!(df, :harm_supply, rev = select_safest)             # safest (high harmonic supply) first
        selected_bg = df.gene[1:min(n_bg, nrow(df))]

        if verbose
            println("\nSelected $(length(selected_bg)) background genes ($(select_safest ? "SAFEST" : "WORST")):")
            show(first(df, length(selected_bg)), allrows = true, allcols = true); println()
        end
    end

    # ---------- write FASTA (primary optionally first, then background) ----------
    if write_fasta
        open(out_fasta, "w") do io
            if include_primary_in_out && haskey(seqs_all, primary_gene)
                println(io, ">", primary_gene); println(io, String(seqs_all[primary_gene].DNA))
            end
            for g in selected_bg
                println(io, ">", g); println(io, String(seqs_all[g].DNA))
            end
        end
        verbose && println("\nWrote → ", out_fasta)
    end

    return (selected_bg = selected_bg, df = df, out_fasta = out_fasta, seqs_all = seqs_all)
end
