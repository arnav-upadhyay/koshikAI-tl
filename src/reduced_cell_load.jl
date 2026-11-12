# ============================================================================
# reduced_cell_load.jl  —  mini-transcriptome with representative mRNA masses
# ----------------------------------------------------------------------------
# Compresses the tabulated E. coli mRNA table into N ≈ 20 genes:
#   HEAD  = the real highest-mRNA genes (keep their xlsx µM)
#   TAIL  = evenly spaced ranks from the remainder; they carry the **lumped**
#           leftover transcriptome mass (sum of every gene not in the head)
#
# That preserves (1) abundance skew and (2) total codon-demand mass without
# integrating 2000 TASEP tracks. The reporter is a small extra fraction, not
# half the cell (Subramaniam: low-copy YFP in a real proteome).
#
# Where to change concentrations (kwargs, or the knobs in make_reduced_background.jl):
#   n_head, n_tail     — how many high vs low tracks (default 10 + 10)
#   f_reporter         — reporter share of the final budget (default 0.05)
#   m_total            — if set, rescale everything to this µM; if `nothing`,
#                        background mass = sum of the xlsx proteome (native)
#   mRNA_multiplier    — extra global scale after the budget is built
#   tail_mode          — :codon_match (default) picks the tail so leftover
#                        expression-weighted codon demand is preserved;
#                        :even_rank spaces genes through the remainder
#   src_fasta          — sequence source (full geneset if you have it; else
#                        the largest bg_safe_*.fa in data/sequences)
#
# Pass the result into call_model as:
#   call_model(cell.out_fasta, cell.het_genes, 0; het_mRNA_conc = cell.het_mRNA_conc, ...)
# Scientist entry point: julia --project=. src/make_reduced_background.jl  (cwd = repo root
# or run from src/ with --project=..)
# ============================================================================

"""
    reduced_cell_load(; n_head=10, n_tail=10, reporter="rplL", kwargs...)

Returns a NamedTuple:
- `het_genes`        reporter first, then head, then tail (pass as heterologous_genes)
- `background_genes` head ∪ tail
- `het_mRNA_conc`    Dict gene → µM  (the knob `call_model` actually uses)
- `df`               per-gene table (class, rank, m_xlsx, m_assigned, L)
- `m_total`, `m_reporter`, `m_head`, `m_tail`
- `out_fasta`        written under data/sequences unless `write_fasta=false`
"""
function reduced_cell_load(;
        n_head::Int = 10,
        n_tail::Int = 10,
        reporter::AbstractString = "rplL",
        f_reporter::Float64 = 0.05,
        m_total::Union{Nothing,Float64} = nothing,
        mRNA_multiplier::Float64 = 1.0,
        src_fasta::Union{Nothing,AbstractString} = nothing,
        mrna_xlsx::AbstractString = joinpath(@__DIR__, "..", "data", "input", "mRNA_concentrations.xlsx"),
        mrna_sheet::AbstractString = "mRNA_concentration",
        mrna_range::AbstractString = "A1:B2402",
        out_fasta::Union{Nothing,AbstractString} = nothing,
        write_fasta::Bool = true,
        include_reporter_in_fasta::Bool = true,
        tail_mode::Symbol = :codon_match,
        verbose::Bool = true,
    )

    n_head >= 0 && n_tail >= 0 || error("reduced_cell_load: n_head and n_tail must be ≥ 0")
    0.0 <= f_reporter < 1.0 || error("reduced_cell_load: need 0 ≤ f_reporter < 1")
    reporter = String(reporter)

    fasta_name = src_fasta === nothing ? _default_src_fasta() : String(src_fasta)
    nrec = length(readfasta(joinpath(sequence_data_dir(), fasta_name)))
    seqs_all, _ = seqs_import_transcripts(fasta_name, (1, nrec))
    fasta_names = Set(String.(keys(seqs_all)))

    mrna = to_DataFrame(readmatrix(mrna_xlsx, mrna_sheet, mrna_range))
    mrna = filter(row -> isa(row.mRNA, Number), mrna)
    mrna.gene = String.(strip.(string.(mrna.gene)))
    mrna.mRNA = Float64.(mrna.mRNA)
    mrna = filter(row -> row.gene in fasta_names && row.gene != reporter && row.mRNA > 0, mrna)
    nrow(mrna) > 0 || error("reduced_cell_load: no overlapping genes between $fasta_name and the mRNA table")
    sort!(mrna, :mRNA, rev = true)

    n_head_e = min(n_head, nrow(mrna))
    rest = mrna[(n_head_e + 1):end, :]
    n_tail_e = min(n_tail, nrow(rest))
    head = mrna[1:n_head_e, :]
    tail = if tail_mode == :even_rank
        _even_rank_sample(rest, n_tail_e)
    elseif tail_mode == :codon_match
        _codon_match_tail(rest, n_tail_e, seqs_all)
    else
        error("reduced_cell_load: tail_mode must be :codon_match or :even_rank, got $(repr(tail_mode))")
    end

    m_xlsx_all_rest = nrow(rest) == 0 ? 0.0 : sum(rest.mRNA)
    m_head_xlsx = n_head_e == 0 ? 0.0 : sum(head.mRNA)
    tail_w = n_tail_e == 0 ? Float64[] : max.(tail.mRNA, 0.0)
    wsum = sum(tail_w)
    m_tail_each = if n_tail_e == 0
        Float64[]
    elseif wsum <= 0
        fill(m_xlsx_all_rest / n_tail_e, n_tail_e)
    else
        m_xlsx_all_rest .* (tail_w ./ wsum)
    end
    m_head_each = n_head_e == 0 ? Float64[] : Vector{Float64}(head.mRNA)

    m_bg = m_head_xlsx + m_xlsx_all_rest
    if m_total !== nothing
        budget_bg = m_total * (1.0 - f_reporter)
        scale_bg = m_bg > 0 ? budget_bg / m_bg : 0.0
        m_head_each .*= scale_bg
        m_tail_each .*= scale_bg
        m_rep = m_total * f_reporter
        m_bg *= scale_bg
    else
        m_rep = m_bg > 0 ? f_reporter / (1.0 - f_reporter) * m_bg : 0.0
    end

    m_head_each .*= mRNA_multiplier
    m_tail_each .*= mRNA_multiplier
    m_rep *= mRNA_multiplier
    m_bg *= mRNA_multiplier

    background_genes = vcat(Vector{String}(head.gene), Vector{String}(tail.gene))
    het_genes = vcat([reporter], background_genes)
    conc = Dict{String,Float64}(reporter => m_rep)
    for (g, m) in zip(head.gene, m_head_each); conc[String(g)] = m; end
    for (g, m) in zip(tail.gene, m_tail_each); conc[String(g)] = m; end

    df = DataFrame(
        gene = String[],
        class = String[],
        rank_xlsx = Int[],
        L = Int[],
        n_CTG = Int[],
        n_CTA = Int[],
        m_xlsx = Float64[],
        m_assigned = Float64[],
    )
    function _len(g)
        haskey(seqs_all, g) ? length(seqs_all[g].Codon) : 0
    end
    function _leu(g)
        haskey(seqs_all, g) ? _leu_codon_counts(seqs_all[g]) : (0, 0)
    end
    nC, nA = haskey(seqs_all, reporter) ? _leu(reporter) : (0, 0)
    push!(df, (reporter, "reporter", 0, haskey(seqs_all, reporter) ? _len(reporter) : 0,
               nC, nA, NaN, m_rep))
    for i in 1:n_head_e
        g = String(head.gene[i])
        nC, nA = _leu(g)
        push!(df, (g, "head", i, _len(g), nC, nA, head.mRNA[i], m_head_each[i]))
    end
    rank0 = n_head_e
    for i in 1:n_tail_e
        g = String(tail.gene[i])
        r = findfirst(==(g), mrna.gene)
        nC, nA = _leu(g)
        push!(df, (g, "tail", r === nothing ? rank0 + i : Int(r), _len(g), nC, nA,
                   tail.mRNA[i], m_tail_each[i]))
    end

    if out_fasta === nothing
        out_fasta = "reduced_cell_h$(n_head_e)_t$(n_tail_e).fa"
    end
    out_path = joinpath(sequence_data_dir(), out_fasta)
    if write_fasta
        open(out_path, "w") do io
            if include_reporter_in_fasta && haskey(seqs_all, reporter)
                println(io, ">", reporter)
                println(io, String(seqs_all[reporter].DNA))
            end
            for g in background_genes
                println(io, ">", g)
                println(io, String(seqs_all[g].DNA))
            end
        end
    end

    m_head = sum(m_head_each; init = 0.0)
    m_tail = sum(m_tail_each; init = 0.0)
    m_tot = m_rep + m_head + m_tail
    nt_bg = sum(df.L[df.class .!= "reporter"]; init = 0)
    audit = codon_demand_audit(mrna, head, tail, conc, reporter, seqs_all)

    if verbose
        println("reduced_cell_load")
        @printf("  source FASTA     = %s  (%d genes with mRNA)\n", fasta_name, nrow(mrna))
        @printf("  n_head / n_tail  = %d / %d    reporter = %s    tail_mode = :%s\n",
                n_head_e, n_tail_e, reporter, tail_mode)
        @printf("  m_reporter       = %.4f µM  (f_reporter = %.3f)\n", m_rep, f_reporter)
        @printf("  m_head           = %.4f µM  (%d genes, xlsx values%s)\n",
                m_head, n_head_e, m_total === nothing ? "" : " rescaled")
        @printf("  m_tail (lumped)  = %.4f µM  (%d tracks carry the leftover proteome)\n",
                m_tail, n_tail_e)
        println("  note: tail µM per gene can exceed head — those tracks are bulk stand-ins,")
        println("        not the 10 lowest genes. Head stays the real highest individual CDSs.")
        @printf("  m_total          = %.4f µM    background CDS length ≈ %d codons\n",
                m_tot, nt_bg)
        println("  change mRNA here: f_reporter, m_total, n_head, n_tail, mRNA_multiplier")
        println("  then pass cell.het_mRNA_conc into call_model(; het_mRNA_conc=...)")
        if startswith(fasta_name, "bg_safe_")
            @warn "source is $fasta_name, not a full geneset — tail mass is only the leftover of these $(nrow(mrna)) genes. Put Ecoli_K12_MG1655_geneset_with_welch_2009.fa in data/sequences to lump the real proteome tail."
        end
        _print_leu_audit(audit)
        if write_fasta
            println("  FASTA            → ", out_path)
        end
        show(df, allrows = true, allcols = true); println()
    end

    return (
        het_genes = het_genes,
        background_genes = background_genes,
        het_mRNA_conc = conc,
        df = df,
        m_total = m_tot,
        m_reporter = m_rep,
        m_head = m_head,
        m_tail = m_tail,
        out_fasta = out_fasta,
        src_fasta = fasta_name,
        reporter = reporter,
        f_reporter = f_reporter,
        n_head = n_head_e,
        n_tail = n_tail_e,
        tail_mode = tail_mode,
        seqs_all = seqs_all,
        codon_audit = audit,
    )
end

function _even_rank_sample(rest::DataFrame, n::Int)
    n <= 0 && return rest[1:0, :]
    nrow(rest) == 0 && return rest
    n >= nrow(rest) && return rest
    idx = unique(round.(Int, range(1, nrow(rest); length = n)))
    while length(idx) < n
        extra = setdiff(1:nrow(rest), idx)
        isempty(extra) && break
        push!(idx, extra[1])
    end
    return rest[sort(idx[1:n]), :]
end

function _default_src_fasta()
    seqdir = sequence_data_dir()
    for name in ("Ecoli_K12_MG1655_geneset_with_welch_2009.fa",
                 "Ecoli_K12_MG1655_geneset.fa")
        isfile(joinpath(seqdir, name)) && return name
    end
    safes = filter(f -> startswith(f, "bg_safe_") && endswith(f, ".fa"), readdir(seqdir))
    isempty(safes) && error("reduced_cell_load: no source FASTA in $seqdir")
    function _nsafe(f)
        m = match(r"bg_safe_(\d+)", f)
        m === nothing ? 0 : parse(Int, m.captures[1])
    end
    return safes[argmax(_nsafe.(safes))]
end

const LEU_DNA = ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")
const STOP_DNA = Set(["TAA", "TAG", "TGA"])

function _codon_dna(c)::String
    return replace(uppercase(string(c)), "U" => "T")
end

function _codon_counts_dna(sd)::Dict{String,Int}
    d = Dict{String,Int}()
    cods = sd.Codon
    n = length(cods)
    n < 2 && return d
    for c in @view cods[1:n-1]
        k = _codon_dna(c)
        length(k) == 3 || continue
        k in STOP_DNA && continue
        d[k] = get(d, k, 0) + 1
    end
    return d
end

function _leu_codon_counts(sd)::Tuple{Int,Int}
    cc = _codon_counts_dna(sd)
    return (get(cc, "CTG", 0), get(cc, "CTA", 0))
end

function _demand_add!(acc::Dict{String,Float64}, sd, mass::Float64)
    for (c, n) in _codon_counts_dna(sd)
        acc[c] = get(acc, c, 0.0) + mass * n
    end
    return acc
end

function _codon_match_tail(rest::DataFrame, n::Int, seqs_all)::DataFrame
    n <= 0 && return rest[1:0, :]
    nrow(rest) == 0 && return rest
    n >= nrow(rest) && return rest
    leftover = Dict{String,Float64}()
    for r in eachrow(rest)
        haskey(seqs_all, r.gene) || continue
        _demand_add!(leftover, seqs_all[r.gene], Float64(r.mRNA))
    end
    s0 = sum(values(leftover); init = 0.0)
    s0 <= 0 && return _even_rank_sample(rest, n)
    target = Dict(k => v / s0 for (k, v) in leftover)
    chosen = String[]
    chosen_set = Set{String}()
    acc = Dict{String,Float64}()
    for _ in 1:n
        best_g = ""
        best_dev = Inf
        for r in eachrow(rest)
            g = String(r.gene)
            g in chosen_set && continue
            haskey(seqs_all, g) || continue
            trial = copy(acc)
            _demand_add!(trial, seqs_all[g], Float64(r.mRNA))
            st = sum(values(trial); init = 0.0)
            st <= 0 && continue
            dev = 0.0
            keys_u = union(keys(target), keys(trial))
            for k in keys_u
                dev += abs(get(trial, k, 0.0) / st - get(target, k, 0.0))
            end
            if dev < best_dev
                best_dev = dev
                best_g = g
            end
        end
        isempty(best_g) && break
        push!(chosen, best_g)
        push!(chosen_set, best_g)
        rw = rest[findfirst(==(best_g), rest.gene), :]
        _demand_add!(acc, seqs_all[best_g], Float64(rw.mRNA))
    end
    idx = [findfirst(==(g), rest.gene) for g in chosen]
    return rest[idx, :]
end

function _leu_frac(demand::Dict{String,Float64})
    leu = sum(get(demand, c, 0.0) for c in LEU_DNA)
    ctg = get(demand, "CTG", 0.0)
    cta = get(demand, "CTA", 0.0)
    return (leu = leu, ctg = ctg, cta = cta,
            ctg_of_leu = leu > 0 ? ctg / leu : NaN,
            cta_of_leu = leu > 0 ? cta / leu : NaN)
end

function codon_demand_audit(mrna, head, tail, conc, reporter, seqs_all)
    proteome = Dict{String,Float64}()
    reduced = Dict{String,Float64}()
    for r in eachrow(mrna)
        haskey(seqs_all, r.gene) || continue
        _demand_add!(proteome, seqs_all[r.gene], Float64(r.mRNA))
    end
    for r in eachrow(vcat(head, tail))
        g = String(r.gene)
        haskey(seqs_all, g) || continue
        _demand_add!(reduced, seqs_all[g], get(conc, g, 0.0))
    end
    sp = sum(values(proteome); init = 0.0)
    sr = sum(values(reduced); init = 0.0)
    l1 = 0.0
    if sp > 0 && sr > 0
        for k in union(keys(proteome), keys(reduced))
            l1 += abs(get(proteome, k, 0.0) / sp - get(reduced, k, 0.0) / sr)
        end
    end
    return (proteome_leu = _leu_frac(proteome),
            reduced_leu = _leu_frac(reduced),
            codon_l1 = l1,
            reporter = reporter)
end

function _print_leu_audit(audit)
    p = audit.proteome_leu
    r = audit.reduced_leu
    println("  expression-weighted Leu demand (background, no reporter):")
    @printf("    proteome  Leu flux = %.3f   CTG/Leu = %.3f   CTA/Leu = %.3f\n",
            p.leu, p.ctg_of_leu, p.cta_of_leu)
    @printf("    reduced   Leu flux = %.3f   CTG/Leu = %.3f   CTA/Leu = %.3f\n",
            r.leu, r.ctg_of_leu, r.cta_of_leu)
    @printf("    codon L1 (reduced vs proteome fractions) = %.4f\n", audit.codon_l1)
    if r.cta_of_leu > 0.05
        @warn "reduced background CTA/Leu is $(round(r.cta_of_leu; digits=3)) — Subramaniam needs a CTG-heavy proteome (CTA rare)."
    end
end

"""
    write_yfp_background_fastas(cell; pnas_fasta=..., out_subdir="reduced_yfp")

Write three FASTAs (CTG / C1A / CTA YFP CDS under `cell.reporter`, same background
genes). Returns Dict label => relative fasta name under data/sequences.
"""
function write_yfp_background_fastas(cell;
        pnas_fasta::Union{Nothing,AbstractString} = nothing,
        out_subdir::AbstractString = "reduced_yfp")
    pnas = pnas_fasta === nothing ? _pnas_sd02_path() : String(pnas_fasta)
    recs = _read_fa_ordered(pnas)
    getrec(p) = last(recs[findfirst(r -> startswith(first(r), p), recs)])
    yfp0, yfp1 = getrec("yfp0"), getrec("yfp1")
    cods0 = [yfp0[i:i+2] for i in 1:3:length(yfp0)]
    cods1 = [yfp1[i:i+2] for i in 1:3:length(yfp1)]
    sens = findall(k -> cods0[k] != cods1[k], eachindex(cods0))
    cods_1 = copy(cods0); cods_1[sens[1]] = "CTA"
    cdss = Dict("CTG" => yfp0, "C1A" => join(cods_1), "CTA" => yfp1)
    src = joinpath(sequence_data_dir(), cell.out_fasta)
    template = isfile(src) ? _read_fa_ordered(src) : _cell_records(cell)
    dest_dir = joinpath(sequence_data_dir(), out_subdir)
    mkpath(dest_dir)
    out = Dict{String,String}()
    for (lab, cds) in cdss
        rel = joinpath(out_subdir, "bg_$(lab).fa")
        _write_fa_replace(template, cell.reporter, cds,
                          joinpath(sequence_data_dir(), rel))
        out[lab] = rel
    end
    @printf("YFP FASTAs (reporter=%s substituted): sites %s\n",
            cell.reporter, string(sens))
    for lab in ("CTG", "C1A", "CTA")
        println("  ", lab, " → ", joinpath(sequence_data_dir(), out[lab]))
    end
    return out
end

function _pnas_sd02_path()
    for p in (joinpath(@__DIR__, "data", "PNAS_2013", "sd02.txt"),
              joinpath(@__DIR__, "..", "data", "PNAS_2013", "sd02.txt"))
        isfile(p) && return p
    end
    error("write_yfp_background_fastas: sd02.txt not found")
end

function _cell_records(cell)
    recs = Pair{String,String}[]
    if haskey(cell.seqs_all, cell.reporter)
        push!(recs, cell.reporter => String(cell.seqs_all[cell.reporter].DNA))
    end
    for g in cell.background_genes
        push!(recs, g => String(cell.seqs_all[g].DNA))
    end
    return recs
end

function _read_fa_ordered(path::AbstractString)
    records = Pair{String,String}[]
    header = ""
    buf = IOBuffer()
    open(path, "r") do io
        for line in eachline(io)
            ln = strip(line)
            isempty(ln) && continue
            if startswith(ln, ">")
                if !isempty(header)
                    push!(records, header => uppercase(String(take!(buf))))
                end
                header = String(strip(ln[2:end]))
                buf = IOBuffer()
            else
                write(buf, ln)
            end
        end
        !isempty(header) && push!(records, header => uppercase(String(take!(buf))))
    end
    return records
end

function _write_fa_replace(template, replace_name, replace_seq, out_path; linewrap::Int = 70)
    mkpath(dirname(out_path))
    found = false
    rseq = uppercase(String(replace_seq))
    open(out_path, "w") do io
        for (name, seq) in template
            if !found && name == replace_name
                println(io, ">", name)
                for i in 1:linewrap:length(rseq)
                    println(io, rseq[i:min(i + linewrap - 1, length(rseq))])
                end
                found = true
            else
                println(io, ">", name)
                for i in 1:linewrap:length(seq)
                    println(io, seq[i:min(i + linewrap - 1, length(seq))])
                end
            end
        end
        if !found
            println(io, ">", replace_name)
            for i in 1:linewrap:length(rseq)
                println(io, rseq[i:min(i + linewrap - 1, length(rseq))])
            end
        end
    end
    return out_path
end

