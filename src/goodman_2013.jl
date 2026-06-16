# ─────────────────────────────────────────────────────────────────────────────
# Goodman, Church, Kosuri 2013  Science 342:475-479  (PMC3909563)
#   "Causes and effects of N-terminal codon bias in bacterial genes"
#
# Library: ~14 000 synonymous variants of the first 11 codons across ~137
# parent E. coli genes, cassetted behind a constitutive promoter (e.g.
# BBaJ23100) + one of a panel of RBSs. mRNA and protein both measured per
# variant ⇒ Trans = Prot / RNA is the per-construct translation efficiency.
#
# This module supports the rplS pilot: load the supplementary table, filter
# to clean rows, splice each variant's 11-codon N-terminus onto the WT rplS
# body (so only the first 11 codons differ across simulations), and produce
# minimal FASTA files for `call_model`.
#
# Plain-script style (no `module`) so a single `include` exposes everything.
# ─────────────────────────────────────────────────────────────────────────────

using Statistics

# ── CSV-safe numeric coercion ────────────────────────────────────────────────
"""
    coerce_float!(df, cols) -> df

Force the named columns to `Float64` (with `missing` for unparseable values).
CSV.jl will fall back to `String` when a numeric column contains `"NA"`,
blanks, or other tokens; this helper rewrites such columns in place so they
can be used with `log10`, `max`, etc.
"""
function coerce_float!(df, cols::Vector{Symbol})
    for c in cols
        hasproperty(df, c) || continue
        col = df[!, c]
        eltype(col) <: Union{Missing, Real} && continue
        new = Vector{Union{Missing, Float64}}(undef, length(col))
        for (i, v) in enumerate(col)
            if ismissing(v) || v === nothing
                new[i] = missing
            elseif v isa Real
                new[i] = Float64(v)
            else
                s = strip(String(v))
                if isempty(s) || uppercase(s) in ("NA", "NAN", "NULL", "#N/A")
                    new[i] = missing
                else
                    parsed = tryparse(Float64, s)
                    new[i] = parsed === nothing ? missing : parsed
                end
            end
        end
        df[!, c] = new
    end
    return df
end

# ── DNA helpers ──────────────────────────────────────────────────────────────
"""
    revcomp_dna(s) -> String

Reverse-complement DNA. Non-ACGTU letters become `N`.
"""
function revcomp_dna(s::AbstractString)
    s = uppercase(s)
    out = IOBuffer()
    for c in Iterators.reverse(s)
        b = c == 'A' ? 'T' :
            c == 'T' ? 'A' :
            c == 'G' ? 'C' :
            c == 'C' ? 'G' :
            c == 'U' ? 'A' : 'N'
        write(out, b)
    end
    return String(take!(out))
end

const _CODON_TABLE = Dict{String,Char}(
    "TTT"=>'F',"TTC"=>'F',"TTA"=>'L',"TTG"=>'L',
    "CTT"=>'L',"CTC"=>'L',"CTA"=>'L',"CTG"=>'L',
    "ATT"=>'I',"ATC"=>'I',"ATA"=>'I',"ATG"=>'M',
    "GTT"=>'V',"GTC"=>'V',"GTA"=>'V',"GTG"=>'V',
    "TCT"=>'S',"TCC"=>'S',"TCA"=>'S',"TCG"=>'S',
    "CCT"=>'P',"CCC"=>'P',"CCA"=>'P',"CCG"=>'P',
    "ACT"=>'T',"ACC"=>'T',"ACA"=>'T',"ACG"=>'T',
    "GCT"=>'A',"GCC"=>'A',"GCA"=>'A',"GCG"=>'A',
    "TAT"=>'Y',"TAC"=>'Y',"TAA"=>'*',"TAG"=>'*',
    "CAT"=>'H',"CAC"=>'H',"CAA"=>'Q',"CAG"=>'Q',
    "AAT"=>'N',"AAC"=>'N',"AAA"=>'K',"AAG"=>'K',
    "GAT"=>'D',"GAC"=>'D',"GAA"=>'E',"GAG"=>'E',
    "TGT"=>'C',"TGC"=>'C',"TGA"=>'*',"TGG"=>'W',
    "CGT"=>'R',"CGC"=>'R',"CGA"=>'R',"CGG"=>'R',
    "AGT"=>'S',"AGC"=>'S',"AGA"=>'R',"AGG"=>'R',
    "GGT"=>'G',"GGC"=>'G',"GGA"=>'G',"GGG"=>'G',
)

"""
    translate_codons(dna) -> String

Naive translation using the standard genetic code. Non-ACGT codons → 'X'.
Stop codons render as '*'.
"""
function translate_codons(dna::AbstractString)
    d = uppercase(dna)
    n = length(d) - mod(length(d), 3)
    out = IOBuffer()
    for i in 1:3:n
        c = d[i:i+2]
        write(out, get(_CODON_TABLE, c, 'X'))
    end
    return String(take!(out))
end

# ── Variant mRNA assembly ────────────────────────────────────────────────────
"""
    build_variant_mrna(variant_cds_head, wt_mrna; n_replace_codons=11) -> String

Glue the variant's first `n_replace_codons` codons onto the wild-type body so
only the N-terminus differs across constructs. Inputs are 5'→3' mRNA-strand
DNA (or RNA-spelled as DNA). Lengths are validated.
"""
function build_variant_mrna(variant_cds_head::AbstractString,
                            wt_mrna::AbstractString;
                            n_replace_codons::Int = 11)
    n_replace_nt = n_replace_codons * 3
    vh = uppercase(variant_cds_head)
    wt = uppercase(wt_mrna)
    length(vh) == n_replace_nt ||
        error("variant_cds_head length $(length(vh)) ≠ expected $(n_replace_nt) nt ($(n_replace_codons) codons)")
    length(wt) > n_replace_nt ||
        error("wt_mrna length $(length(wt)) ≤ replacement region $(n_replace_nt) nt")
    return string(vh, wt[n_replace_nt + 1:end])
end

"""
    write_variant_fasta(path, mrna; gene_name="rplL", linewrap=70)

Write a single-record FASTA. `gene_name` defaults to "rplL" so existing
hard-coded elongation extraction in `call_model` continues to work — the
biological identity of the sequence is conveyed by the parent CSV row, not
the FASTA header.
"""
function write_variant_fasta(path::AbstractString, mrna::AbstractString;
                             gene_name::String = "rplL", linewrap::Int = 70)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, ">", gene_name)
        L = length(mrna)
        for i in 1:linewrap:L
            println(io, mrna[i:min(i + linewrap - 1, L)])
        end
    end
    return path
end

# ── Ordered-FASTA helpers (preserve gene order, needed when call_model uses
#    `seqs_to_use = (1, N)` to select background genes by position) ───────────
"""
    read_fasta_ordered(path) -> Vector{Pair{String,String}}

Parse a FASTA, preserving record order. Uppercase'd sequences. Each entry is
`header => sequence`.
"""
function read_fasta_ordered(path::AbstractString)
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

"""
    write_fasta_with_replacement(template_records, replace_name, replace_seq,
                                 out_path; linewrap=70) -> out_path

Write a copy of `template_records` to `out_path`, but with the entry whose
header equals `replace_name` swapped to use `replace_seq`. If the entry is
not found, the variant is appended at the end. Order is otherwise preserved.
"""
function write_fasta_with_replacement(template_records::Vector{Pair{String,String}},
                                      replace_name::AbstractString,
                                      replace_seq::AbstractString,
                                      out_path::AbstractString;
                                      linewrap::Int = 70)
    mkpath(dirname(out_path))
    found = false
    rseq  = uppercase(String(replace_seq))
    open(out_path, "w") do io
        for (name, seq) in template_records
            if name == replace_name && !found
                _write_fasta_record(io, name, rseq, linewrap)
                found = true
            else
                _write_fasta_record(io, name, seq, linewrap)
            end
        end
        if !found
            _write_fasta_record(io, replace_name, rseq, linewrap)
        end
    end
    return out_path
end

function _write_fasta_record(io::IO, name::AbstractString, seq::AbstractString,
                             linewrap::Int)
    println(io, ">", name)
    L = length(seq)
    for i in 1:linewrap:L
        println(io, seq[i:min(i + linewrap - 1, L)])
    end
    return nothing
end

# ── Variant selection helpers ────────────────────────────────────────────────
"""
    select_variants_quintile(df, n_pick; sort_col=:Trans) -> indices::Vector{Int}

Return `n_pick` row indices into `df` whose `sort_col` values are evenly
spaced through the rank distribution (1st, 2nd, … `n_pick`-th quantile-like
positions). Spans the dynamic range of the column so a small pilot batch
maximises chances of seeing model discrimination.
"""
function select_variants_quintile(df, n_pick::Int; sort_col::Symbol = :Trans)
    N = nrow(df)
    n = min(n_pick, N)
    n < 1 && return Int[]
    perm = sortperm(df[!, sort_col])
    picks = Int[]
    for k in 1:n
        q_idx = round(Int, (k - 0.5) * N / n)
        q_idx = clamp(q_idx, 1, N)
        push!(picks, perm[q_idx])
    end
    return picks
end

# ── Filtering ────────────────────────────────────────────────────────────────
"""
    _truthy(v) -> Bool

Treat both `Bool` (`true`) and `String` ("TRUE", "true") as truthy. `missing`
and everything else become `false`. This is needed because CSV.jl will infer
either a `Bool` or `String` column for the `TRUE`/`FALSE` flags depending on
heuristics.
"""
function _truthy(v)
    v === missing && return false
    v isa Bool && return v
    v isa AbstractString && return uppercase(strip(v)) == "TRUE"
    return false
end

"""
    is_clean_goodman_row(row) -> Bool

Apply the paper's per-row quality flags (`Fltr.SetGood == TRUE`, no insufficient
counts, not out of linear range). Returns false if any of the flags are
missing or fail.
"""
function is_clean_goodman_row(row)
    _truthy(getproperty(row, Symbol("Fltr.SetGood"))) || return false
    for sym in (Symbol("Insuff.Prot"), Symbol("Insuff.DNA"), Symbol("Insuff.RNA"),
                Symbol("Fltr.BelowRange"), Symbol("Fltr.AboveRange"))
        _truthy(getproperty(row, sym)) && return false
    end
    return true
end

# ── Statistics: rank correlation ─────────────────────────────────────────────
"""
    tied_rank_g(v) -> Vector{Float64}

Average-rank for ties (used by Spearman).
"""
function tied_rank_g(v::AbstractVector{<:Real})
    n = length(v)
    p = sortperm(v)
    r = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && v[p[j + 1]] == v[p[i]]
            j += 1
        end
        avg = (i + j) / 2
        for k in i:j
            r[p[k]] = avg
        end
        i = j + 1
    end
    return r
end

"""
    spearman_rho_g(x, y) -> Float64

Spearman rank correlation, robust to NaN/Inf. Returns NaN if <2 valid pairs.
"""
function spearman_rho_g(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || error("x and y must be same length")
    keep = Int[]
    for i in eachindex(x)
        xi = float(x[i]); yi = float(y[i])
        if isfinite(xi) && isfinite(yi)
            push!(keep, i)
        end
    end
    length(keep) < 2 && return NaN
    return cor(tied_rank_g(Float64.(x[keep])), tied_rank_g(Float64.(y[keep])))
end

"""
    countmap_to_str(v) -> String

Pretty short summary of categorical value counts, e.g. "WT=42, Strong=10".
Skips `missing`. Sorted by descending count.
"""
function countmap_to_str(v)
    d = Dict{Any,Int}()
    for x in v
        ismissing(x) && continue
        d[x] = get(d, x, 0) + 1
    end
    pairs_sorted = sort(collect(d); by = kv -> -kv[2])
    return join(["$(k)=$(c)" for (k, c) in pairs_sorted], ", ")
end
