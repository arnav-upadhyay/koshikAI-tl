# ─────────────────────────────────────────────────────────────────────────────
# Mohammad / Buskirk 2019 ribo-seq utilities  (eLife 8:e42591, GSE119104)
#
# Reads the fixedStep .wig files produced by the lab's bacterial Ribo-seq
# pipeline (github.com/greenlabjhmi/2018_Bacterial_Pipeline_riboseq), slices
# them at gene coordinates, applies the 3'-end → A-site offset, and converts
# per-nucleotide footprint density into per-codon dwell time so it can be
# compared head-to-head with the model's τ profile in elongation_rate_log.csv.
#
# Plain-script style (no `module`) so it can be `include`-d from main.ipynb
# without having to fully qualify every call.
# ─────────────────────────────────────────────────────────────────────────────

using Statistics

# ── WIG parser ───────────────────────────────────────────────────────────────
"""
    parse_mohammad_wig(path) -> (chrom::String, density::Vector{Float64})

Parse a fixedStep WIG produced by Mohammad's `ribo_util.countstowig`:

    track type=wiggle_0 ...
    fixedStep  chrom="..."  start=1  step=1
    <value at pos 1>
    <value at pos 2>
    ...

`density[p]` gives the read count (or rpm) at 1-based chromosomal position `p`.
Lines that are not a numeric value are treated as headers and skipped.
"""
function parse_mohammad_wig(path::AbstractString)
    chrom = ""
    density = Float64[]
    sizehint!(density, 5_000_000)
    open(path, "r") do io
        for line in eachline(io)
            ln = strip(line)
            isempty(ln) && continue
            c1 = ln[1]
            if c1 == 't' || c1 == 'f' || c1 == 'v' || c1 == '#'
                m = match(r"chrom=\"?([^\s\"]+)\"?"i, ln)
                m !== nothing && (chrom = String(m.captures[1]))
                continue
            end
            push!(density, parse(Float64, ln))
        end
    end
    return chrom, density
end

# ── FASTA helpers ────────────────────────────────────────────────────────────
"""
    read_fasta(path) -> Dict{String, String}

Tiny FASTA loader. Returns `Dict{header_without_>, uppercase_sequence}`.
"""
function read_fasta(path::AbstractString)
    records = Dict{String, String}()
    header = ""
    buf = IOBuffer()
    open(path, "r") do io
        for line in eachline(io)
            ln = strip(line)
            isempty(ln) && continue
            if startswith(ln, ">")
                if !isempty(header)
                    records[header] = uppercase(String(take!(buf)))
                end
                header = String(strip(ln[2:end]))
                buf = IOBuffer()
            else
                write(buf, ln)
            end
        end
        !isempty(header) && (records[header] = uppercase(String(take!(buf))))
    end
    return records
end

"""
    revcomp(s) -> String

Reverse-complement DNA. Non-ACGTU letters become `N`.
"""
function revcomp(s::AbstractString)
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

# ── Gene slicing ─────────────────────────────────────────────────────────────
"""
    extract_gene_sequence(genome, start, stop, strand) -> String

`genome[start:stop]` (1-based inclusive). If `strand == :minus`, reverse-
complements the slice so the returned sequence runs 5'→3' along the mRNA.
"""
function extract_gene_sequence(genome::AbstractString,
                               gene_start::Int, gene_stop::Int,
                               strand::Symbol)
    @assert 1 <= gene_start <= gene_stop <= length(genome) "coords out of range"
    plus = genome[gene_start:gene_stop]
    return strand == :minus ? revcomp(plus) : String(plus)
end

"""
    extract_gene_density(wig, start, stop, strand) -> Vector{Float64}

`wig[start:stop]` reversed if `strand == :minus`, so the returned vector runs
5'→3' along the mRNA (parallel to `extract_gene_sequence`).
"""
function extract_gene_density(wig::Vector{Float64},
                              gene_start::Int, gene_stop::Int,
                              strand::Symbol)
    @assert 1 <= gene_start <= gene_stop <= length(wig) "coords out of range"
    seg = wig[gene_start:gene_stop]
    return strand == :minus ? reverse(seg) : seg
end

# ── Per-codon dwell-time conversion ──────────────────────────────────────────
"""
    mrna_density_to_codon_dwell(mrna_density;
                                a_site_shift = 14,
                                start_trim_codons = 20,
                                stop_trim_codons = 20)

Convert per-nucleotide 3'-end ribosome footprint density along an mRNA (5'→3')
into per-codon dwell time.

For each codon `i` (1-indexed, mRNA positions `3i-2, 3i-1, 3i`), the
A-site footprint reads are summed at mRNA positions `3i-2+s, 3i-1+s, 3i+s`
with `s = a_site_shift`. This matches Mohammad's `A_site shift = 14`
convention for 23–25 nt footprints (3' end of footprint lies ~14 nt
downstream of the A-site codon on the mRNA).

`relative_dwell[i] = codon_reads[i] / mean(codon_reads on body)`, where the
body excludes the first `start_trim_codons` (initiation pileup) and last
`stop_trim_codons` (termination pileup) codons.

Returns `(codon_reads, relative_dwell, body_idx)`. Codons whose A-site window
falls past the mRNA end get `NaN`.
"""
function mrna_density_to_codon_dwell(mrna_density::Vector{Float64};
                                     a_site_shift::Int = 14,
                                     start_trim_codons::Int = 20,
                                     stop_trim_codons::Int = 20)
    n_nt = length(mrna_density)
    n_codons = div(n_nt, 3)
    codon_reads = fill(NaN, n_codons)
    for i in 1:n_codons
        p1 = 3*(i - 1) + 1 + a_site_shift
        p3 = 3*i           + a_site_shift
        (p1 < 1 || p3 > n_nt) && continue
        codon_reads[i] = sum(@view mrna_density[p1:p3])
    end

    lo = start_trim_codons + 1
    hi = n_codons - stop_trim_codons
    body_idx = Int[i for i in lo:hi if !isnan(codon_reads[i])]
    rel_dwell = fill(NaN, n_codons)
    if !isempty(body_idx)
        mean_body = mean(codon_reads[body_idx])
        if mean_body > 0
            for i in 1:n_codons
                isnan(codon_reads[i]) || (rel_dwell[i] = codon_reads[i] / mean_body)
            end
        end
    end
    return codon_reads, rel_dwell, body_idx
end

# ── Statistics: rank correlation + top-K overlap ─────────────────────────────
"""
    tied_rank(v) -> Vector{Float64}

Average-rank for ties (used by Spearman).
"""
function tied_rank(v::AbstractVector{<:Real})
    n = length(v)
    p = sortperm(v)
    r = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && v[p[j+1]] == v[p[i]]
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
    spearman_rho(x, y) -> Float64

Spearman rank correlation, ignoring index positions where either value is NaN.
Returns NaN if fewer than 2 valid pairs remain.
"""
function spearman_rho(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y) "x and y must be same length"
    keep = [i for i in eachindex(x) if !isnan(x[i]) && !isnan(y[i])]
    length(keep) < 2 && return NaN
    return cor(tied_rank(x[keep]), tied_rank(y[keep]))
end

"""
    jaccard_topk(x, y, k; idx_pool=eachindex(x))
        -> (jaccard, top_x_indices, top_y_indices)

Top-k overlap (Jaccard) between two ranked profiles. Restricts to `idx_pool`,
ignores NaNs, picks the `k` largest values from each, returns `|A∩B| / |A∪B|`.
"""
function jaccard_topk(x::AbstractVector, y::AbstractVector, k::Int;
                      idx_pool = eachindex(x))
    keep = [i for i in idx_pool if !isnan(x[i]) && !isnan(y[i])]
    isempty(keep) && return (NaN, Int[], Int[])
    kk = min(k, length(keep))
    top_x = sort(keep; by = i -> -x[i])[1:kk]
    top_y = sort(keep; by = i -> -y[i])[1:kk]
    inter = intersect(Set(top_x), Set(top_y))
    uni   = union(Set(top_x), Set(top_y))
    return (length(inter) / length(uni), top_x, top_y)
end

# ── Assembly-version offset helper ───────────────────────────────────────────
"""
    pick_best_wig_offset(wig, base_start, base_stop, strand;
                        candidates = [0, -776, -1977, +776, +1977])

The user's WIG file is aligned to NC_000913.2, but their FASTA is .3 (v3 is
~1977 bp longer due to small insertions, mostly the CRISPR expansion at
~3.65 Mb plus a few smaller indels). For a v3 coordinate downstream of those
insertions, the equivalent v2 coordinate is `v3_coord - 1977` (or thereabouts).

Rather than make the user look this up, we try a few candidate shifts of the
given coordinates and pick the one that yields the most total reads in the
extracted region. Returns `(best_offset, totals::Vector{Tuple{Int,Float64}})`.
"""
function pick_best_wig_offset(wig::Vector{Float64},
                              base_start::Int, base_stop::Int, strand::Symbol;
                              candidates::Vector{Int} = [0, -776, -1977, +776, +1977])
    L = length(wig)
    totals = Tuple{Int, Float64}[]
    for off in candidates
        s = base_start + off
        e = base_stop  + off
        if 1 <= s <= e <= L
            push!(totals, (off, sum(view(wig, s:e))))
        else
            push!(totals, (off, -Inf))
        end
    end
    best = totals[argmax(map(t -> t[2], totals))]
    return best[1], totals
end
