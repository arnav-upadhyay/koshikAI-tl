# ============================================================================
# make_reduced_background.jl
# ----------------------------------------------------------------------------
# Scientist entry point for a Subramaniam-style mini-transcriptome.
# Does NOT run the ODE. Does NOT touch the production equal-m_bg N=20 cell.
#
# Run from src/:
#   julia --project=.. make_reduced_background.jl
#
# Change concentrations HERE (only):
# ============================================================================

n_head          = 10
n_tail          = 10
reporter        = "rplL"          # YFP CDS is swapped onto this name
f_reporter      = 0.05            # small extra; not the 0.46 production share
m_total         = 2.86            # µM cell-sized pool; `nothing` = native xlsx mass
mRNA_multiplier = 1.0
tail_mode       = :codon_match    # :codon_match | :even_rank
src_fasta       = nothing         # nothing → full geneset if present, else bg_safe_*.fa
write_yfp       = true            # three FASTAs: CTG / C1A / CTA + same background
out_csv         = "../data/output/reduced_cell_load.csv"

# ============================================================================

include("packages.jl")
using CSV

cell = reduced_cell_load(;
    n_head = n_head,
    n_tail = n_tail,
    reporter = reporter,
    f_reporter = f_reporter,
    m_total = m_total,
    mRNA_multiplier = mRNA_multiplier,
    src_fasta = src_fasta,
    tail_mode = tail_mode,
    verbose = true,
)

mkpath(dirname(out_csv))
CSV.write(out_csv, cell.df)
println("Wrote gene table → ", out_csv)

yfp_fastas = write_yfp ? write_yfp_background_fastas(cell) : Dict{String,String}()

println("\ncall_model wiring (freeze off, all genes heterologous):")
println("  fasta = cell.out_fasta          # chromosomal reporter, or yfp_fastas[\"CTG\"] etc.")
println("  genes = cell.het_genes")
println("  conc  = cell.het_mRNA_conc")
println("  call_model(fasta, genes, 0; het_mRNA_conc=conc, freeze_stalled=false, ...)")
if !isempty(yfp_fastas)
    println("\nYFP relative FASTAs (same conc dict, only CDS of $(reporter) changes):")
    for lab in ("CTG", "C1A", "CTA")
        println("  ", lab, "  ", yfp_fastas[lab])
    end
end
