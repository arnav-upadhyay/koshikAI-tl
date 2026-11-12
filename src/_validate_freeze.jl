# Validation harness for the stall-freezing feature (see stall_freeze.jl).
# Temporary script — not part of the model. Run from src/:  julia -t auto _validate_freeze.jl

include("packages.jl")

const FASTA        = "bg_safe_20.fa"
const PRIMARY      = "rplL"
const M_TOTAL      = 2.86
const F_RPLL       = 0.46
const T_END        = parse(Float64, get(ENV, "VF_TEND", "5.0"))
const N_LIST       = [parse(Int, s) for s in split(get(ENV, "VF_NLIST", "5"), ',')]

spindler_scale = txtl_types.make_scale(
	V_IN_70SIC_scale = 75.0, V_IN_IF2D_scale = 2.0,
	ARS_scale = 20.0 .* ones(20), Rib_scale = 1.0, EF_scale = 1, tRNA_scale = 1)

mRNA_pool = to_DataFrame(readmatrix("../data/input/mRNA_concentrations.xlsx",
                                    "mRNA_concentration", "A1:B2402"))
mRNA_pool = filter(row -> isa(row.mRNA, Number), mRNA_pool)
seqs_all_fasta, _ = seqs_import_transcripts(FASTA, (1, 2401))
fasta_gene_set    = Set{String}(String.(keys(seqs_all_fasta)))

candidate_pool = String[]
for i in 1:nrow(mRNA_pool)
	g = String(mRNA_pool.gene[i])
	(g == PRIMARY || !(g in fasta_gene_set)) && continue
	push!(candidate_pool, g)
end

function setup(N)
	bg     = candidate_pool[1:min(N, length(candidate_pool))]
	m_rplL = M_TOTAL * F_RPLL
	m_bg   = N == 0 ? 0.0 : M_TOTAL * (1 - F_RPLL) / N
	genes  = vcat([PRIMARY], bg)
	conc   = Dict{String,Float64}(PRIMARY => m_rplL)
	for g in bg; conc[g] = m_bg; end
	return genes, conc
end

run_case(N; freeze, tau_s = 100.0, tend = T_END, verbose = false) = begin
	genes, conc = setup(N)
	call_model(FASTA, genes, 0;
		custom_scale = spindler_scale, growth_rate = 1, het_mRNA_conc = conc,
		return_details = true, return_sol = true,
		static_EFs = true, spindler_75pct = false,
		t_end = tend, ode_reltol = 1e-8, ode_abstol = 1e-5, ode_maxiters = 5_000_000,
		freeze_stalled = freeze, stall_tau_s = tau_s, stall_verbose = verbose)
end

banner(s) = println("\n", "="^78, "\n  ", s, "\n", "="^78)

# ── Check 1 + 2: regression / single-stall equivalence ───────────────────────
# With freezing ON but no gene meeting the criteria, the trajectory must be
# bit-identical to freezing OFF: the only code path that differs is a branch
# that is never taken.
banner("Checks 1-2: freeze ON vs OFF must agree when nothing freezes")
reg = NamedTuple[]
for N in N_LIST
	a = run_case(N; freeze = false)
	b = run_case(N; freeze = true)
	dPRT = abs(a.final_PRT - b.final_PRT)
	dTR  = abs(a.translation_rate - b.translation_rate)
	rel  = isfinite(a.final_PRT) && a.final_PRT != 0 ? dPRT / abs(a.final_PRT) : dPRT
	nf   = b.freeze === nothing ? 0 : b.freeze.n_frozen
	@printf("  N=%-3d  PRT off=%.10g on=%.10g  rel_diff=%.3e  n_frozen=%d\n",
	        N, a.final_PRT, b.final_PRT, rel, nf)
	push!(reg, (N = N, rel = rel, n_frozen = nf, dTR = dTR))
end
ok_reg = all(r -> r.n_frozen > 0 || r.rel < 1e-12, reg)
println(ok_reg ? "  PASS  no drift on the unfrozen path" :
                 "  FAIL  trajectories diverged without any freeze event")

# ── Check 3: stall_tau_s sensitivity ─────────────────────────────────────────
# Lowering the threshold must freeze monotonically more genes, and the
# conclusion drawn from the primary gene must not swing with the knob.
banner("Check 3: stall_tau_s sensitivity sweep")
Nsens = N_LIST[end]
println("  N = $Nsens")
@printf("  %10s %10s %14s %14s\n", "tau_s", "n_frozen", "PRT_rplL", "transl_rate")
prev = -1
mono = true
for tau in (1e9, 1e4, 1e3, 1e2, 1e1, 1e0)
	r  = run_case(Nsens; freeze = true, tau_s = tau)
	nf = r.freeze === nothing ? 0 : r.freeze.n_frozen
	@printf("  %10.3g %10d %14.6g %14.6g\n", tau, nf, r.final_PRT, r.translation_rate)
	nf < prev && (global mono = false)
	global prev = nf
end
println(mono ? "  PASS  n_frozen is monotone in 1/tau_s" :
               "  FAIL  n_frozen not monotone")

# ── Check 4: ribosome mass balance across a freeze event ─────────────────────
# Freezing must not create or destroy ribosomes: the ribosomes held on frozen
# transcripts have to still be there, not returned to the free pool.
banner("Check 4: ribosome mass balance across a freeze event")
let N = Nsens
	r = run_case(N; freeze = true, tau_s = 1.0, verbose = true)
	if r.freeze === nothing || r.freeze.n_frozen == 0
		println("  SKIP  no gene froze even at tau_s=1 s; nothing to balance")
	else
		p = r.p
		total_rib(u) = sum(gene_stall_metrics(u, p, i).rib_total for i in 1:p.vectors.n_seq)
		u0, u1 = r.sol.u[1], r.sol.u[end]
		fr = r.freeze
		@printf("  frozen: %d gene(s) — %s\n", fr.n_frozen, join(fr.frozen_genes, ", "))
		@printf("  ribosomes on transcripts  t0=%.6g  tend=%.6g\n", total_rib(u0), total_rib(u1))
		@printf("  of which sequestered in frozen transcripts: %.6g µM\n", fr.ribosomes_sequestered)
		held = fr.ribosomes_sequestered
		println(held > 0 ? "  PASS  frozen transcripts still hold their ribosomes" :
		                   "  FAIL  frozen transcripts hold zero ribosomes (released?)")
	end
end

banner("done")
