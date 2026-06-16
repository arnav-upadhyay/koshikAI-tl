# Simple RHS prints for ODE_parallel_distributed! — largest |du| components (convergence / stiffness hints).
#
# Manual: print_problematic_rhs(sol.u[end], sol.prob.p, sol.t[end]; top_k=20)

"""
	print_convergence_root_hint(u, p)

Uses the same EF-Tu / EFG formulas as `rhs_metabolites.jl` (metabolite block only) and prints
a short interpretation: negative EFT/EGD, tiny Michaelis denominator for d(EFT)/dt, etc.
Called automatically from `print_problematic_rhs` (and thus from `call_model` on solver failure).
"""
function print_convergence_root_hint(u, p)
	tags = p.local_tags_list[1]
	idx = tags.idx
	n_metab = p.vectors.n_metab
	@views um = u[1:n_metab]
	EFG = um[idx.EFG]
	EFD = um[idx.EFD]
	EGT = um[idx.EGT]
	EGD = um[idx.EGD]
	EFT = um[idx.EFT]

	@unpack k_EFGGDP_Diss, k_EFG_GDP_Ass, k_EFG_GTP_Ass, k_EFGGTP_Diss,
		k_EFTsf, k_EFTsr, Keq_EFTu, KM_GTP, KM_GDP, KM_EFTuGDP, KM_EFTuGTP, Ki_EFTuGDP, Ki_EFTuGTP = tags.p_met.EF
	cv = tags.init_c.constant_vars
	GTP = cv[2]
	GDP = cv[6]
	Est = cv[14]

	println()
	println("========== Convergence root hint (elongation-factor pool) ==========")
	println("  EFG=", EFG, "  EFD=", EFD, "  EGT=", EGT, "  EGD=", EGD, "  EFT=", EFT)
	println("  (names as in rhs_metabolites.jl / Tag; du[EFD] comment = EF-Tu·GDP mass balance)")

	neg = String[]
	for (nm, v) in (("EFG", EFG), ("EFD", EFD), ("EGT", EGT), ("EGD", EGD), ("EFT", EFT))
		v < 0 && push!(neg, nm)
	end
	if !isempty(neg)
		println("  >>> ISSUE: negative pool(s): ", join(neg, ", "), " — concentrations should stay ≥ 0.")
		println("      EFT<0: d(EFT)/dt uses a Michaelis-like denom with several terms ∝ EFT (rhs_metabolites.jl);")
		println("              negative EFT shrinks the denominator → huge synthesis contribution to du[EFT].")
		println("      EGD<0: d(EGD) includes -k_Diss*EGD → large positive spike in du[EGD].")
	else
		println("  EF pools are nonnegative at this u (if failure persists, check tiny denominator below).")
	end

	r_egt = k_EFG_GTP_Ass * EFG * GTP - k_EFGGTP_Diss * EGT
	r_egd = -(k_EFGGDP_Diss * EGD - k_EFG_GDP_Ass * EFG * GDP)
	println("  Metabolite-only piece of d(EGT) ≈ ", r_egt)
	println("  Metabolite-only piece of d(EGD) ≈ ", r_egd)

	if abs(Est) < 1e-24
		println("  (skip EF-Tu rational flux: Est ≈ 0)")
	else
		numer = (k_EFTsf * Est) * (EFD * GTP - (GDP * EFT / Keq_EFTu))
		den = (KM_GTP * EFD + KM_EFTuGDP * GTP
			+ (k_EFTsf * Est * KM_EFTuGTP * GDP) / (k_EFTsr * Est * Keq_EFTu)
			+ (k_EFTsf * Est * KM_GDP * EFT) / (k_EFTsr * Est * Keq_EFTu) + EFD * GTP
			+ (k_EFTsf * Est * KM_EFTuGTP * EFD * GDP) / (k_EFTsr * Est * Keq_EFTu * Ki_EFTuGDP)
			+ (k_EFTsf * Est * GDP * EFT) / (k_EFTsr * Est * Keq_EFTu)
			+ (KM_EFTuGDP * GTP * EFT) / Ki_EFTuGTP)
		println("  EF-Tu synthesis: numer=", numer, "  denom=", den)
		if !isfinite(den) || abs(den) < 1e-6
			println("  >>> ISSUE: denominator tiny or non-finite → Michaelis-like d(EFT) term blows up.")
		elseif den < 0
			println("  >>> ISSUE: negative denominator → sign conflict in d(EFT)/dt rational form.")
		else
			println("  synthesis numer/denom ≈ ", numer / den, "  (before subtracting T3 consumption in full RHS)")
		end
	end

	println("  Note: rhs_sequence.jl also adds elongation/recycling terms to d(EGT), d(EGD); full du sums threads.")
	println("========== end convergence hint ==========")
	println()
	return nothing
end

function build_ode_state_labels(p)
	@unpack y_codons_lengths, offsets, n_metab, n_seq = p.vectors
	lt = p.local_tags_list
	tag_full = lt[1].species_and_sequences.Tag
	(length(tag_full) - 4 == n_metab) ||
		error("build_ode_state_labels: n_metab mismatch with Tag length")
	names_metab = String.(tag_full[1:end-4])
	ntot = n_metab + sum(y_codons_lengths)
	labels = Vector{String}(undef, ntot)
	labels[1:n_metab] .= names_metab
	for i in 1:n_seq
		gene = lt[i].sequence_name
		base = n_metab + offsets[i]
		len = y_codons_lengths[i]
		labels[base + 1] = gene * "|S70"
		labels[base + 2] = gene * "|rna"
		labels[base + 3] = gene * "|PRT"
		labels[base + 4] = gene * "|PTC"
		for j in 1:(len - 4)
			labels[base + 4 + j] = gene * "|codon_" * string(j)
		end
	end
	return labels
end

function _seq_block_sublabel(j_within_block::Int)
	j_within_block == 1 && return "S70"
	j_within_block == 2 && return "rna"
	j_within_block == 3 && return "PRT"
	j_within_block == 4 && return "PTC"
	return "codon_" * string(j_within_block - 4)
end

function print_problematic_rhs(u, p, t; top_k::Int = 20)
	du = similar(u)
	fill!(du, 0.0)
	ODE_parallel_distributed!(du, u, p, t)
	labels = build_ode_state_labels(p)

	nfu = 0
	nfd = 0
	@inbounds for i in eachindex(u)
		isfinite(u[i]) || (nfu += 1)
	end
	@inbounds for i in eachindex(du)
		isfinite(du[i]) || (nfd += 1)
	end

	println()
	println("========== Problematic RHS (largest |du|)  t = ", t, " ==========")
	println("nonfinite count: u = ", nfu, "   du = ", nfd)

	absdu = abs.(du)
	perm = sortperm(absdu; rev = true)
	k = min(top_k, length(u))
	println("top ", k, " components by |du|:")
	for r in 1:k
		i = perm[r]
		println("  ", r, "  idx=", i, "  |du|=", absdu[i], "  du=", du[i], "  u=", u[i], "  :: ", labels[i])
	end

	@unpack y_codons_lengths, offsets, n_metab, n_seq = p.vectors
	println("per-gene max |du|:")
	for si in 1:n_seq
		base = n_metab + offsets[si]
		len = y_codons_lengths[si]
		mx = 0.0
		jm = 1
		@inbounds for j in 1:len
			v = abs(du[base + j])
			if v > mx
				mx = v
				jm = j
			end
		end
		g = p.local_tags_list[si].sequence_name
		println("  ", g, "  max|du|=", mx, "  at ", _seq_block_sublabel(jm))
	end

	print_convergence_root_hint(u, p)

	println("========== end RHS report ==========")
	println()
	return nothing
end

"""
	ode_bottleneck_report(res, t_sample; primary_gene="rplL", top_k=15,
	                       n_codon_show=8, n_t3_show=8)

Post-solve "where is the model stiff and where is rplL stalling?" report. Designed
to be called once after a successful (or partial) `call_model(..., return_details=true,
return_sol=true)` so we can answer:

1. **Solver-internal**: how many steps / rejects / Jacobians did FBDF need? What was
   the smallest dt the solver had to take, and where in time? — these flag stiffness.
2. **State-level**: which ODE component(s) have the largest |du| or the largest
   relative rate |du| / max(|u|, ε) at `t_sample`? Names are decoded via
   `build_ode_state_labels` so per-gene codon positions show up by gene name and
   codon index. The relative-rate ranking exposes *stiff* variables (small u, large
   |du|) that the raw |du| ranking sometimes hides behind the big-pool species.
3. **Biology**: which codon positions on `primary_gene` are the slowest, and which
   charged-tRNA (T3) species are most depleted at `t_sample`? These are the
   actual translation-elongation bottlenecks — typically the cause of the
   `T3_min ≈ 0 → codon stall → Σ 1/v_i blows up` pattern we documented in the
   2026-05-19 CHANGELOG entry.

`res` must have fields `.sol` (the ODESolution returned by DifferentialEquations.jl)
and `.p` (the parameter struct passed to ODE_parallel_distributed!).
"""
function ode_bottleneck_report(res, t_sample::Float64;
		primary_gene::String = "rplL",
		top_k::Int = 15,
		n_codon_show::Int = 8,
		n_t3_show::Int = 8)

	sol = res.sol
	p   = res.p
	t_e = clamp(t_sample, sol.t[1], sol.t[end])

	println()
	println("════════════════════════════════════════════════════════════════════════════")
	println("  ODE bottleneck report  —  t_sample=", t_sample, " min  (clamped to t_e=", t_e, ")")
	println("════════════════════════════════════════════════════════════════════════════")

	# ─── 1. Solver statistics (FBDF accept/reject/jacobian counts) ──────────────
	stats = nothing
	if hasproperty(sol, :destats)
		stats = sol.destats
	elseif hasproperty(sol, :stats)
		stats = sol.stats
	end
	if stats !== nothing
		na = hasproperty(stats, :naccept) ? stats.naccept : (hasproperty(stats, :accept) ? stats.accept : missing)
		nr = hasproperty(stats, :nreject) ? stats.nreject : (hasproperty(stats, :reject) ? stats.reject : missing)
		nj = hasproperty(stats, :njacs)   ? stats.njacs   : (hasproperty(stats, :nw)      ? stats.nw      : missing)
		nf = hasproperty(stats, :nf)      ? stats.nf      : missing
		println("  -- Solver statistics --")
		println("     naccept = ", na, "   nreject = ", nr,
			"   njacs = ", nj, "   nf = ", nf)
		if !ismissing(na) && !ismissing(nr) && (na + nr) > 0
			@printf("     reject rate = %.2f%%  (stiff/oscillatory if > ~25%%)\n",
				100.0 * nr / (na + nr))
		end
	else
		println("  -- Solver statistics unavailable on sol --")
	end

	# ─── 2. dt history (where did the solver have to step tiny?) ───────────────
	ts = sol.t
	dts = diff(ts)
	if !isempty(dts)
		i_min = argmin(dts)
		@printf("  -- dt history --\n")
		@printf("     smallest dt = %.3e min  at t = %.4f min  (step %d / %d)\n",
			dts[i_min], ts[i_min], i_min, length(ts))
		@printf("     median   dt = %.3e min   |   max dt = %.3e min   |   n_steps = %d\n",
			median(dts), maximum(dts), length(dts))
	end

	# ─── 3. RHS at t_e — top components by |du| and by |du|/|u| ────────────────
	u = sol(t_e)
	du = similar(u)
	fill!(du, 0.0)
	try
		ODE_parallel_distributed!(du, u, p, t_e)
	catch e
		println("  -- RHS evaluation at t_e threw: ", e, "  (skipping state-level ranking)")
		println("════════════════════════════════════════════════════════════════════════════")
		return nothing
	end

	labels = build_ode_state_labels(p)
	absdu = abs.(du)
	perm_abs = sortperm(absdu; rev = true)

	println("  -- Top ", top_k, " states by |du(t_e)| --")
	@printf("     %-3s  %-30s  %12s  %14s\n", "rk", "state", "|du|", "u")
	for r in 1:min(top_k, length(absdu))
		i = perm_abs[r]
		@printf("     %-3d  %-30s  %.3e   %.4e\n", r, labels[i], absdu[i], u[i])
	end

	# Relative rate |du|/max(|u|, eps) — surfaces stiff small-pool variables that
	# raw |du| ranking misses behind big-pool species like mRNA / S30 / S50.
	ratio = similar(absdu)
	floor_u = 1e-12
	@inbounds for i in eachindex(absdu)
		ratio[i] = absdu[i] / max(abs(u[i]), floor_u)
	end
	perm_rel = sortperm(ratio; rev = true)
	println("\n  -- Top ", top_k, " states by |du|/max(|u|, 1e-12)  (relative-rate / stiffness proxy) --")
	@printf("     %-3s  %-30s  %14s  %12s  %14s\n", "rk", "state", "|du|/|u|", "|du|", "u")
	for r in 1:min(top_k, length(ratio))
		i = perm_rel[r]
		@printf("     %-3d  %-30s  %.3e      %.3e   %.4e\n",
			r, labels[i], ratio[i], absdu[i], u[i])
	end

	# ─── 4. Slowest codons on primary gene + T3_min argmin ─────────────────────
	tags = p.local_tags_list[1]
	idx  = tags.idx
	t3_vec    = u[idx.a1t : idx.y1t]
	t3_labels = String.(strip.(tags.species_and_sequences.Tag[idx.a1t : idx.y1t]))
	KM_T3j    = tags.p_seq.EL.KM_T3j

	perm_t3 = sortperm(t3_vec)
	println("\n  -- Most-depleted T3 (charged tRNA·EF-Tu·GTP) pools at t_e --")
	@printf("     (KM_T3j on primary block = %.4f µM — pools below this are saturating downward)\n", KM_T3j)
	@printf("     %-3s  %-12s  %14s  %14s\n", "rk", "species", "u (µM)", "KM_T3j / u")
	for r in 1:min(n_t3_show, length(t3_vec))
		i = perm_t3[r]
		ratio_km = t3_vec[i] > 0 ? KM_T3j / t3_vec[i] : Inf
		@printf("     %-3d  %-12s  %.3e        %.3e\n", r, t3_labels[i], t3_vec[i], ratio_km)
	end

	# Slowest codons on primary gene — use spindler_codon_translation_time which
	# already gives smoothed per-codon v + τ.
	try
		sp = spindler_codon_translation_time(sol, p, primary_gene; t_sample = t_e)
		tau = sp.tau
		v   = sp.v_smooth
		perm_slow = sortperm(tau; rev = true)
		println("\n  -- Slowest codons on ", primary_gene, "  (τ = 1 / v_smooth) --")
		@printf("     L = %d codons   |   Σ τ = %.3f s   |   mean v = %.3f codons/s\n",
			sp.L, sp.total_time_s, sp.mean_v)
		@printf("     %-3s  %-6s  %14s  %14s\n", "rk", "pos", "τ (s)", "v (codons/s)")
		for r in 1:min(n_codon_show, length(tau))
			i = perm_slow[r]
			@printf("     %-3d  %-6d  %.4e        %.4e\n", r, i, tau[i], v[i])
		end
		# Decode codon identity at the slowest position into RNA letters → tells the
		# user *which* codon (and therefore which T3 species) is the local bottleneck.
		try
			si = _find_seq_index(p, primary_gene)
			tags_p = p.local_tags_list[si]
			i_top = perm_slow[1]
			cod_rna = string(tags_p.species_and_sequences.Seq[i_top + 1])
			println("     slowest codon at pos ", i_top, " reads ", cod_rna,
				"  (decode T3 letter via Tag[a1t..y1t] above)")
		catch
			# Non-fatal; the codon-position diagnostic is the actionable part.
		end
	catch e
		println("\n  -- spindler_codon_translation_time threw at t_e=", t_e, ": ", e)
	end

	println("════════════════════════════════════════════════════════════════════════════")
	return nothing
end

"""
	print_timeout_snapshot_diagnostic(u_snap, t_snap, p; primary_gene="rplL",
	                                   top_k=15, n_codon_show=8, n_t3_show=8)

Tier-2 partner of `ode_bottleneck_report` for the case where `_run_with_timeout`
killed the wait on a slow integration and we therefore do NOT have a full
`sol` object — only the last `(u, t, p)` snapshot that `call_model`'s
`PeriodicCallback` wrote into the caller's `progress_ref`.

Produces the same actionable rankings as `ode_bottleneck_report` minus the
solver-stats / dt-history sections that need `sol.destats` and `sol.t`:
1. Top-k state components by raw `|du(snap)|`, labelled via `build_ode_state_labels`.
2. Top-k state components by `|du| / max(|u|, 1e-12)` — stiffness proxy that
   surfaces small-pool species (T3, codon positions) hidden by big-pool RHS.
3. Most-depleted T3 (charged tRNA·EF-Tu·GTP) species with `KM_T3j / u`.
4. Slowest codons on `primary_gene` via `spindler_codon_translation_time_at_u`.

Returns a NamedTuple with the same key rankings (plus the snapshot t) so the
caller can persist them to a CSV without rerunning the analysis.
"""
function print_timeout_snapshot_diagnostic(u_snap, t_snap::Float64, p;
		primary_gene::String = "rplL",
		top_k::Int = 15,
		n_codon_show::Int = 8,
		n_t3_show::Int = 8)

	println()
	println("════════════════════════════════════════════════════════════════════════════")
	println("  TIMEOUT snapshot diagnostic  —  last captured t = ", t_snap, " min")
	println("  (integration did NOT reach t_sample; using PeriodicCallback snapshot)")
	println("════════════════════════════════════════════════════════════════════════════")

	# 1. RHS at the snapshot
	du = similar(u_snap)
	fill!(du, 0.0)
	try
		ODE_parallel_distributed!(du, u_snap, p, t_snap)
	catch e
		println("  -- RHS evaluation at snapshot threw: ", e, "  (returning empty result)")
		println("════════════════════════════════════════════════════════════════════════════")
		return nothing
	end

	labels = build_ode_state_labels(p)
	absdu  = abs.(du)
	perm_abs = sortperm(absdu; rev = true)

	println("  -- Top ", top_k, " states by |du(snap)| --")
	@printf("     %-3s  %-30s  %12s  %14s\n", "rk", "state", "|du|", "u")
	for r in 1:min(top_k, length(absdu))
		i = perm_abs[r]
		@printf("     %-3d  %-30s  %.3e   %.4e\n", r, labels[i], absdu[i], u_snap[i])
	end

	# 2. Relative-rate stiffness proxy
	ratio = similar(absdu)
	@inbounds for i in eachindex(absdu)
		ratio[i] = absdu[i] / max(abs(u_snap[i]), 1e-12)
	end
	perm_rel = sortperm(ratio; rev = true)
	println("\n  -- Top ", top_k, " states by |du|/max(|u|, 1e-12)  (stiffness proxy) --")
	@printf("     %-3s  %-30s  %14s  %12s  %14s\n", "rk", "state", "|du|/|u|", "|du|", "u")
	for r in 1:min(top_k, length(ratio))
		i = perm_rel[r]
		@printf("     %-3d  %-30s  %.3e      %.3e   %.4e\n",
			r, labels[i], ratio[i], absdu[i], u_snap[i])
	end

	# 3. T3 depletion
	tags = p.local_tags_list[1]
	idx  = tags.idx
	t3_vec    = u_snap[idx.a1t : idx.y1t]
	t3_labels = String.(strip.(tags.species_and_sequences.Tag[idx.a1t : idx.y1t]))
	KM_T3j    = tags.p_seq.EL.KM_T3j
	perm_t3   = sortperm(t3_vec)

	println("\n  -- Most-depleted T3 (charged tRNA·EF-Tu·GTP) pools at snapshot --")
	@printf("     (KM_T3j on primary block = %.4f µM)\n", KM_T3j)
	@printf("     %-3s  %-12s  %14s  %14s\n", "rk", "species", "u (µM)", "KM_T3j / u")
	for r in 1:min(n_t3_show, length(t3_vec))
		i = perm_t3[r]
		ratio_km = t3_vec[i] > 0 ? KM_T3j / t3_vec[i] : Inf
		@printf("     %-3d  %-12s  %.3e        %.3e\n", r, t3_labels[i], t3_vec[i], ratio_km)
	end

	# 4. Slowest codons on primary gene — uses u-only variant so no sol needed.
	slow_pos = Int[]
	slow_tau = Float64[]
	slow_v   = Float64[]
	try
		sp = spindler_codon_translation_time_at_u(u_snap, p, primary_gene, t_snap)
		tau = sp.tau
		v   = sp.v_smooth
		perm_slow = sortperm(tau; rev = true)
		println("\n  -- Slowest codons on ", primary_gene, "  (snapshot, τ = 1 / v_smooth) --")
		@printf("     L = %d codons   |   Σ τ = %.3f s   |   mean v = %.3f codons/s\n",
			sp.L, sp.total_time_s, sp.mean_v)
		@printf("     %-3s  %-6s  %14s  %14s\n", "rk", "pos", "τ (s)", "v (codons/s)")
		for r in 1:min(n_codon_show, length(tau))
			i = perm_slow[r]
			@printf("     %-3d  %-6d  %.4e        %.4e\n", r, i, tau[i], v[i])
			push!(slow_pos, i); push!(slow_tau, tau[i]); push!(slow_v, v[i])
		end
		try
			si = _find_seq_index(p, primary_gene)
			tags_p = p.local_tags_list[si]
			i_top = perm_slow[1]
			cod_rna = string(tags_p.species_and_sequences.Seq[i_top + 1])
			println("     slowest codon at pos ", i_top, " reads ", cod_rna,
				"  (decode T3 letter via Tag[a1t..y1t] above)")
		catch
		end
	catch e
		println("\n  -- spindler_codon_translation_time_at_u threw at snapshot: ", e)
	end

	println("════════════════════════════════════════════════════════════════════════════")

	return (
		t_snap            = t_snap,
		top_du_idx        = perm_abs[1:min(top_k, length(absdu))],
		top_du_labels     = labels[perm_abs[1:min(top_k, length(absdu))]],
		top_du_vals       = absdu[perm_abs[1:min(top_k, length(absdu))]],
		top_du_u          = u_snap[perm_abs[1:min(top_k, length(absdu))]],
		top_ratio_idx     = perm_rel[1:min(top_k, length(ratio))],
		top_ratio_labels  = labels[perm_rel[1:min(top_k, length(ratio))]],
		top_ratio_vals    = ratio[perm_rel[1:min(top_k, length(ratio))]],
		top_ratio_du      = absdu[perm_rel[1:min(top_k, length(ratio))]],
		top_ratio_u       = u_snap[perm_rel[1:min(top_k, length(ratio))]],
		t3_labels         = t3_labels[perm_t3[1:min(n_t3_show, length(t3_vec))]],
		t3_vals           = t3_vec[perm_t3[1:min(n_t3_show, length(t3_vec))]],
		slow_codon_pos    = slow_pos,
		slow_codon_tau    = slow_tau,
		slow_codon_v      = slow_v,
	)
end

"""
	print_ef_pool_literature_audit(y_state, initial_concentrations, p_met, idx)

After `EF_balance!`, print `Eut`, `Est`, `S50`, `c_EFG_t`, and simple ratios for comparison
to fast-growing *E. coli* (Tu/Rib ~5–9×, Ts/Rib ~1×, Tu/Ts ~5–9×, EF-G/Rib ~1.5–2.5×).
"""
function print_ef_pool_literature_audit(y_state, initial_concentrations, p_met, idx)
	cv = initial_concentrations.constant_vars
	Eut = cv[13]
	Est = cv[14]
	S50 = Float64(y_state[idx.S50])
	cG = p_met.EF.c_EFG_t
	s_efg = Float64(y_state[idx.EFG] + y_state[idx.EGT] + y_state[idx.EGD])

	println()
	println("========== EF pool audit (literature bands: Tu/Rib ~5–9×, Ts/Rib ~1×, Tu/Ts ~5–9×, EF-G/Rib ~1.5–2.5×) ==========")
	println("  Eut (total EF-Tu, constant_vars[13])     = ", Eut, " µM")
	println("  Est (EF-Ts, constant_vars[14])          = ", Est, " µM")
	println("  S50 (50S subunit, state after balance)  = ", S50, " µM")
	println("  c_EFG_t (DAE total EF-G, parameters.jl) = ", cG, " µM")
	println("  EFG+EGT+EGD (state after EF_balance!)   = ", s_efg, " µM  (should ≈ c_EFG_t)")
	if S50 > 1e-12
		println("  Eut/S50      = ", Eut / S50, "   [~5–9]")
		println("  Est/S50      = ", Est / S50, "   [~~1]")
		println("  Eut/Est      = ", abs(Est) > 1e-24 ? Eut / Est : NaN, "   [~5–9]")
		println("  c_EFG_t/S50  = ", cG / S50, "   [~1.5–2.5]")
	else
		println("  (S50 ≈ 0 — ratios skipped)")
	end
	println("========== end EF pool audit ==========")
	println()
	return nothing
end
