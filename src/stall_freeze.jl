# ============================================================================
# stall_freeze.jl  —  freeze hard-stalled genes instead of integrating them
# ----------------------------------------------------------------------------
# When a gene's cognate ternary-complex pool is driven toward zero, elongation
# flux goes to zero (t3_efg_saturation) and per-codon dwell diverges. The gene
# is physically stationary: its ribosomes are jammed. The RHS itself is no
# longer singular (A); freezing is then only a way to skip a dead block.
#
# We therefore LATCH it as "frozen": its state is held constant (du = 0) and
# its block is skipped by the RHS. Because the state is held rather than
# deleted, its ribosomes, mRNA and post-termination complex stay sequestered
# and are NOT returned to the shared free pools — the resource-competition
# physics is preserved. This is the key distinction from simply dropping a
# gene, which would release the very demand that caused the stall.
#
# Freezing is only applied to genes that are already quasi-stationary, so the
# discontinuity introduced at the freeze instant is near-zero by construction.
#
# Entry points:
#   gene_stall_metrics(u_full, p, i)     — per-gene dwell / pool diagnostics
#   make_stall_callback(p; ...)          — PeriodicCallback that latches freezes
#   freeze_report(p, u_full)             — summary of what froze, when, and
#                                          how many ribosomes it locked up
# ============================================================================

"""
	gene_stall_metrics(u_full, p, i; with_du=false, t=0.0)

Per-gene stall diagnostics computed from a raw state vector `u_full`.

Uses the same intrinsic per-codon velocity as `spindler_codon_translation_time_at_u`
(`v = k_TLEmaxj * eff * t3_efg_saturation(cc, EGT, …) / 60`) but WITHOUT the
sliding-window smoothing: a stall confined to two or three codons must not be
diluted by its fast neighbours when we are deciding whether to freeze.

Returns a NamedTuple:
- `tau_max`   largest raw per-codon dwell time (s); `Inf` if a pool is empty
- `pos_max`   codon position carrying `tau_max`
- `codon_max` codon triplet at `pos_max`
- `cc_min`    smallest cognate ternary-complex pool across the ORF (µM)
- `KM_T3j`    ternary-complex Michaelis constant for this gene (µM)
- `rib_total` ribosomes held on this transcript (Σ occupancy + S70 + PTC, µM)
- `du_shared` ‖this gene's contribution to the shared metabolite rates‖∞ when
              `with_du`, else `NaN`. This is the size of the discontinuity that
              freezing the gene would introduce, and it is what gates freezing.
- `du_rib`    ‖d(ribosome occupancy)/dt‖∞ when `with_du`, else `NaN`. Diagnostic
              only — see the note in the body on why it is not the gate.
"""
function gene_stall_metrics(u_full, p, i::Int; with_du::Bool = false, t::Float64 = 0.0)
	tags  = p.local_tags_list[i]
	u_seq = extract_u_sequence(u_full, p, i)

	@unpack k_TLEmaxj, KM_T3j, KM_EFGGTP = tags.p_seq.EL
	EGT = max(u_seq[tags.idx.EGT], 0.0)

	num_species  = tags.num_species
	u_metab_safe = max.(@view(u_seq[1:num_species]), 0.0)
	C_vect       = reshape(u_metab_safe, 1, num_species)
	c_Codon      = tags.t3_Codon_struc * C_vect'

	efac = tags.flags.dynamic_efficiency ?
	       dynEffficiencyCalculation(tags, EGT, u_seq) :
	       tags.species_and_sequences.efac

	L       = tags.len_sequence - 2
	tau_max = 0.0
	cc_min  = Inf
	pos_max = 0

	@inbounds for j in 1:L
		cod_next = tags.species_and_sequences.Seq[j+1]
		cod_curr = tags.species_and_sequences.Seq[j]
		jn = tags.species_and_sequences.cod_idx[cod_next]
		kc = tags.species_and_sequences.cod_idx[cod_curr]

		cc = c_Codon[jn]
		cc < cc_min && (cc_min = cc)

		eff = efac[kc] * min(tags.efac_pos_aa[j], tags.efac_int_sd[j])

		tau = if eff > 0.0
			sat = t3_efg_saturation(cc, EGT, KM_T3j, KM_EFGGTP)
			v   = (k_TLEmaxj * eff * sat) / 60.0            # codons/s
			v > 0.0 ? 1.0 / v : Inf
		else
			Inf
		end

		if tau > tau_max
			tau_max = tau
			pos_max = j
		end
	end

	codon_max = pos_max > 0 ? String(tags.species_and_sequences.Seq[pos_max]) : ""

	# Ribosomes held on this transcript: elongating occupancies + the 70S
	# initiation complex + the post-termination complex. These stay out of the
	# shared free pool while the gene is frozen.
	S70       = max(u_seq[tags.idx.S70], 0.0)
	PTC       = max(u_seq[tags.idx.PTC], 0.0)
	c_RIB     = @view u_seq[tags.idx_first_codon:end]
	rib_total = sum(x -> max(x, 0.0), c_RIB) + S70 + PTC

	# Quasi-stationarity guards.
	#
	# du_shared is the one that gates freezing. It is this gene's own
	# contribution to the shared metabolite balance, isolated by subtracting a
	# second RHS call with the per-gene states zeroed (the cascade terms are
	# identical for every gene, so they cancel). That contribution is EXACTLY
	# what disappears from du[1:n_metab] when the gene is frozen, so its size is
	# the size of the discontinuity we are about to introduce.
	#
	# du_rib — the ribosome-occupancy norm — is reported for diagnostics but is
	# NOT the gate: rhs_sequence.jl zeroes derivatives at empty codons as a
	# positivity guard, so du_rib can read as ~0 for a gene that is still
	# drawing appreciably on the shared pools.
	n_metab   = p.vectors.n_metab
	du_rib    = NaN
	du_shared = NaN
	if with_du
		du_full = zeros(length(u_seq))
		ODE!(du_full, u_seq, tags, t)
		du_rib = maximum(abs, @view du_full[tags.idx_first_codon:end])

		u_bare = copy(u_seq)
		@views u_bare[n_metab+1:end] .= 0.0
		du_bare = zeros(length(u_seq))
		ODE!(du_bare, u_bare, tags, t)

		du_shared = 0.0
		@inbounds for j in 1:n_metab
			d = abs(du_full[j] - du_bare[j])
			d > du_shared && (du_shared = d)
		end
	end

	return (tau_max = tau_max, pos_max = pos_max, codon_max = codon_max,
	        cc_min = cc_min, KM_T3j = KM_T3j, rib_total = rib_total,
	        du_rib = du_rib, du_shared = du_shared)
end

"""
	make_stall_callback(p; tau_s, pool_frac, du_tol, check_every_min, protect, verbose)

Build a `PeriodicCallback` that latches hard-stalled genes as frozen.

The check runs every `check_every_min` simulation minutes (not every accepted
step) so its cost is bounded, and it runs BETWEEN steps so the right-hand side
stays continuous within a step.

A gene is frozen when all three hold:
1. `tau_max   > tau_s`               — a hard stall, not merely a slow codon
2. `cc_min    < pool_frac * KM_T3j`  — caused by cognate tRNA exhaustion
3. `du_shared < du_tol`              — its draw on the shared pools has ceased,
                                       so removing it is a near-no-op rather
                                       than an artificial jump

Freezing is monotone: a frozen gene is never thawed, which prevents on/off
chattering that would defeat the implicit solver.

`protect` is a set of gene indices that must never be frozen — always include
the gene under study, otherwise the measurement is destroyed. At least one gene
is always left active.
"""
function make_stall_callback(p;
		tau_s::Float64           = 100.0,
		pool_frac::Float64       = 0.01,
		du_tol::Float64          = 1e-6,
		check_every_min::Float64 = 1.0,
		protect                  = Set{Int}(),
		verbose::Bool            = true)

	v        = p.vectors
	protect_ = Set{Int}(protect)
	dt_check = max(check_every_min, 0.05)

	return PeriodicCallback(dt_check; initial_affect = false) do integrator
		u = integrator.u
		t = integrator.t

		newly = 0
		for i in 1:v.n_seq
			v.frozen[i] && continue
			i in protect_ && continue

			# Never freeze the last active gene — the shared-cascade correction
			# and every downstream metric assume at least one live sequence.
			(v.n_active[] - newly) <= 1 && break

			# Two-stage: the cheap dwell/pool screen first, and only for a
			# candidate do we pay for an extra RHS evaluation to confirm it
			# has actually stopped moving.
			m = gene_stall_metrics(u, p, i)
			(m.tau_max > tau_s && m.cc_min < pool_frac * m.KM_T3j) || continue

			m = gene_stall_metrics(u, p, i; with_du = true, t = t)

			if m.du_shared < du_tol
				v.frozen[i]      = true
				v.freeze_time[i] = t
				newly += 1

				if verbose
					gname = String(p.local_tags_list[1].species_and_sequences.seq_names[i])
					@info "stall freeze" gene = gname t_min = round(t, digits = 3) codon = m.codon_max position = m.pos_max tau_s = m.tau_max pool_uM = m.cc_min ribosomes = m.rib_total du_shared = m.du_shared du_rib = m.du_rib
				end
			end
		end

		if newly > 0
			v.n_active[] = count(!, v.frozen)
		end

		# A freeze changes the RHS discontinuously: the frozen block's
		# derivatives vanish, and the shared-cascade over-count factor drops by
		# one, which shifts the shared metabolite rates by exactly that gene's
		# own contribution. The integrator's cached history and Jacobian are
		# therefore stale and it MUST reinitialise — telling it otherwise makes
		# the Newton solve stop converging and the step size collapse. When
		# nothing froze we only read u, so the reinit can be skipped.
		u_modified!(integrator, newly > 0)
		return nothing
	end
end

"""
	freeze_report(p, u_full)

Summarise the outcome of stall freezing for a completed (or aborted) run.

Returns a NamedTuple:
- `n_frozen`              number of genes frozen
- `n_active`              number still integrating
- `frozen_genes`          names of frozen genes
- `freeze_times`          sim-time (min) at which each froze
- `freeze_codons`         codon triplet carrying the stall at freeze detection
- `freeze_positions`      codon position of that stall
- `ribosomes_sequestered` total ribosomes held on frozen transcripts (µM)
- `rows`                  per-gene NamedTuples, convenient for `DataFrame(...)`
"""
function freeze_report(p, u_full)
	v     = p.vectors
	names = p.local_tags_list[1].species_and_sequences.seq_names

	rows = NamedTuple[]
	rib_seq = 0.0

	for i in 1:v.n_seq
		v.frozen[i] || continue
		m = gene_stall_metrics(u_full, p, i)
		rib_seq += m.rib_total
		push!(rows, (gene        = String(names[i]),
		             seq_index   = i,
		             freeze_time = v.freeze_time[i],
		             position    = m.pos_max,
		             codon       = m.codon_max,
		             tau_s       = m.tau_max,
		             pool_uM     = m.cc_min,
		             ribosomes   = m.rib_total))
	end

	sort!(rows, by = r -> r.freeze_time)

	return (n_frozen              = length(rows),
	        n_active              = v.n_active[],
	        frozen_genes          = String[r.gene for r in rows],
	        freeze_times          = Float64[r.freeze_time for r in rows],
	        freeze_codons         = String[r.codon for r in rows],
	        freeze_positions      = Int[r.position for r in rows],
	        ribosomes_sequestered = rib_seq,
	        rows                  = rows)
end

"""
	reset_freeze_state!(p)

Clear all latched freeze flags. Called before each solve so a reused parameter
object cannot carry freezes over from a previous run.
"""
function reset_freeze_state!(p)
	v = p.vectors
	fill!(v.frozen, false)
	fill!(v.freeze_time, NaN)
	v.n_active[] = v.n_seq
	return nothing
end

"""
	print_freeze_summary(rep)

Compact console table of frozen genes, printed after a run.
"""
function print_freeze_summary(rep)
	if rep.n_frozen == 0
		println("Stall freezing: no genes frozen (", rep.n_active, " active).")
		return nothing
	end

	println("\nStall freezing: ", rep.n_frozen, " gene(s) frozen, ",
	        rep.n_active, " still active")
	@printf("  %-18s %10s %8s %8s %12s %12s\n",
	        "gene", "t_freeze", "pos", "codon", "tau (s)", "ribosomes")
	for r in rep.rows
		@printf("  %-18s %10.2f %8d %8s %12.3g %12.4g\n",
		        r.gene, r.freeze_time, r.position, r.codon, r.tau_s, r.ribosomes)
	end
	@printf("  ribosomes sequestered in frozen transcripts: %.4g µM\n",
	        rep.ribosomes_sequestered)
	return nothing
end
