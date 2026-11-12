# like: MOD_seq_rhs.m 
# Start: 14.04.2025
# Last edit: 
# By: Cedric and Arnav 
# Purpose: ODE system of the sequence module; this ODE deviates from the matlab model by indexing of the efficiency factors and tags.tags.q_R_elo

function rhs_sequence!(du, u, p)

    # du = zeros(Float64, length(u))
    tags = p
    # Unpack parameters 
    @unpack µ, k_proc, k_cat_Pol, C_RNP, V_TCmax, K_M_D, k_cat_E, C_RNE, V_degmax, K_M_deg, k_D, n_DNA_t, n_DNA_s = tags.p_seq.TX
    @unpack k_TLEmaxj, KM_T3j, KM_EFGGTP = tags.p_seq.EL
    @unpack k_TE, KM_RK, KM_GTP, KM_RF3GTP, k_RF3_diss_GDP, k_RF1_diss_GTP, k_GTP_hydr, k_RF1_ass = tags.p_seq.TER
    @unpack k_RF3_GDP_ass, k_RF3_GTP_ass, k_RF3_GDP_diss, k_RF3_GTP_diss = tags.p_seq.RF3
    @unpack A1, A2, A3, A4, A5, k_cat_K_M_G1, k_cat_K_M_G2, k_RRF, q_RRF, one_q_G2_k_split = tags.p_seq.PTC 
    
    # Unpack constant variables 
    ATP, GTP, UTP, CTP, ADP, GDP, UDP, CDP, AMP, GMP, UMP, CMP, Eut, Est, RFt, RRF, Ala, Cys, 
    Asp, Glu, Phe, Gly, His, Ile, Lys, Leu, Met, Asn, Pro, Gln, Arg, Ser, Thr, val, Trp, Tyr, C_D = tags.init_c.constant_vars 

    # Unpack dynamic variables 
    S30, S50, IF1, IF2, IF3, SC1, SC2, SC3, SC4, PIC, EFG, EFD, EGT, EGD, EFT, RF3, RGT, RGD, b1t, ta1, ta2, tc1, td1, te2, tf1, tg2, tg3, th1, ti1, tk1, t1l, t2l, t3l, t4l, t5l, tm1, tn1, tp1, tp2, tp3, 
    tq1, tq2, tr2, tr3, tr4, tr5, ts1, ts2, ts3, ts5, tt1, tt2, tt4, tv1, tv2, tw1, ty1, tA1, tA2, tC1, tD1, tE2, tF1, tG2, tG3, tH1, tI1, tK1, tL1, tL2, tL3, tL4, tL5, tM1, tN1, tP1, tP2, tP3, tQ1, tQ2, 
    tR2, tR3, tR4, tR5, tS1, tS2, tS3, tS5, tT1, tT2, tT4, tV1, tV2, tW1, tY1, a1t, a2t, c1t, d1t, e2t, f1t, g2t, g3t, h1t, i1t, k1t, l1t, l2t, l3t, l4t, l5t, m1t, n1t, p1t, p2t, p3t, q1t, q2t, r2t, r3t, 
    r4t, r5t, s1t, s2t, s3t, s5t, t1t, t2t, t4t, v1t, v2t, w1t, y1t, S70, rna, PRT, PTC = u 


    # Define RIB1

    fill!(tags.rate_RIB, 0.0)

    @views c_RIB	= u[tags.idx_first_codon:end]     

    # Positivity guard: the implicit solver can produce tiny negative ribosome occupancies.
    # Using negative occupancies directly in flux formulas can amplify the error and create runaway negatives.
    # We compute fluxes using non-negative ribosome occupancies.
    total_RNA = max(rna, 0.0)  # Dynamic mRNA from state vector (was tags.RNA)
    @inbounds for k in 1:12
        total_RNA += max(c_RIB[k], 0.0)
    end

    fill!(tags.q_R_elo, 1.0)
        
    @views for i in 2:(tags.len_sequence - 13)

        # sum of the window, treating negative occupancies as 0
        window_sum = 0.0
        @inbounds for k in i:(i+12)
            window_sum += max(c_RIB[k], 0.0)
        end

        if total_RNA - window_sum == 0

            tags.q_R_elo[i] 	= 0;
        else
            denom = total_RNA
            @inbounds for k in i:(i+11)
                denom -= max(c_RIB[k], 0.0)
            end
            # Clamp to [0,1]: the occupancy factor is a probability.
            # Without clamping, high local ribosome density makes denom tiny,
            # c_RIB[i+12]/denom >> 1, q_R_elo becomes hugely negative,
            # producing unphysical reverse elongation flux (TASEP mean-field breakdown).
            tags.q_R_elo[i] = clamp(1 - max(c_RIB[i+12], 0.0) / denom, 0.0, 1.0)
        
        end
    end 


    @views for i in (tags.len_sequence-13) : (tags.len_sequence-1)
        tail_sum = 0.0
        @inbounds for k in i:length(c_RIB)
            tail_sum += max(c_RIB[k], 0.0)
        end
        if total_RNA - tail_sum == 0
            tags.q_R_elo[i] = 0
        else 
            denom = total_RNA
            @inbounds for k in i:(length(c_RIB)-1)
                denom -= max(c_RIB[k], 0.0)
            end
            tags.q_R_elo[i] = clamp(1 - (max(c_RIB[end], 0.0) + max(PTC, 0.0)) / denom, 0.0, 1.0)
        end
    end

    # dynamic efficiency factor 
    if tags.flags.dynamic_efficiency == true

        efac = dynEffficiencyCalculation(tags, EGT, u);

    end
        

    # V_EL
    # Clamp metabolite concentrations to ≥ 0 before computing effective T3 pools.
    # Negative T3 complexes (numerical drift) would make c_Codon ≈ 0 → division
    # blow-up in matrix_free_KT3!. Elongation itself uses t3_efg_saturation
    # (finite at empty pools); this clamp is still required for the T3 mapping.
    u_metab_safe = max.(@view(u[1:tags.num_species]), 0.0)
    C_vect = reshape(u_metab_safe, 1, tags.num_species)
    c_Codon = tags.t3_Codon_struc * C_vect'

    V_TL_EL	= zeros(tags.len_sequence-1);

    EGT_pos = max(EGT, 0.0)

    @inbounds for i in 1:tags.len_sequence-2

        cod_next = tags.species_and_sequences.Seq[i+1]
        cod_curr = tags.species_and_sequences.Seq[i]

        j = tags.species_and_sequences.cod_idx[cod_next]    # O(1) Lookup
        k = tags.species_and_sequences.cod_idx[cod_curr]    # O(1) Lookup

        cc = c_Codon[j]
        sat = t3_efg_saturation(cc, EGT_pos, KM_T3j, KM_EFGGTP)
        if sat == 0.0
            V_TL_EL[i+1] = 0
        else
            eff = efac[k] * min(tags.efac_pos_aa[i], tags.efac_int_sd[i])
            V_TL_EL[i+1] =
                tags.q_R_elo[i] * k_TLEmaxj * max(c_RIB[i], 0.0) * eff * sat
        end
    end

    V_TL_EL_sum = sum(V_TL_EL);
    c_codon_seq = c_Codon[tags.matrix_N]

    matrix_free_KT3!(tags, V_TL_EL, c_codon_seq, C_vect)

     
    ## translation termination
    te = tags.len_sequence-1 # Number of last coding codon
    # ke = findfirst(x -> x == tags["Seq1"]["Codon"][te], tags["cod"])

    k_RF1_diss_GTP = k_RF1_diss_GTP#*tags.scale.termination

    efficiency_term = 1;
    k_ter_RF1 = 1 / ( 1/(k_RF1_ass*RFt) + 1/(k_RF1_diss_GTP) );
    k_ter_RF3 = 1 / ( 1/(k_RF1_ass*RGT) + 1/(k_RF3_diss_GDP) + 1/(k_GTP_hydr) );

    if k_ter_RF1 < k_ter_RF3
        k_ter = k_ter_RF1;
    else
        k_ter = k_ter_RF3;
    end


    V_TE_max = k_ter * RFt   # RFt is RF1/2 total
    c_te = max(c_RIB[te], 0.0)
    if c_te == 0.0 || GTP == 0
        V_TL_TE = 0.0
    else
        V_TL_TE = V_TE_max * efficiency_term / (1 + KM_RK / c_te + KM_RF3GTP / RGT + KM_RK * KM_RF3GTP / (c_te * RGT))
    end
        
    V_RF3_GDP_ass  = k_RF3_GDP_ass  * RF3 * GDP;
    V_RF3_GTP_ass  = k_RF3_GTP_ass  * RF3 * GTP;
    V_RF3_GDP_diss = k_RF3_GDP_diss * RGD;
    V_RF3_GTP_diss = k_RF3_GTP_diss * RGT 

    ## post termination complex recycling
    t_rec = A1 + A2 * 1/RRF + A3 * EGT/RRF + A4 * 1/EGT + A5 * 1/(RRF*EGT); # time/s for recycling according to Borg et al. 2016
    k_rec = 1/t_rec; # rate/s^-1 for recycling according to Borg et al. 2016
    k_rec_min = k_rec*60; # rate/min^-1 for recycling according to Borg et al. 2016
    f_s = k_cat_K_M_G1/(k_RRF * RRF)*(q_RRF/k_cat_K_M_G2*one_q_G2_k_split + EGT)+ one_q_G2_k_split; #EGT hydrolised per splitting event according to Borg et al 2016

    # translation termination - calculation of du 

    du[tags.alpha_one[end]] += V_TL_TE

    tags.rate_RIB[tags.len_sequence-1] -= V_TL_TE

    du[tags.idx.PTC] += V_TL_TE - k_rec_min * PTC;	            # PTC
    du[tags.idx.S30] += k_rec_min * PTC;                         # rate of 30S release due to recycling
    du[tags.idx.S50] +=  k_rec_min * PTC;                        # rate of 50S release due to recycling
    du[tags.idx.PRT] = V_TL_TE;                                  # rate of PRT (Protein)
    du[tags.idx.RGT] = -V_TL_TE -V_RF3_GTP_diss +V_RF3_GTP_ass;  # rate of RGT
    du[tags.idx.RGD] =  V_TL_TE -V_RF3_GDP_diss +V_RF3_GDP_ass;  # rate of RGD
    du[tags.idx.RF3] =          +V_RF3_GDP_diss +V_RF3_GTP_diss -V_RF3_GTP_ass -V_RF3_GDP_ass; # rate of RF3


    # translation elongation - calculation of du 

    # concentration of the ribosome at the specific codon
    tags.rate_RIB[1] = - V_TL_EL[2];  

    for i in 2: tags.len_sequence-2
        tags.rate_RIB[i] = V_TL_EL[i] - V_TL_EL[i+1];	
    end

    tags.rate_RIB[tags.len_sequence-1] = tags.rate_RIB[tags.len_sequence-1] + V_TL_EL[end];

    # concentration of the unloaded tRNA
    du[tags.idx.ta1:tags.idx.ty1] .+= tags.K_tRNA[end-tags.species_and_sequences.no_of_tRNAs+1:end]

    # concentration of the ternary-complex
    du[1: tags.num_species] .+= - tags.K_T3;
    if !tags.flags.static_EFs
        du[tags.idx.EGT]  += - V_TL_EL_sum - f_s * k_rec_min * PTC;	# Rate of EFG_GTP 
        du[tags.idx.EGD]  += + V_TL_EL_sum + f_s * k_rec_min * PTC;	# Rate of EFG_GDP
    end

    # mRNA dynamics: ribosome clearing the RBS frees mRNA for new initiation
    # (matches MATLAB: rate.RNA += V_TL_EL(12), rate.rna -= V_TL_EL(12))
    if length(V_TL_EL) >= 12
        du[tags.idx.rna] += V_TL_EL[12]
    end

    rna_conc = max(u[tags.idx.rna], 0.0)

    if tags.flags.spindler_tx
        # Spindler et al. 2025 transcription / degradation with EXPLICIT competition for
        # the shared RNA-polymerase (C_RNP) and RNase E (C_RNE) pools — faithful port of
        # MOD_seq_rhs.m L22-43.  The polymerase pool in the denominator is what caps total
        # transcription (and lets the active-ribosome fraction settle below 100 %); the old
        # formula dropped both C_D and C_RNP, so mRNA over-produced and saturated the pool.
        #
        # Per-gene DNA dosage (Spindler's V_TC_scale(ID)) and the total over all genes are
        # carried on the local tags, so the RNA-polymerase competition term is now EXACT for
        # heterogeneous dosing (e.g. rplL=340 vs eCU=1133) rather than the uniform-count
        # approximation.  V_TCmax is nuc·µM/min; divide by L_nuc to get µM_mRNA/min.
        L_nuc    = Float64(3 * tags.len_sequence)
        n_seq    = tags.species_and_sequences.no_of_seqs
        s_dose   = tags.V_TC_dose                                 # this gene's DNA dosage
        tot_dose = tags.V_TC_dose_total                           # Σ dosage over all genes
        # n_DNA_t − n_DNA_s : competing DNA (Σ_other − self); 0 for a single gene (MATLAB L22-24)
        comp_DNA = n_seq <= 1 ? 0.0 : (tot_dose - s_dose) - s_dose
        comp_RNA = n_seq <= 1 ? 0.0 : (n_seq - 1) * rna_conc      # Σ_other RNA ≈ (n_seq−1)·rna
        V_TC_s  = (V_TCmax / L_nuc) * C_D * s_dose /
                  (K_M_D * (1.0 + comp_DNA * 0.001 / K_M_D) + C_RNP + C_D * s_dose)
        V_deg_s = (rna_conc > 0.0) ?
                  V_degmax * rna_conc /
                  (K_M_deg * (1.0 + comp_RNA / K_M_deg) + C_RNE + rna_conc) : 0.0
        du[tags.idx.rna] += V_TC_s - V_deg_s
    else
        # Michaelis-Menten production/degradation (matches MATLAB V_TC / V_deg time scale).
        # Steady state: rna → tags.RNA.  Effective δ ≈ V_degmax/(K_M_deg + C_RNE + RNA)
        # ≈ 5 /min, fast enough to balance V_IN_70SIC consumption from rhs_metabolites.
        K_eff = K_M_deg + C_RNE
        production_rate = V_degmax * tags.RNA / (K_eff + tags.RNA)
        deg_flux = (rna_conc > 0.0) ? V_degmax * rna_conc / (K_eff + rna_conc) : 0.0
        du[tags.idx.rna] += production_rate - deg_flux
    end

    # save RIB in du
    du[tags.idx_first_codon:end] .+= tags.rate_RIB;

    # Final positivity guard for ribosome-position states:
    # if a ribosome position is at/below zero, do not allow its derivative to drive it further negative.
    @inbounds for j in tags.idx_first_codon:length(u)
        if u[j] <= 0.0 && du[j] < 0.0
            du[j] = 0.0
        end
    end

    return du
end

