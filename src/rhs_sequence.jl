# like: MOD_seq_rhs.m 
# Start: 14.04.2025
# Last edit: 
# By: Cedric Meyer 
# Purpose: ODE system of the sequence module; this ODE deviates from the matlab model by indexing of the efficiency factors and q_R_elo

function rhs_sequence!(du, u, p)


    # du = zeros(Float64, length(u))
    tags = p
    # Unpack parameters 
    k_cat_Pol, k_D, C_RNP, V_TCmax, K_M_D, C_RNE, k_cat_E, K_M_deg, k_TLEmaxj, KM_T3j, KM_EFGGTP, KM_RK, KM_RF3GTP, k_RF1_diss_GTP, k_RF3_diss_GDP, k_GTP_hydr, k_RF1_ass, k_RF3_GDP_ass, k_RF3_GTP_ass, k_RF3_GDP_diss, k_RF3_GTP_diss, A1, A2, A3, A4, A5, k_cat_K_M_G1, k_cat_K_M_G2, k_RRF, q_RRF, one_q_G2_k_split, k_TE = tags.p_seq

    # Unpack constant variables 
    ATP, GTP, UTP, CTP, ADP, GDP, UDP, CDP, AMP, GMP, UMP, CMP, Eut, Est, RFt, RRF, Ala, Cys, Asp, Glu, Phe, Gly, His, Ile, Lys, Leu, Met, Asn, Pro, Gln, Arg, Ser, Thr, Val, Trp, Tyr, C_D = tags.constant_vars 

    # Unpack dynamic variables 
    S30, S50, IF1, IF2, IF3, SC1, SC2, SC3, SC4, PIC, EFG, EFD, EGT, EGD, EFT, RF3, RGT, RGD, b1t, ta1, ta2, tc1, td1, te2, tf1, tg2, tg3, th1, ti1, tk1, t1l, t2l, t3l, t4l, t5l, tm1, tn1, tp1, tp2, tp3, tq1, tq2, tr2, tr3, tr4, tr5, ts1, ts2, ts3, ts5, tt1, tt2, tt4, tv1, tv2, tw1, ty1, tA1, tA2, tC1, tD1, tE2, tF1, tG2, tG3, tH1, tI1, tK1, tL1, tL2, tL3, tL4, tL5, tM1, tN1, tP1, tP2, tP3, tQ1, tQ2, tR2, tR3, tR4, tR5, tS1, tS2, tS3, tS5, tT1, tT2, tT4, tV1, tV2, tW1, tY1, a1t, a2t, c1t, d1t, e2t, f1t, g2t, g3t, h1t, i1t, k1t, l1t, l2t, l3t, l4t, l5t, m1t, n1t, p1t, p2t, p3t, q1t, q2t, r2t, r3t, r4t, r5t, s1t, s2t, s3t, s5t, t1t, t2t, t4t, v1t, v2t, w1t, y1t, RNA, S70, rna, PRT, PTC = u 

    # Containers 
    num_species = length(tags.Name)
    len_sequence = length(tags.Seq) # Codon length with stop codon

    # Define RIB1
    idx_first_codon = length(tags.Tag)+1
    rate_RIB        = zeros(Float64, len_sequence-1);
    c_RIB	 		= u[idx_first_codon:end]

    # Transcription -----------------

    # Free mRNA (RNA) levels are given by a starting value to implement mRNA concentrations from transcriptomics measurements

    n_DNA_t = 1;
    n_DNA_s = 1;

    # # # tQSSR model used in Spindler

    V_TCmax 	= k_cat_Pol * C_RNP
    V_degmax    = k_cat_E * C_RNE

    du[tags.idx["RNA"]] += (V_TCmax/(len_sequence*3) * C_D * tags.scale.V_TC_scale / ( K_M_D * (1 + (n_DNA_t-n_DNA_s)* 0.001/K_M_D)+ C_RNP + C_D * tags.scale.V_TC_scale )) - 
    (V_degmax * RNA / (K_M_deg  + C_RNE + RNA ))	# 170

    ## translation elongation ----------------------------------------

    #= queuing factor q_R_elo
        Nieß 2017:
        Any influence of ribosomal hindrance was incorporated by applying a queuing factor q_R_elo,
        which reduces the elongation rate with increasing ribosome density. For an elongation step from codon 
        j to j+1, the queuing factor corresponds to the density of ribosomes on  the 12 codons further downstream. =#         

    total_RNA   = RNA + sum(c_RIB[1:12])

    q_R_elo     = ones(len_sequence-1);
        
    for i in 2:(len_sequence - 13)

        if total_RNA-sum(c_RIB[i:i+12]) == 0
            q_R_elo[i] 	= 0;
        else
            q_R_elo[i] 	= 1 - c_RIB[i+12] / (total_RNA-sum(c_RIB[i:i+11]));
        
        end
    end 

        
    for i in (len_sequence-13) : (len_sequence-1)
        if total_RNA - sum(c_RIB[i : end]) == 0
            q_R_elo[i] = 0
        else 
            q_R_elo[i] = 1 - (c_RIB[end] + PTC) / (total_RNA - sum(c_RIB[i:end-1]))
        end
    end


        # dynamic efficiency factor 
    if tags.dynamic_efficiency == true

        c_T3 = u[tags.idx["a1t"] : tags.idx["y1t"]];
        efac = zeros(Float64, 62,1); 

        WC = tags.T_WC * c_T3;
        WB = tags.T_WB * c_T3;
        NC = tags.T_NC * c_T3;
        NN = tags.T_NN * c_T3;

        efac = dynEffficiencyCalculation(WC,WB,NC,NN,EGT);

        tags.efac = efac
    end
        

    # V_EL
    C_vect = reshape(u[1:num_species], 1, num_species)  # Directly assign
    c_Codon = tags.t3_Codon_struc * C_vect'

    V_TL_EL	= zeros(len_sequence-1);

    EGT_inv = inv(EGT)

    # rna_free = u[170]
    # s70_ini = u[171]
    # # this calculation is done according to the analyze_parllel script and is fitting to the description in the Nieß paper
    # q_R_elo[1] = 1 - (sum(c_RIB[1:12]) +s70_ini)/(rna_free+sum(c_RIB[1:12])) 


    for i in 1: len_sequence-2

        j = findfirst(x -> x == String(tags.Seq[i+1]), tags.cod)
        k = findfirst(x -> x == String(tags.Seq[i]), tags.cod)

        if c_Codon[j] == 0
            V_TL_EL[i+1] = 0;  
        else 
            efficiency = tags.efac[k] * min(tags.efac_pos_aa[i], tags.efac_int_sd[i]) ############### The indexing is was altered in comparison to the matlab model from i+1 to i. This is in agreement with the spindler paper and the descriptions in Arnold and Nieß. This is also how it is calculated in analyze_parallel_data.m.
            V_TL_EL[i+1] =  q_R_elo[i] * k_TLEmaxj * c_RIB[i] *  efficiency /(1+ KM_T3j/c_Codon[j]+ KM_EFGGTP*EGT_inv); ########### Indexing of q_R_elo was altered from i+1 to i. Now its same as in analyze_paralle_data data and the description in Nieß paper.  
        end

    end
       
    

    # V_TL_EL     .= V_TL_EL .* q_R_elo; # original version that is like in the mod_seq matlab script

    V_TL_EL_sum = sum(V_TL_EL);

    c_codon_seq = tags.matrix_N * c_Codon
    N = C_vect[1,:] * (1 ./c_codon_seq)';

    # Replace Inf values with 0
    N[isinf.(N)] .= 0
    # Replace NaN values with 0
    N[isnan.(N)] .= 0

    M = zeros(size(N))

    for i in 1:size(N,1)
        for k in 1:size(N,2)
            M[i,k] = tags.alpha_struc[i,k] * N[i,k]
        end
    end


    K_T3 = M*V_TL_EL;

    mtRNA 		= zeros(size(M,1),2);
    M_tRNA 		= hcat(mtRNA,M)
    M_tRNA		= M_tRNA[:, Not(3)]
    M_tRNA		= M_tRNA[:, 1:end-1]

    K_tRNA 		= M_tRNA*V_TL_EL   

    ## translation termination
    te = len_sequence-1 # Number of last coding codon
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
    V_TL_TE = V_TE_max * efficiency_term / (1 + KM_RK/ c_RIB[te] + KM_RF3GTP/RGT + KM_RK * KM_RF3GTP /(c_RIB[te] * RGT))
        
    if c_RIB[te] == 0 || GTP == 0
        V_TL_TE=0;
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
        
    ## translation termination - calculation of du 
    tRNA_K = tags.alpha_struc[:,end] * V_TL_TE

    du[1:num_species] .+= tRNA_K;

    rate_RIB[len_sequence-1] += - V_TL_TE

    ######### !!!!!!!!!!! Indexing of output variables needs to be adjusted to single sequence ??????????????? Or is it????????????
    du[tags.idx["PTC"]] += V_TL_TE - k_rec_min * PTC;	            # PTC
    du[tags.idx["S30"]] += k_rec_min * PTC;                         # rate of 30S release due to recycling
    du[tags.idx["S50"]] +=  k_rec_min * PTC;                        # rate of 50S release due to recycling
    du[tags.idx["PRT"]] = V_TL_TE;                                  # rate of PRT (Protein)
    du[tags.idx["RGT"]] = -V_TL_TE -V_RF3_GTP_diss +V_RF3_GTP_ass;  # rate of RGT
    du[tags.idx["RGD"]] =  V_TL_TE -V_RF3_GDP_diss +V_RF3_GDP_ass;  # rate of RGD
    du[tags.idx["RF3"]] =          +V_RF3_GDP_diss +V_RF3_GTP_diss -V_RF3_GTP_ass -V_RF3_GDP_ass; # rate of RF3

    # translation elongation - calculation of du 

    if RNA == 0 # concentration of free mRNA
        du = du
    else

        # concentration of the ribosome at the specific codon
        rate_RIB[1] = - V_TL_EL[2];  

        for i in 2: len_sequence-2
            rate_RIB[i] = V_TL_EL[i] - V_TL_EL[i+1];	
        end

        rate_RIB[len_sequence-1] = rate_RIB[len_sequence-1] + V_TL_EL[end];

        # concentration of the unloaded tRNA
        du[tags.idx["ta1"]:tags.idx["ty1"]] .+= K_tRNA[end-tags.no_of_tRNAs+1:end]

        # concentration of the ternary-complex
        du[1: num_species] .+= - K_T3;
        du[tags.idx["EGT"]]  += - V_TL_EL_sum - f_s * k_rec_min * PTC;	# Rate of EFG_GTP 
        du[tags.idx["EGD"]]  += + V_TL_EL_sum + f_s * k_rec_min * PTC;	# Rate of EFG_GDP
        du[tags.idx["RNA"]] += + V_TL_EL[12];
        du[tags.idx["rna"]] += - V_TL_EL[12];

    end

    # save RIB in du
    du[idx_first_codon:end] .+= rate_RIB;

    return du
end