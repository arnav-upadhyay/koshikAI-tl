# long_met_rhs_01
# like: MOD_met_rhs.m 
# Start: 13.03.2025
# By: Cedric Meyer 
# Purpose: ODE system of the metabolite module, long form without the loops 

# function ODE_met_long(du, u, p::data_struct, t)

# Can mod_met_rhs and mod_seq_rhs be combined so that du can be directly saved into du?

function rhs_metabolites!(du, u, p)


	# du =zeros(Float64, length(u))
	tags = p

	# Unpack parameters
	# @unpack q_RBS, k_IN_70SIC, KM_mRNA, KM_50S, k_IN_IF2D, k_on_30S_IF3, k_off_30SIF3, k_on_30SIF3_IF1, k_off_30SIF3IF1, k_on_30SIF3_IF2, k_off_30SIF3IF2, k_on_30SIF3IF1_IF2, k_off_30SIF3IF1IF2, k_on_30SIF3IF2_IF1, k_off_30SIF3IF2IF1, k_on_30SIF123_fMet, k_off_30SIF123fMet, k_EFGGDP_Diss, k_EFG_GDP_Ass, k_EFG_GTP_Ass, k_EFGGTP_Diss, k_EFTsf, k_EFTsr, Keq_EFTu, KM_GTP, KM_GDP, KM_EFTuGDP, KM_EFTuGTP, Ki_EFTuGDP, Ki_EFTuGTP, c_EFG_t, k_cat, KM_ARS_tRNA, KM_ARS_aa, KM_ARS_ATP, k_T3j_f, k_T3j_d = tags.p_met

	@unpack k_IN_70SIC, KM_mRNA, KM_50S, k_IN_IF2D, k_on_30S_IF3, k_off_30SIF3, k_on_30SIF3_IF1, k_off_30SIF3IF1, k_on_30SIF3_IF2, k_off_30SIF3IF2, k_on_30SIF3IF1_IF2, k_off_30SIF3IF1IF2, k_on_30SIF3IF2_IF1, k_off_30SIF3IF2IF1, k_on_30SIF123_fMet, k_off_30SIF123fMet = tags.p_met.IN
	@unpack k_EFGGDP_Diss, k_EFG_GDP_Ass, k_EFG_GTP_Ass, k_EFGGTP_Diss, k_EFTsf, k_EFTsr, Keq_EFTu, KM_GTP, KM_GDP, KM_EFTuGDP, KM_EFTuGTP, Ki_EFTuGDP, Ki_EFTuGTP, c_EFG_t = tags.p_met.EF
	@unpack k_cat, KM_ARS_tRNA, KM_ARS_aa, KM_ARS_ATP = tags.p_met.tRNA
	@unpack k_T3j_f, k_T3j_d = tags.p_met.t3 
	q_RBS = tags.p_met.q_RBS

	# Unpack constant variables 
	ATP, GTP, UTP, CTP, ADP, GDP, UDP, CDP, AMP, GMP, UMP, CMP, Eut, Est, RFt, RRF, Ala, Cys, Asp, Glu, Phe, Gly, His, Ile, Lys, Leu, Met, Asn, Pro, Gln, Arg, Ser, Thr, val, Trp, Tyr, C_D = tags.init_c.constant_vars

	# Unpack dynamic variables 
	S30, S50, IF1, IF2, IF3, SC1, SC2, SC3, SC4, PIC, EFG, EFD, EGT, EGD, EFT, RF3, RGT, RGD, b1t, ta1, ta2, tc1, td1, te2, tf1, tg2, tg3, th1, ti1, tk1, t1l, t2l, t3l, t4l, t5l, tm1, tn1, tp1, tp2, tp3, tq1, tq2, tr2, tr3, tr4, tr5, ts1, ts2, ts3, ts5, tt1, tt2, tt4, tv1, tv2, tw1, ty1, tA1, tA2, tC1, tD1, tE2, tF1, tG2, tG3, tH1, tI1, tK1, tL1, tL2, tL3, tL4, tL5, tM1, tN1, tP1, tP2, tP3, tQ1, tQ2, tR2, tR3, tR4, tR5, tS1, tS2, tS3, tS5, tT1, tT2, tT4, tV1, tV2, tW1, tY1, a1t, a2t, c1t, d1t, e2t, f1t, g2t, g3t, h1t, i1t, k1t, l1t, l2t, l3t, l4t, l5t, m1t, n1t, p1t, p2t, p3t, q1t, q2t, r2t, r3t, r4t, r5t, s1t, s2t, s3t, s5t, t1t, t2t, t4t, v1t, v2t, w1t, y1t, S70, rna, PRT, PTC, RIB1 = u 

	RNA = max(rna, 0.0)  # Dynamic mRNA from state vector (was tags.RNA)

	# Guard physically nonnegative species against numerical drift (local copies only — solver state untouched).
	# Unloaded tRNAs appear as KM/tRNA in denominators: negative value → denominator shrinks → rate explodes.
	_eps = 1e-12
	ta1 = max(ta1, _eps); ta2 = max(ta2, _eps); tc1 = max(tc1, _eps); td1 = max(td1, _eps)
	te2 = max(te2, _eps); tf1 = max(tf1, _eps); tg2 = max(tg2, _eps); tg3 = max(tg3, _eps)
	th1 = max(th1, _eps); ti1 = max(ti1, _eps); tk1 = max(tk1, _eps)
	t1l = max(t1l, _eps); t2l = max(t2l, _eps); t3l = max(t3l, _eps); t4l = max(t4l, _eps); t5l = max(t5l, _eps)
	tm1 = max(tm1, _eps); tn1 = max(tn1, _eps)
	tp1 = max(tp1, _eps); tp2 = max(tp2, _eps); tp3 = max(tp3, _eps)
	tq1 = max(tq1, _eps); tq2 = max(tq2, _eps)
	tr2 = max(tr2, _eps); tr3 = max(tr3, _eps); tr4 = max(tr4, _eps); tr5 = max(tr5, _eps)
	ts1 = max(ts1, _eps); ts2 = max(ts2, _eps); ts3 = max(ts3, _eps); ts5 = max(ts5, _eps)
	tt1 = max(tt1, _eps); tt2 = max(tt2, _eps); tt4 = max(tt4, _eps)
	tv1 = max(tv1, _eps); tv2 = max(tv2, _eps); tw1 = max(tw1, _eps); ty1 = max(ty1, _eps)
	# EF-Tu·GTP (EFT): negative → Michaelis denom shrinks, T3 terms flip sign
	EFT = max(EFT, 0.0)
	# EFG·GDP (EGD): negative → dissociation term flips, drives huge |du|
	EGD = max(EGD, 0.0)
	# Ribosomal subunits and initiation intermediates
	S30 = max(S30, 0.0); S50 = max(S50, 0.0); S70 = max(S70, 0.0); PIC = max(PIC, 0.0)
	# aa-tRNAs: appear as k_T3j_f * EFT * tXX; negative flips T3 formation sign
	tA1 = max(tA1, 0.0); tA2 = max(tA2, 0.0); tC1 = max(tC1, 0.0); tD1 = max(tD1, 0.0)
	tE2 = max(tE2, 0.0); tF1 = max(tF1, 0.0); tG2 = max(tG2, 0.0); tG3 = max(tG3, 0.0)
	tH1 = max(tH1, 0.0); tI1 = max(tI1, 0.0); tK1 = max(tK1, 0.0)
	tL1 = max(tL1, 0.0); tL2 = max(tL2, 0.0); tL3 = max(tL3, 0.0); tL4 = max(tL4, 0.0); tL5 = max(tL5, 0.0)
	tM1 = max(tM1, 0.0); tN1 = max(tN1, 0.0)
	tP1 = max(tP1, 0.0); tP2 = max(tP2, 0.0); tP3 = max(tP3, 0.0)
	tQ1 = max(tQ1, 0.0); tQ2 = max(tQ2, 0.0)
	tR2 = max(tR2, 0.0); tR3 = max(tR3, 0.0); tR4 = max(tR4, 0.0); tR5 = max(tR5, 0.0)
	tS1 = max(tS1, 0.0); tS2 = max(tS2, 0.0); tS3 = max(tS3, 0.0); tS5 = max(tS5, 0.0)
	tT1 = max(tT1, 0.0); tT2 = max(tT2, 0.0); tT4 = max(tT4, 0.0)
	tV1 = max(tV1, 0.0); tV2 = max(tV2, 0.0); tW1 = max(tW1, 0.0); tY1 = max(tY1, 0.0)
	# T3 complexes: appear as k_T3j_d * a1t (dissociation); negative flips to production
	a1t = max(a1t, 0.0); a2t = max(a2t, 0.0); c1t = max(c1t, 0.0); d1t = max(d1t, 0.0)
	e2t = max(e2t, 0.0); f1t = max(f1t, 0.0); g2t = max(g2t, 0.0); g3t = max(g3t, 0.0)
	h1t = max(h1t, 0.0); i1t = max(i1t, 0.0); k1t = max(k1t, 0.0)
	l1t = max(l1t, 0.0); l2t = max(l2t, 0.0); l3t = max(l3t, 0.0); l4t = max(l4t, 0.0); l5t = max(l5t, 0.0)
	m1t = max(m1t, 0.0); n1t = max(n1t, 0.0)
	p1t = max(p1t, 0.0); p2t = max(p2t, 0.0); p3t = max(p3t, 0.0)
	q1t = max(q1t, 0.0); q2t = max(q2t, 0.0)
	r2t = max(r2t, 0.0); r3t = max(r3t, 0.0); r4t = max(r4t, 0.0); r5t = max(r5t, 0.0)
	s1t = max(s1t, 0.0); s2t = max(s2t, 0.0); s3t = max(s3t, 0.0); s5t = max(s5t, 0.0)
	t1t = max(t1t, 0.0); t2t = max(t2t, 0.0); t4t = max(t4t, 0.0)
	v1t = max(v1t, 0.0); v2t = max(v2t, 0.0); w1t = max(w1t, 0.0); y1t = max(y1t, 0.0)

	# Precompute inverse values (division is slow)
	inv_ATP = 1/ATP

	inv_Ala = 1/Ala 
	inv_Cys = 1/Cys
	inv_Asp = 1/Asp
	inv_Glu = 1/Glu
	inv_Phe = 1/Phe
	inv_Gly = 1/Gly
	inv_His = 1/His
	inv_Ile = 1/Ile
	inv_Lys = 1/Lys
	inv_Leu = 1/Leu
	inv_Met = 1/Met
	inv_Asn = 1/Asn
	inv_Pro = 1/Pro
	inv_Gln = 1/Gln
	inv_Arg = 1/Arg
	inv_Ser = 1/Ser
	inv_Thr = 1/Thr
	inv_Val = 1/val
	inv_Trp = 1/Trp
	inv_Tyr = 1/Tyr

	# Define RIB1
	idx_first_codon = length(tags.species_and_sequences.Tag)+1
	@views c_RIB	 = u[idx_first_codon:end];

	# ODEs ####################################################################################

	if RNA <= 1e-8
		q_RBS = 0.0
	else
		sum_RIB12 = sum(@view(c_RIB[1:12]))
		PTC_safe = max(PTC, 0.0)
		q_RBS = clamp(1.0 - (sum_RIB12 + S70 + PTC_safe) / (RNA + sum_RIB12), 0.0, 1.0)
	end


	# Initiation ###########################################
	# Standard 30S cascade (always used; spindler_75pct only affects initial conditions)

	V_IN_IF2D 		= k_IN_IF2D * S70 * tags.scale.V_IN_IF2D_scale

	V_IN_70SIC_max  = k_IN_70SIC * PIC * tags.scale.V_IN_70SIC_scale

	V_IN_70SIC 		= ((V_IN_70SIC_max)/  ( 1+(KM_mRNA/RNA+(KM_50S/S50)) ).* q_RBS)


	du[tags.idx.S30] = - k_on_30S_IF3*S30*IF3 + k_off_30SIF3*SC1 # 13
	du[tags.idx.S50] = - V_IN_70SIC # 14
	du[tags.idx.IF1] = - k_on_30SIF3_IF1*SC1*IF1    + k_off_30SIF3IF1*SC2 - k_on_30SIF3IF2_IF1*SC3*IF1 + k_off_30SIF3IF2IF1*SC4 + V_IN_IF2D # 15
	du[tags.idx.IF2] = - k_on_30SIF3_IF2*SC1*IF2    + k_off_30SIF3IF2*SC3 - k_on_30SIF3IF1_IF2*SC2*IF2 + k_off_30SIF3IF1IF2*SC4 + V_IN_IF2D# 16
	du[tags.idx.IF3] = - k_on_30S_IF3*S30*IF3 + k_off_30SIF3*SC1 + V_IN_IF2D # 17
	du[tags.idx.SC1] = + k_on_30S_IF3*S30*IF3 - k_off_30SIF3*SC1 + k_off_30SIF3IF1*SC2  - k_on_30SIF3_IF1*SC1*IF1 + k_off_30SIF3IF2*SC3    - k_on_30SIF3_IF2*SC1*IF2 #18
	du[tags.idx.SC2] = + k_on_30SIF3_IF1*SC1*IF1 - k_off_30SIF3IF1*SC2 + k_off_30SIF3IF1IF2*SC4    - k_on_30SIF3IF1_IF2*SC2*IF2 # 19
	du[tags.idx.SC3] = + k_on_30SIF3_IF2*SC1*IF2 - k_off_30SIF3IF2*SC3 + k_off_30SIF3IF2IF1*SC4    - k_on_30SIF3IF2_IF1*SC3*IF1 # 20
	du[tags.idx.SC4] = + k_off_30SIF123fMet*PIC       - k_on_30SIF123_fMet*SC4*b1t + k_on_30SIF3IF1_IF2*SC2*IF2 - k_off_30SIF3IF1IF2*SC4 + k_on_30SIF3IF2_IF1*SC3*IF1 - k_off_30SIF3IF2IF1*SC4 # 21
	du[tags.idx.PIC] =  + k_on_30SIF123_fMet*SC4*b1t - k_off_30SIF123fMet*PIC - V_IN_70SIC # 22 30S PIC
	du[tags.idx.b1t] = + k_off_30SIF123fMet*PIC - k_on_30SIF123_fMet*SC4*b1t + V_IN_IF2D # 35, db1t
	du[tags.idx.S70] = V_IN_70SIC - V_IN_IF2D # 171, dS70

	du[tags.idx.rna] = -V_IN_70SIC  # mRNA consumed by ribosome binding RBS; V_TL_EL[12] + production - degradation added in rhs_sequence.jl

	du[idx_first_codon] = V_IN_IF2D # 175, dRIB



	# tRNA charging - concentration of unloaded tRNA ###############################################
	# V_ARS[j] = k_cat[j] * c_ARS[j] * ARS_scale[j]  (Spindler: ARS_scale=20)
	V_ARS = k_cat .* tags.init_c.c_ARS .* tags.scale.ARS_scale

	du[tags.idx.ta1] = -(V_ARS[1] / (1 + KM_ARS_aa[1] *inv_Ala + KM_ARS_ATP[1]*inv_ATP + KM_ARS_tRNA[1]/ta1))	# 36
	du[tags.idx.ta2] = -(V_ARS[1] / (1 + KM_ARS_aa[1] *inv_Ala + KM_ARS_ATP[1]*inv_ATP + KM_ARS_tRNA[2]/ta2))	# 37
	du[tags.idx.tc1] = -(V_ARS[2] / (1 + KM_ARS_aa[2] *inv_Cys + KM_ARS_ATP[2]*inv_ATP + KM_ARS_tRNA[3]/tc1))	# 38
	du[tags.idx.td1] = -(V_ARS[3] / (1 + KM_ARS_aa[3] *inv_Asp + KM_ARS_ATP[3]*inv_ATP + KM_ARS_tRNA[4]/td1))	# 39
	du[tags.idx.te2] = -(V_ARS[4] / (1 + KM_ARS_aa[4] *inv_Glu + KM_ARS_ATP[4]*inv_ATP + KM_ARS_tRNA[5]/te2))	# 40
	du[tags.idx.tf1] = -(V_ARS[5] / (1 + KM_ARS_aa[5] *inv_Phe + KM_ARS_ATP[5]*inv_ATP + KM_ARS_tRNA[6]/tf1))	# 41
	du[tags.idx.tg2] = -(V_ARS[6] / (1 + KM_ARS_aa[6] *inv_Gly + KM_ARS_ATP[6]*inv_ATP + KM_ARS_tRNA[7]/tg2))	# 42
	du[tags.idx.tg3] = -(V_ARS[6] / (1 + KM_ARS_aa[6] *inv_Gly + KM_ARS_ATP[6]*inv_ATP + KM_ARS_tRNA[8]/tg3))	# 43
	du[tags.idx.th1] = -(V_ARS[7] / (1 + KM_ARS_aa[7] *inv_His + KM_ARS_ATP[7]*inv_ATP + KM_ARS_tRNA[9]/th1))	# 44
	du[tags.idx.ti1] = -(V_ARS[8] / (1 + KM_ARS_aa[8] *inv_Ile + KM_ARS_ATP[8]*inv_ATP + KM_ARS_tRNA[10]/ti1))		# 45
	du[tags.idx.tk1] = -(V_ARS[9] / (1 + KM_ARS_aa[9] *inv_Lys + KM_ARS_ATP[9]*inv_ATP + KM_ARS_tRNA[11]/tk1))		# 46
	du[tags.idx.t1l] = -(V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[12]/t1l))	# 47
	du[tags.idx.t2l] = -(V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[13]/t2l))	# 48
	du[tags.idx.t3l] = -(V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[14]/t3l))	# 49
	du[tags.idx.t4l] = -(V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[15]/t4l))	# 50
	du[tags.idx.t5l] = -(V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[16]/t5l))	# 51
	du[tags.idx.tm1] = -(V_ARS[11] / (1 + KM_ARS_aa[11] *inv_Met + KM_ARS_ATP[11]*inv_ATP + KM_ARS_tRNA[17]/tm1))	# 52
	du[tags.idx.tn1] = -(V_ARS[12] / (1 + KM_ARS_aa[12] *inv_Asn + KM_ARS_ATP[12]*inv_ATP + KM_ARS_tRNA[18]/tn1))	# 53
	du[tags.idx.tp1] = -(V_ARS[13] / (1 + KM_ARS_aa[13] *inv_Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[19]/tp1))	# 54
	du[tags.idx.tp2] = -(V_ARS[13] / (1 + KM_ARS_aa[13] *inv_Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[20]/tp2))	# 55
	du[tags.idx.tp3] = -(V_ARS[13] / (1 + KM_ARS_aa[13] *inv_Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[21]/tp3))	# 56
	du[tags.idx.tq1] = -(V_ARS[14] / (1 + KM_ARS_aa[14] *inv_Gln + KM_ARS_ATP[14]*inv_ATP + KM_ARS_tRNA[22]/tq1))	# 57
	du[tags.idx.tq2] = -(V_ARS[14] / (1 + KM_ARS_aa[14] *inv_Gln + KM_ARS_ATP[14]*inv_ATP + KM_ARS_tRNA[23]/tq2))	# 58
	du[tags.idx.tr2] = -(V_ARS[15] / (1 + KM_ARS_aa[15] *inv_Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[24]/tr2))	# 59
	du[tags.idx.tr3] = -(V_ARS[15] / (1 + KM_ARS_aa[15] *inv_Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[25]/tr3))	# 60
	du[tags.idx.tr4] = -(V_ARS[15] / (1 + KM_ARS_aa[15] *inv_Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[26]/tr4)) 	# 61
	du[tags.idx.tr5] = -(V_ARS[15] / (1 + KM_ARS_aa[15] *inv_Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[27]/tr5))	# 62
	du[tags.idx.ts1] = -(V_ARS[16] / (1 + KM_ARS_aa[16] *inv_Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[28]/ts1))	# 63
	du[tags.idx.ts2] = -(V_ARS[16] / (1 + KM_ARS_aa[16] *inv_Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[29]/ts2))	# 64
	du[tags.idx.ts3] = -(V_ARS[16] / (1 + KM_ARS_aa[16] *inv_Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[30]/ts3))	# 65
	du[tags.idx.ts5] = -(V_ARS[16] / (1 + KM_ARS_aa[16] *inv_Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[31]/ts5))	# 66
	du[tags.idx.tt1] = -(V_ARS[17] / (1 + KM_ARS_aa[17] *inv_Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[32]/tt1))	# 67
	du[tags.idx.tt2] = -(V_ARS[17] / (1 + KM_ARS_aa[17] *inv_Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[33]/tt2))	# 68
	du[tags.idx.tt4] = -(V_ARS[17] / (1 + KM_ARS_aa[17] *inv_Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[34]/tt4))	# 69
	du[tags.idx.tv1] = -(V_ARS[18] / (1 + KM_ARS_aa[18] *inv_Val + KM_ARS_ATP[18]*inv_ATP + KM_ARS_tRNA[35]/tv1))	# 70
	du[tags.idx.tv2] = -(V_ARS[18] / (1 + KM_ARS_aa[18] *inv_Val + KM_ARS_ATP[18]*inv_ATP + KM_ARS_tRNA[36]/tv2))	# 71
	du[tags.idx.tw1] = -(V_ARS[19] / (1 + KM_ARS_aa[19] *inv_Trp + KM_ARS_ATP[19]*inv_ATP + KM_ARS_tRNA[37]/tw1))	# 72
	du[tags.idx.ty1] = -(V_ARS[20] / (1 + KM_ARS_aa[20] *inv_Tyr + KM_ARS_ATP[20]*inv_ATP + KM_ARS_tRNA[38]/ty1))	# 73

	# tRNA charging - concentration of aa_tRNA ###############################################

	du[tags.idx.tA1]  = (V_ARS[1] / (1 + KM_ARS_aa[1] *inv_Ala + KM_ARS_ATP[1]*inv_ATP + KM_ARS_tRNA[1]/ta1))	 - (k_T3j_f * EFT * tA1 - k_T3j_d * a1t)	# 94
	du[tags.idx.tA2]  = (V_ARS[1] / (1 + KM_ARS_aa[1] *inv_Ala + KM_ARS_ATP[1]*inv_ATP + KM_ARS_tRNA[2]/ta2))	 - (k_T3j_f * EFT * tA2 - k_T3j_d * a2t)	# 95
	du[tags.idx.tC1]  = (V_ARS[2] / (1 + KM_ARS_aa[2] *inv_Cys + KM_ARS_ATP[2]*inv_ATP + KM_ARS_tRNA[3]/tc1))	 - (k_T3j_f * EFT * tC1 - k_T3j_d * c1t)	# 96
	du[tags.idx.tD1]  = (V_ARS[3] / (1 + KM_ARS_aa[3] *inv_Asp + KM_ARS_ATP[3]*inv_ATP + KM_ARS_tRNA[4]/td1))	 - (k_T3j_f * EFT * tD1 - k_T3j_d * d1t)	# 97
	du[tags.idx.tE2]  = (V_ARS[4] / (1 + KM_ARS_aa[4] *inv_Glu + KM_ARS_ATP[4]*inv_ATP + KM_ARS_tRNA[5]/te2))	 - (k_T3j_f * EFT * tE2 - k_T3j_d * e2t)	# 98
	du[tags.idx.tF1]  = (V_ARS[5] / (1 + KM_ARS_aa[5] *inv_Phe + KM_ARS_ATP[5]*inv_ATP + KM_ARS_tRNA[6]/tf1))	 - (k_T3j_f * EFT * tF1 - k_T3j_d * f1t)	# 99
	du[tags.idx.tG2] = (V_ARS[6] / (1 + KM_ARS_aa[6] *inv_Gly + KM_ARS_ATP[6]*inv_ATP + KM_ARS_tRNA[7]/tg2))	 - (k_T3j_f * EFT * tG2 - k_T3j_d * g2t)	# 100
	du[tags.idx.tG3] = (V_ARS[6] / (1 + KM_ARS_aa[6] *inv_Gly + KM_ARS_ATP[6]*inv_ATP + KM_ARS_tRNA[8]/tg3))	 - (k_T3j_f * EFT * tG3 - k_T3j_d * g3t)	# 101
	du[tags.idx.tH1] = (V_ARS[7] / (1 + KM_ARS_aa[7] *inv_His + KM_ARS_ATP[7]*inv_ATP + KM_ARS_tRNA[9]/th1))	 - (k_T3j_f * EFT * tH1 - k_T3j_d * h1t)	# 102
	du[tags.idx.tI1] = (V_ARS[8] / (1 + KM_ARS_aa[8] *inv_Ile + KM_ARS_ATP[8]*inv_ATP + KM_ARS_tRNA[10]/ti1))	 - (k_T3j_f * EFT * tI1 - k_T3j_d * i1t)	# 103	
	du[tags.idx.tK1] = (V_ARS[9] / (1 + KM_ARS_aa[9] *inv_Lys + KM_ARS_ATP[9]*inv_ATP + KM_ARS_tRNA[11]/tk1))	 - (k_T3j_f * EFT * tK1 - k_T3j_d * k1t)	# 104	
	du[tags.idx.tL1] = (V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[12]/t1l)) - (k_T3j_f * EFT * tL1 - k_T3j_d * l1t)	# 105
	du[tags.idx.tL2] = (V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[13]/t2l)) - (k_T3j_f * EFT * tL2 - k_T3j_d * l2t)	# 106
	du[tags.idx.tL3] = (V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[14]/t3l)) - (k_T3j_f * EFT * tL3 - k_T3j_d * l3t)	# 107
	du[tags.idx.tL4] = (V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[15]/t4l)) - (k_T3j_f * EFT * tL4 - k_T3j_d * l4t)	# 108
	du[tags.idx.tL5] = (V_ARS[10] / (1 + KM_ARS_aa[10] *inv_Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[16]/t5l)) - (k_T3j_f * EFT * tL5 - k_T3j_d * l5t)	# 109
	du[tags.idx.tM1] = (V_ARS[11] / (1 + KM_ARS_aa[11] *inv_Met + KM_ARS_ATP[11]*inv_ATP + KM_ARS_tRNA[17]/tm1)) - (k_T3j_f * EFT * tM1 - k_T3j_d * m1t)	# 110
	du[tags.idx.tN1] = (V_ARS[12] / (1 + KM_ARS_aa[12] *inv_Asn + KM_ARS_ATP[12]*inv_ATP + KM_ARS_tRNA[18]/tn1)) - (k_T3j_f * EFT * tN1 - k_T3j_d * n1t)	# 111
	du[tags.idx.tP1] = (V_ARS[13] / (1 + KM_ARS_aa[13] *inv_Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[19]/tp1)) - (k_T3j_f * EFT * tP1 - k_T3j_d * p1t)	# 112
	du[tags.idx.tP2] = (V_ARS[13] / (1 + KM_ARS_aa[13] *inv_Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[20]/tp2)) - (k_T3j_f * EFT * tP2 - k_T3j_d * p2t)	# 113
	du[tags.idx.tP3] = (V_ARS[13] / (1 + KM_ARS_aa[13] *inv_Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[21]/tp3)) - (k_T3j_f * EFT * tP3 - k_T3j_d * p3t)	# 114
	du[tags.idx.tQ1] = (V_ARS[14] / (1 + KM_ARS_aa[14] *inv_Gln + KM_ARS_ATP[14]*inv_ATP + KM_ARS_tRNA[22]/tq1)) - (k_T3j_f * EFT * tQ1 - k_T3j_d * q1t)	# 115
	du[tags.idx.tQ2] = (V_ARS[14] / (1 + KM_ARS_aa[14] *inv_Gln + KM_ARS_ATP[14]*inv_ATP + KM_ARS_tRNA[23]/tq2)) - (k_T3j_f * EFT * tQ2 - k_T3j_d * q2t)	# 116
	du[tags.idx.tR2] = (V_ARS[15] / (1 + KM_ARS_aa[15] *inv_Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[24]/tr2)) - (k_T3j_f * EFT * tR2 - k_T3j_d * r2t)	# 117
	du[tags.idx.tR3] = (V_ARS[15] / (1 + KM_ARS_aa[15] *inv_Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[25]/tr3)) - (k_T3j_f * EFT * tR3 - k_T3j_d * r3t)	# 118
	du[tags.idx.tR4] = (V_ARS[15] / (1 + KM_ARS_aa[15] *inv_Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[26]/tr4)) - (k_T3j_f * EFT * tR4 - k_T3j_d * r4t)	# 119
	du[tags.idx.tR5] = (V_ARS[15] / (1 + KM_ARS_aa[15] *inv_Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[27]/tr5)) - (k_T3j_f * EFT * tR5 - k_T3j_d * r5t)	# 120
	du[tags.idx.tS1] = (V_ARS[16] / (1 + KM_ARS_aa[16] *inv_Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[28]/ts1)) - (k_T3j_f * EFT * tS1 - k_T3j_d * s1t)	# 121
	du[tags.idx.tS2] = (V_ARS[16] / (1 + KM_ARS_aa[16] *inv_Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[29]/ts2)) - (k_T3j_f * EFT * tS2 - k_T3j_d * s2t)	# 122
	du[tags.idx.tS3] = (V_ARS[16] / (1 + KM_ARS_aa[16] *inv_Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[30]/ts3)) - (k_T3j_f * EFT * tS3 - k_T3j_d * s3t)	# 123
	du[tags.idx.tS5] = (V_ARS[16] / (1 + KM_ARS_aa[16] *inv_Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[31]/ts5)) - (k_T3j_f * EFT * tS5 - k_T3j_d * s5t)	# 124
	du[tags.idx.tT1] = (V_ARS[17] / (1 + KM_ARS_aa[17] *inv_Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[32]/tt1)) - (k_T3j_f * EFT * tT1 - k_T3j_d * t1t)	# 125
	du[tags.idx.tT2] = (V_ARS[17] / (1 + KM_ARS_aa[17] *inv_Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[33]/tt2)) - (k_T3j_f * EFT * tT2 - k_T3j_d * t2t)	# 126
	du[tags.idx.tT4] = (V_ARS[17] / (1 + KM_ARS_aa[17] *inv_Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[34]/tt4)) - (k_T3j_f * EFT * tT4 - k_T3j_d * t4t)	# 127
	du[tags.idx.tV1] = (V_ARS[18] / (1 + KM_ARS_aa[18] *inv_Val + KM_ARS_ATP[18]*inv_ATP + KM_ARS_tRNA[35]/tv1)) - (k_T3j_f * EFT * tV1 - k_T3j_d * v1t)	# 128
	du[tags.idx.tV2] = (V_ARS[18] / (1 + KM_ARS_aa[18] *inv_Val + KM_ARS_ATP[18]*inv_ATP + KM_ARS_tRNA[36]/tv2)) - (k_T3j_f * EFT * tV2 - k_T3j_d * v2t)	# 129
	du[tags.idx.tW1] = (V_ARS[19] / (1 + KM_ARS_aa[19] *inv_Trp + KM_ARS_ATP[19]*inv_ATP + KM_ARS_tRNA[37]/tw1)) - (k_T3j_f * EFT * tW1 - k_T3j_d * w1t)	# 130
	du[tags.idx.tY1] = (V_ARS[20] / (1 + KM_ARS_aa[20] *inv_Tyr + KM_ARS_ATP[20]*inv_ATP + KM_ARS_tRNA[38]/ty1)) - (k_T3j_f * EFT * tY1 - k_T3j_d * y1t)	# 131

	# tRNA charging - concentration of ternary complexes ###############################################

	du[tags.idx.a1t] = k_T3j_f * EFT * tA1 - k_T3j_d * a1t	# 132
	du[tags.idx.a2t] = k_T3j_f * EFT * tA2 - k_T3j_d * a2t	# 133
	du[tags.idx.c1t] = k_T3j_f * EFT * tC1 - k_T3j_d * c1t	# 134
	du[tags.idx.d1t] = k_T3j_f * EFT * tD1 - k_T3j_d * d1t	# 135
	du[tags.idx.e2t] = k_T3j_f * EFT * tE2 - k_T3j_d * e2t	# 136 
	du[tags.idx.f1t] = k_T3j_f * EFT * tF1 - k_T3j_d * f1t	# 137
	du[tags.idx.g2t] = k_T3j_f * EFT * tG2 - k_T3j_d * g2t	# 138
	du[tags.idx.g3t] = k_T3j_f * EFT * tG3 - k_T3j_d * g3t	# 139
	du[tags.idx.h1t] = k_T3j_f * EFT * tH1 - k_T3j_d * h1t	# 140
	du[tags.idx.i1t] = k_T3j_f * EFT * tI1 - k_T3j_d * i1t	# 141 
	du[tags.idx.k1t] = k_T3j_f * EFT * tK1 - k_T3j_d * k1t	# 142 
	du[tags.idx.l1t] = k_T3j_f * EFT * tL1 - k_T3j_d * l1t	# 143 
	du[tags.idx.l2t] = k_T3j_f * EFT * tL2 - k_T3j_d * l2t	# 144 
	du[tags.idx.l3t] = k_T3j_f * EFT * tL3 - k_T3j_d * l3t	# 145 
	du[tags.idx.l4t] = k_T3j_f * EFT * tL4 - k_T3j_d * l4t	# 146 
	du[tags.idx.l5t] = k_T3j_f * EFT * tL5 - k_T3j_d * l5t	# 147 
	du[tags.idx.m1t] = k_T3j_f * EFT * tM1 - k_T3j_d * m1t	# 148 
	du[tags.idx.n1t] = k_T3j_f * EFT * tN1 - k_T3j_d * n1t	# 149 
	du[tags.idx.p1t] = k_T3j_f * EFT * tP1 - k_T3j_d * p1t	# 150 
	du[tags.idx.p2t] = k_T3j_f * EFT * tP2 - k_T3j_d * p2t	# 151 
	du[tags.idx.p3t] = k_T3j_f * EFT * tP3 - k_T3j_d * p3t	# 152 
	du[tags.idx.q1t] = k_T3j_f * EFT * tQ1 - k_T3j_d * q1t	# 153
	du[tags.idx.q2t] = k_T3j_f * EFT * tQ2 - k_T3j_d * q2t	# 154 
	du[tags.idx.r2t] = k_T3j_f * EFT * tR2 - k_T3j_d * r2t	# 155 
	du[tags.idx.r3t] = k_T3j_f * EFT * tR3 - k_T3j_d * r3t	# 156 
	du[tags.idx.r4t] = k_T3j_f * EFT * tR4 - k_T3j_d * r4t	# 157 
	du[tags.idx.r5t] = k_T3j_f * EFT * tR5 - k_T3j_d * r5t	# 158 
	du[tags.idx.s1t] = k_T3j_f * EFT * tS1 - k_T3j_d * s1t	# 159 
	du[tags.idx.s2t] = k_T3j_f * EFT * tS2 - k_T3j_d * s2t	# 160
	du[tags.idx.s3t] = k_T3j_f * EFT * tS3 - k_T3j_d * s3t	# 161 
	du[tags.idx.s5t] = k_T3j_f * EFT * tS5 - k_T3j_d * s5t	# 162 
	du[tags.idx.t1t] = k_T3j_f * EFT * tT1 - k_T3j_d * t1t	# 163
	du[tags.idx.t2t] = k_T3j_f * EFT * tT2 - k_T3j_d * t2t	# 164 
	du[tags.idx.t4t] = k_T3j_f * EFT * tT4 - k_T3j_d * t4t	# 165 
	du[tags.idx.v1t] = k_T3j_f * EFT * tV1 - k_T3j_d * v1t	# 166 
	du[tags.idx.v2t] = k_T3j_f * EFT * tV2 - k_T3j_d * v2t	# 167
	du[tags.idx.w1t] = k_T3j_f * EFT * tW1 - k_T3j_d * w1t	# 168
	du[tags.idx.y1t] = k_T3j_f * EFT * tY1 - k_T3j_d * y1t	# 169

	if tags.scale.V_TC_scale == 0
		du[tags.idx.ta1:tags.idx.ty1] .= 0
		du[tags.idx.tA1:tags.idx.tY1] .= 0
		du[tags.idx.a1t:tags.idx.y1t] .= 0
	end

	# Elongation factors
	# EFG·GTP / EFG·GDP exchange (first-order in each pool).
	du[tags.idx.EGT] =  k_EFG_GTP_Ass * EFG * GTP - k_EFGGTP_Diss * EGT		# Rate of EFG_GTP # 25
	du[tags.idx.EGD] = -(k_EFGGDP_Diss * EGD - k_EFG_GDP_Ass * EFG * GDP)		# Rate of EFG_GDP # 26
	# EF-Tu (EFT): Michaelis-type rate (numerator / denominator) minus T3 consumption ∝ k_T3j_f*EFT*tRNA.
	# If EFT < 0: terms in the denominator that are linear in EFT become negative → denominator can shrink
	# and |dEFT/dt| explodes; T3 terms also flip sign. Keep EFT ≥ 0 physically (IC / clamp / reformulate).
	du[tags.idx.EFT] = (k_EFTsf * Est) * (EFD * GTP - (GDP * EFT/Keq_EFTu))/ (KM_GTP * EFD + KM_EFTuGDP * GTP + (k_EFTsf * Est * KM_EFTuGTP * GDP)/(k_EFTsr * Est * Keq_EFTu) + (k_EFTsf * Est * KM_GDP * EFT)/(k_EFTsr * Est * Keq_EFTu)+
			EFD * GTP + (k_EFTsf * Est * KM_EFTuGTP * EFD * GDP)/(k_EFTsr * Est * Keq_EFTu * Ki_EFTuGDP) + (k_EFTsf * Est * GDP * EFT) /(k_EFTsr * Est * Keq_EFTu) + 
			(KM_EFTuGDP * GTP * EFT)/ Ki_EFTuGTP) -
			+((k_T3j_f * EFT * tA1 - k_T3j_d * a1t),(k_T3j_f * EFT * tA2 - k_T3j_d * a2t),(k_T3j_f * EFT * tC1 - k_T3j_d * c1t),(k_T3j_f * EFT * tD1 - k_T3j_d * d1t),(k_T3j_f * EFT * tE2 - k_T3j_d * e2t),(k_T3j_f * EFT * tF1 - k_T3j_d * f1t),(k_T3j_f * EFT * tG2 - k_T3j_d * g2t),(k_T3j_f * EFT * tG3 - k_T3j_d * g3t),(k_T3j_f * EFT * tH1 - k_T3j_d * h1t),(k_T3j_f * EFT * tI1 - k_T3j_d * i1t),(k_T3j_f * EFT * tK1 - k_T3j_d * k1t),(k_T3j_f * EFT * tL1 - k_T3j_d * l1t),(k_T3j_f * EFT * tL2 - k_T3j_d * l2t),(k_T3j_f * EFT * tL3 - k_T3j_d * l3t),(k_T3j_f * EFT * tL4 - k_T3j_d * l4t),(k_T3j_f * EFT * tL5 - k_T3j_d * l5t),(k_T3j_f * EFT * tM1 - k_T3j_d * m1t),(k_T3j_f * EFT * tN1 - k_T3j_d * n1t),(k_T3j_f * EFT * tP1 - k_T3j_d * p1t),(k_T3j_f * EFT * tP2 - k_T3j_d * p2t),(k_T3j_f * EFT * tP3 - k_T3j_d * p3t),(k_T3j_f * EFT * tQ1 - k_T3j_d * q1t),(k_T3j_f * EFT * tQ2 - k_T3j_d * q2t),(k_T3j_f * EFT * tR2 - k_T3j_d * r2t),(k_T3j_f * EFT * tR3 - k_T3j_d * r3t),(k_T3j_f * EFT * tR4 - k_T3j_d * r4t),(k_T3j_f * EFT * tR5 - k_T3j_d * r5t),(k_T3j_f * EFT * tS1 - k_T3j_d * s1t),(k_T3j_f * EFT * tS2 - k_T3j_d * s2t),(k_T3j_f * EFT * tS3 - k_T3j_d * s3t),(k_T3j_f * EFT * tS5 - k_T3j_d * s5t),(k_T3j_f * EFT * tT1 - k_T3j_d * t1t),(k_T3j_f * EFT * tT2 - k_T3j_d * t2t),(k_T3j_f * EFT * tT4 - k_T3j_d * t4t),(k_T3j_f * EFT * tV1 - k_T3j_d * v1t),(k_T3j_f * EFT * tV2 - k_T3j_d * v2t),(k_T3j_f * EFT * tW1 - k_T3j_d * w1t),(k_T3j_f * EFT * tY1 - k_T3j_d * y1t)); # 27

	# steady-state of elongation factors

	T3  = zeros(Float64, tags.species_and_sequences.no_of_tRNAs);

	T3[1:end] = u[tags.idx.a1t:tags.idx.y1t]

	T3sum = sum(T3);

	# DAE algebraic constraints — MUST use raw state values (not clamped locals)
	# to preserve mass conservation.  Clamping is only for rate kinetics above.
	raw_EFG = u[tags.idx.EFG]
	raw_EGT = u[tags.idx.EGT]
	raw_EGD = u[tags.idx.EGD]
	raw_EFT = u[tags.idx.EFT]
	raw_EFD = u[tags.idx.EFD]
	du[tags.idx.EFG] = raw_EFG + raw_EGT + raw_EGD - c_EFG_t;   # EF-G conservation
	du[tags.idx.EFD] = raw_EFT + raw_EFD - Eut + T3sum;          # EF-Tu conservation

	# Spindler et al. 2025: EF regeneration reactions are static parameters.
	# Freeze EFG·GTP, EFG·GDP, EFTu·GTP at their preintegration equilibrium values.
	# T3 formation still uses the (constant) EFT value computed above; only the
	# differential equations for the EF pools themselves are zeroed.
	if tags.flags.static_EFs
		du[tags.idx.EGT] = 0.0
		du[tags.idx.EGD] = 0.0
		du[tags.idx.EFT] = 0.0
	end

	nothing
end