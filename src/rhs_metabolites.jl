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
	q_RBS, k_IN_70SIC, KM_mRNA, KM_50S, k_IN_IF2D, k_on_30S_IF3, k_off_30SIF3, k_on_30SIF3_IF1, k_off_30SIF3IF1, k_on_30SIF3_IF2, k_off_30SIF3IF2, k_on_30SIF3IF1_IF2, k_off_30SIF3IF1IF2, k_on_30SIF3IF2_IF1, k_off_30SIF3IF2IF1, k_on_30SIF123_fMet, k_off_30SIF123fMet, k_EFGGDP_Diss, k_EFG_GDP_Ass, k_EFG_GTP_Ass, k_EFGGTP_Diss, k_EFTsf, k_EFTsr, Keq_EFTu, KM_GTP, KM_GDP, KM_EFTuGDP, KM_EFTuGTP, Ki_EFTuGDP, Ki_EFTuGTP, c_EFG_t, k_cat, KM_ARS_tRNA, KM_ARS_aa, KM_ARS_ATP, k_T3j_f, k_T3j_d = tags.p_met

	# Unpack constant variables 
	ATP, GTP, UTP, CTP, ADP, GDP, UDP, CDP, AMP, GMP, UMP, CMP, Eut, Est, RFt, RRF, Ala, Cys, Asp, Glu, Phe, Gly, His, Ile, Lys, Leu, Met, Asn, Pro, Gln, Arg, Ser, Thr, Val, Trp, Tyr, C_D = tags.constant_vars

	# Unpack dynamic variables 
	S30, S50, IF1, IF2, IF3, SC1, SC2, SC3, SC4, PIC, EFG, EFD, EGT, EGD, EFT, RF3, RGT, RGD, b1t, ta1, ta2, tc1, td1, te2, tf1, tg2, tg3, th1, ti1, tk1, t1l, t2l, t3l, t4l, t5l, tm1, tn1, tp1, tp2, tp3, tq1, tq2, tr2, tr3, tr4, tr5, ts1, ts2, ts3, ts5, tt1, tt2, tt4, tv1, tv2, tw1, ty1, tA1, tA2, tC1, tD1, tE2, tF1, tG2, tG3, tH1, tI1, tK1, tL1, tL2, tL3, tL4, tL5, tM1, tN1, tP1, tP2, tP3, tQ1, tQ2, tR2, tR3, tR4, tR5, tS1, tS2, tS3, tS5, tT1, tT2, tT4, tV1, tV2, tW1, tY1, a1t, a2t, c1t, d1t, e2t, f1t, g2t, g3t, h1t, i1t, k1t, l1t, l2t, l3t, l4t, l5t, m1t, n1t, p1t, p2t, p3t, q1t, q2t, r2t, r3t, r4t, r5t, s1t, s2t, s3t, s5t, t1t, t2t, t4t, v1t, v2t, w1t, y1t, RNA, S70, rna, PRT, PTC, RIB1 = u 

	# Precompute inverse values (division is slow)
	inv_ATP = 1/ATP

	# Define RIB1
	idx_first_codon = length(tags.Tag)+1
	c_RIB	 = u[idx_first_codon:end];

	# ODEs ####################################################################################

	if RNA <= 1e-8
		q_RBS = 0
	else 
		q_RBS = 1 - sum(c_RIB[1:12].+S70) / (RNA + sum(c_RIB[1:12])) 
	end


	# Initiation ###########################################

	V_IN_IF2D 		= k_IN_IF2D * S70 * tags.scale.V_IN_IF2D_scale

	V_IN_70SIC_max  = k_IN_70SIC * PIC * tags.scale.V_IN_70SIC_scale

	V_IN_70SIC 		= ((V_IN_70SIC_max)/  ( 1+(KM_mRNA/RNA+(KM_50S/S50)) ).* q_RBS)  # Here is a big difference between Spindler and Nieß model, Investigate!!!!!!!!!!!!!!!


	du[tags.idx["S30"]] = - k_on_30S_IF3*S30*IF3 + k_off_30SIF3*SC1 # 13
	du[tags.idx["S50"]] = - V_IN_70SIC # 14
	du[tags.idx["IF1"]] = - k_on_30SIF3_IF1*SC1*IF1    + k_off_30SIF3IF1*SC2 - k_on_30SIF3IF2_IF1*SC3*IF1 + k_off_30SIF3IF2IF1*SC4 + V_IN_IF2D # 15
	du[tags.idx["IF2"]] = - k_on_30SIF3_IF2*SC1*IF2    + k_off_30SIF3IF2*SC3 - k_on_30SIF3IF1_IF2*SC2*IF2 + k_off_30SIF3IF1IF2*SC4 + V_IN_IF2D# 16
	du[tags.idx["IF3"]] = - k_on_30S_IF3*S30*IF3 + k_off_30SIF3*SC1 + V_IN_IF2D # 17
	du[tags.idx["SC1"]] = + k_on_30S_IF3*S30*IF3 - k_off_30SIF3*SC1 + k_off_30SIF3IF1*SC2  - k_on_30SIF3_IF1*SC1*IF1 + k_off_30SIF3IF2*SC3    - k_on_30SIF3_IF2*SC1*IF2 #18
	du[tags.idx["SC2"]] = + k_on_30SIF3_IF1*SC1*IF1 - k_off_30SIF3IF1*SC2 + k_off_30SIF3IF1IF2*SC4    - k_on_30SIF3IF1_IF2*SC2*IF2 # 19
	du[tags.idx["SC3"]] = + k_on_30SIF3_IF2*SC1*IF2 - k_off_30SIF3IF2*SC3 + k_off_30SIF3IF2IF1*SC4    - k_on_30SIF3IF2_IF1*SC3*IF1 # 20
	du[tags.idx["SC4"]] = + k_off_30SIF123fMet*PIC       - k_on_30SIF123_fMet*SC4*b1t + k_on_30SIF3IF1_IF2*SC2*IF2 - k_off_30SIF3IF1IF2*SC4 + k_on_30SIF3IF2_IF1*SC3*IF1 - k_off_30SIF3IF2IF1*SC4 # 21
	du[tags.idx["PIC"]] =  + k_on_30SIF123_fMet*SC4*b1t - k_off_30SIF123fMet*PIC - (((k_IN_70SIC * PIC)/ ( 1+(KM_mRNA/RNA+(KM_50S/S50)) )) .* q_RBS) # 22 30S PIC


	du[tags.idx["b1t"]] = + k_off_30SIF123fMet*PIC - k_on_30SIF123_fMet*SC4*b1t + V_IN_IF2D # 35, db1t
	du[tags.idx["RNA"]] = - V_IN_70SIC # 170, dRNA
	du[tags.idx["S70"]] = V_IN_70SIC - V_IN_IF2D # 171, dS70
	du[tags.idx["rna"]] = V_IN_70SIC # 172, drna

	# du[tags.idx.GTP] 	= - V_IN_IF2D # 2, dGTP - not  balanced

	# du[tags.idx.GDP] 	= V_IN_IF2D # 6, dGDP - not balanced

	du[idx_first_codon] = V_IN_IF2D # 175, dRIB



	# tRNA charging - concentration of unloaded tRNA ###############################################

	du[tags.idx["ta1"]] = -(k_cat[1] * tags.c_ARS[1] / (1 + KM_ARS_aa[1] /Ala + KM_ARS_ATP[1]*inv_ATP + KM_ARS_tRNA[1]/ta1))	# 36
	du[tags.idx["ta2"]] = -(k_cat[1] * tags.c_ARS[1] / (1 + KM_ARS_aa[1] /Ala + KM_ARS_ATP[1]*inv_ATP + KM_ARS_tRNA[2]/ta2))	# 37
	du[tags.idx["tc1"]] = -(k_cat[2] * tags.c_ARS[2] / (1 + KM_ARS_aa[2] /Cys + KM_ARS_ATP[2]*inv_ATP + KM_ARS_tRNA[3]/tc1))	# 38
	du[tags.idx["td1"]] = -(k_cat[3] * tags.c_ARS[3] / (1 + KM_ARS_aa[3] /Asp + KM_ARS_ATP[3]*inv_ATP + KM_ARS_tRNA[4]/td1))	# 39
	du[tags.idx["te2"]] = -(k_cat[4] * tags.c_ARS[4] / (1 + KM_ARS_aa[4] /Glu + KM_ARS_ATP[4]*inv_ATP + KM_ARS_tRNA[5]/te2))	# 40
	du[tags.idx["tf1"]] = -(k_cat[5] * tags.c_ARS[5] / (1 + KM_ARS_aa[5] /Phe + KM_ARS_ATP[5]*inv_ATP + KM_ARS_tRNA[6]/tf1))	# 41
	du[tags.idx["tg2"]] = -(k_cat[6] * tags.c_ARS[6] / (1 + KM_ARS_aa[6] /Gly + KM_ARS_ATP[6]*inv_ATP + KM_ARS_tRNA[7]/tg2))	# 42
	du[tags.idx["tg3"]] = -(k_cat[6] * tags.c_ARS[6] / (1 + KM_ARS_aa[6] /Gly + KM_ARS_ATP[6]*inv_ATP + KM_ARS_tRNA[8]/tg3))	# 43
	du[tags.idx["th1"]] = -(k_cat[7] * tags.c_ARS[7] / (1 + KM_ARS_aa[7] /His + KM_ARS_ATP[7]*inv_ATP + KM_ARS_tRNA[9]/th1))	# 44
	du[tags.idx["ti1"]] = -(k_cat[8] * tags.c_ARS[8] / (1 + KM_ARS_aa[8] /Ile + KM_ARS_ATP[8]*inv_ATP + KM_ARS_tRNA[10]/ti1))		# 45
	du[tags.idx["tk1"]] = -(k_cat[9] * tags.c_ARS[9] / (1 + KM_ARS_aa[9] /Lys + KM_ARS_ATP[9]*inv_ATP + KM_ARS_tRNA[11]/tk1))		# 46
	du[tags.idx["t1l"]] = -(k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[12]/t1l))	# 47
	du[tags.idx["t2l"]] = -(k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[13]/t2l))	# 48
	du[tags.idx["t3l"]] = -(k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[14]/t3l))	# 49
	du[tags.idx["t4l"]] = -(k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[15]/t4l))	# 50
	du[tags.idx["t5l"]] = -(k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[16]/t5l))	# 51
	du[tags.idx["tm1"]] = -(k_cat[11] * tags.c_ARS[11] / (1 + KM_ARS_aa[11] /Met + KM_ARS_ATP[11]*inv_ATP + KM_ARS_tRNA[17]/tm1))	# 52
	du[tags.idx["tn1"]] = -(k_cat[12] * tags.c_ARS[12] / (1 + KM_ARS_aa[12] /Asn + KM_ARS_ATP[12]*inv_ATP + KM_ARS_tRNA[18]/tn1))	# 53
	du[tags.idx["tp1"]] = -(k_cat[13] * tags.c_ARS[13] / (1 + KM_ARS_aa[13] /Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[19]/tp1))	# 54
	du[tags.idx["tp2"]] = -(k_cat[13] * tags.c_ARS[13] / (1 + KM_ARS_aa[13] /Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[20]/tp2))	# 55
	du[tags.idx["tp3"]] = -(k_cat[13] * tags.c_ARS[13] / (1 + KM_ARS_aa[13] /Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[21]/tp3))	# 56
	du[tags.idx["tq1"]] = -(k_cat[14] * tags.c_ARS[14] / (1 + KM_ARS_aa[14] /Gln + KM_ARS_ATP[14]*inv_ATP + KM_ARS_tRNA[22]/tq1))	# 57
	du[tags.idx["tq2"]] = -(k_cat[14] * tags.c_ARS[14] / (1 + KM_ARS_aa[14] /Gln + KM_ARS_ATP[14]*inv_ATP + KM_ARS_tRNA[23]/tq2))	# 58
	du[tags.idx["tr2"]] = -(k_cat[15] * tags.c_ARS[15] / (1 + KM_ARS_aa[15] /Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[24]/tr2))	# 59
	du[tags.idx["tr3"]] = -(k_cat[15] * tags.c_ARS[15] / (1 + KM_ARS_aa[15] /Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[25]/tr3))	# 60
	du[tags.idx["tr4"]] = -(k_cat[15] * tags.c_ARS[15] / (1 + KM_ARS_aa[15] /Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[26]/tr4)) 	# 61
	du[tags.idx["tr5"]] = -(k_cat[15] * tags.c_ARS[15] / (1 + KM_ARS_aa[15] /Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[27]/tr5))	# 62
	du[tags.idx["ts1"]] = -(k_cat[16] * tags.c_ARS[16] / (1 + KM_ARS_aa[16] /Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[28]/ts1))	# 63
	du[tags.idx["ts2"]] = -(k_cat[16] * tags.c_ARS[16] / (1 + KM_ARS_aa[16] /Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[29]/ts2))	# 64
	du[tags.idx["ts3"]] = -(k_cat[16] * tags.c_ARS[16] / (1 + KM_ARS_aa[16] /Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[30]/ts3))	# 65
	du[tags.idx["ts5"]] = -(k_cat[16] * tags.c_ARS[16] / (1 + KM_ARS_aa[16] /Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[31]/ts5))	# 66
	du[tags.idx["tt1"]] = -(k_cat[17] * tags.c_ARS[17] / (1 + KM_ARS_aa[17] /Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[32]/tt1))	# 67
	du[tags.idx["tt2"]] = -(k_cat[17] * tags.c_ARS[17] / (1 + KM_ARS_aa[17] /Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[33]/tt2))	# 68
	du[tags.idx["tt4"]] = -(k_cat[17] * tags.c_ARS[17] / (1 + KM_ARS_aa[17] /Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[34]/tt4))	# 69
	du[tags.idx["tv1"]] = -(k_cat[18] * tags.c_ARS[18] / (1 + KM_ARS_aa[18] /Val + KM_ARS_ATP[18]*inv_ATP + KM_ARS_tRNA[35]/tv1))	# 70
	du[tags.idx["tv2"]] = -(k_cat[18] * tags.c_ARS[18] / (1 + KM_ARS_aa[18] /Val + KM_ARS_ATP[18]*inv_ATP + KM_ARS_tRNA[36]/tv2))	# 71
	du[tags.idx["tw1"]] = -(k_cat[19] * tags.c_ARS[19] / (1 + KM_ARS_aa[19] /Trp + KM_ARS_ATP[19]*inv_ATP + KM_ARS_tRNA[37]/tw1))	# 72
	du[tags.idx["ty1"]] = -(k_cat[20] * tags.c_ARS[20] / (1 + KM_ARS_aa[20] /Tyr + KM_ARS_ATP[20]*inv_ATP + KM_ARS_tRNA[38]/ty1))	# 73

	# tRNA charging - concentration of aa_tRNA ###############################################

	du[tags.idx["tA1"]]  = (k_cat[1] * tags.c_ARS[1] / (1 + KM_ARS_aa[1] /Ala + KM_ARS_ATP[1]*inv_ATP + KM_ARS_tRNA[1]/ta1))	 - (k_T3j_f * EFT * tA1 - k_T3j_d * a1t)	# 94
	du[tags.idx["tA2"]]  = (k_cat[1] * tags.c_ARS[1] / (1 + KM_ARS_aa[1] /Ala + KM_ARS_ATP[1]*inv_ATP + KM_ARS_tRNA[2]/ta2))	 - (k_T3j_f * EFT * tA2 - k_T3j_d * a2t)	# 95
	du[tags.idx["tC1"]]  = (k_cat[2] * tags.c_ARS[2] / (1 + KM_ARS_aa[2] /Cys + KM_ARS_ATP[2]*inv_ATP + KM_ARS_tRNA[3]/tc1))	 - (k_T3j_f * EFT * tC1 - k_T3j_d * c1t)	# 96
	du[tags.idx["tD1"]]  = (k_cat[3] * tags.c_ARS[3] / (1 + KM_ARS_aa[3] /Asp + KM_ARS_ATP[3]*inv_ATP + KM_ARS_tRNA[4]/td1))	 - (k_T3j_f * EFT * tD1 - k_T3j_d * d1t)	# 97
	du[tags.idx["tE2"]]  = (k_cat[4] * tags.c_ARS[4] / (1 + KM_ARS_aa[4] /Glu + KM_ARS_ATP[4]*inv_ATP + KM_ARS_tRNA[5]/te2))	 - (k_T3j_f * EFT * tE2 - k_T3j_d * e2t)	# 98
	du[tags.idx["tF1"]]  = (k_cat[5] * tags.c_ARS[5] / (1 + KM_ARS_aa[5] /Phe + KM_ARS_ATP[5]*inv_ATP + KM_ARS_tRNA[6]/tf1))	 - (k_T3j_f * EFT * tF1 - k_T3j_d * f1t)	# 99
	du[tags.idx["tG2"]] = (k_cat[6] * tags.c_ARS[6] / (1 + KM_ARS_aa[6] /Gly + KM_ARS_ATP[6]*inv_ATP + KM_ARS_tRNA[7]/tg2))	 - (k_T3j_f * EFT * tG2 - k_T3j_d * g2t)	# 100
	du[tags.idx["tG3"]] = (k_cat[6] * tags.c_ARS[6] / (1 + KM_ARS_aa[6] /Gly + KM_ARS_ATP[6]*inv_ATP + KM_ARS_tRNA[8]/tg3))	 - (k_T3j_f * EFT * tG3 - k_T3j_d * g3t)	# 101
	du[tags.idx["tH1"]] = (k_cat[7] * tags.c_ARS[7] / (1 + KM_ARS_aa[7] /His + KM_ARS_ATP[7]*inv_ATP + KM_ARS_tRNA[9]/th1))	 - (k_T3j_f * EFT * tH1 - k_T3j_d * h1t)	# 102
	du[tags.idx["tI1"]] = (k_cat[8] * tags.c_ARS[8] / (1 + KM_ARS_aa[8] /Ile + KM_ARS_ATP[8]*inv_ATP + KM_ARS_tRNA[10]/ti1))	 - (k_T3j_f * EFT * tI1 - k_T3j_d * i1t)	# 103	
	du[tags.idx["tK1"]] = (k_cat[9] * tags.c_ARS[9] / (1 + KM_ARS_aa[9] /Lys + KM_ARS_ATP[9]*inv_ATP + KM_ARS_tRNA[11]/tk1))	 - (k_T3j_f * EFT * tK1 - k_T3j_d * k1t)	# 104	
	du[tags.idx["tL1"]] = (k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[12]/t1l)) - (k_T3j_f * EFT * tL1 - k_T3j_d * l1t)	# 105
	du[tags.idx["tL2"]] = (k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[13]/t2l)) - (k_T3j_f * EFT * tL2 - k_T3j_d * l2t)	# 106
	du[tags.idx["tL3"]] = (k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[14]/t3l)) - (k_T3j_f * EFT * tL3 - k_T3j_d * l3t)	# 107
	du[tags.idx["tL4"]] = (k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[15]/t4l)) - (k_T3j_f * EFT * tL4 - k_T3j_d * l4t)	# 108
	du[tags.idx["tL5"]] = (k_cat[10] * tags.c_ARS[10] / (1 + KM_ARS_aa[10] /Leu + KM_ARS_ATP[10]*inv_ATP + KM_ARS_tRNA[16]/t5l)) - (k_T3j_f * EFT * tL5 - k_T3j_d * l5t)	# 109
	du[tags.idx["tM1"]] = (k_cat[11] * tags.c_ARS[11] / (1 + KM_ARS_aa[11] /Met + KM_ARS_ATP[11]*inv_ATP + KM_ARS_tRNA[17]/tm1)) - (k_T3j_f * EFT * tM1 - k_T3j_d * m1t)	# 110
	du[tags.idx["tN1"]] = (k_cat[12] * tags.c_ARS[12] / (1 + KM_ARS_aa[12] /Asn + KM_ARS_ATP[12]*inv_ATP + KM_ARS_tRNA[18]/tn1)) - (k_T3j_f * EFT * tN1 - k_T3j_d * n1t)	# 111
	du[tags.idx["tP1"]] = (k_cat[13] * tags.c_ARS[13] / (1 + KM_ARS_aa[13] /Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[19]/tp1)) - (k_T3j_f * EFT * tP1 - k_T3j_d * p1t)	# 112
	du[tags.idx["tP2"]] = (k_cat[13] * tags.c_ARS[13] / (1 + KM_ARS_aa[13] /Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[20]/tp2)) - (k_T3j_f * EFT * tP2 - k_T3j_d * p2t)	# 113
	du[tags.idx["tP3"]] = (k_cat[13] * tags.c_ARS[13] / (1 + KM_ARS_aa[13] /Pro + KM_ARS_ATP[13]*inv_ATP + KM_ARS_tRNA[21]/tp3)) - (k_T3j_f * EFT * tP3 - k_T3j_d * p3t)	# 114
	du[tags.idx["tQ1"]] = (k_cat[14] * tags.c_ARS[14] / (1 + KM_ARS_aa[14] /Gln + KM_ARS_ATP[14]*inv_ATP + KM_ARS_tRNA[22]/tq1)) - (k_T3j_f * EFT * tQ1 - k_T3j_d * q1t)	# 115
	du[tags.idx["tQ2"]] = (k_cat[14] * tags.c_ARS[14] / (1 + KM_ARS_aa[14] /Gln + KM_ARS_ATP[14]*inv_ATP + KM_ARS_tRNA[23]/tq2)) - (k_T3j_f * EFT * tQ2 - k_T3j_d * q2t)	# 116
	du[tags.idx["tR2"]] = (k_cat[15] * tags.c_ARS[15] / (1 + KM_ARS_aa[15] /Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[24]/tr2)) - (k_T3j_f * EFT * tR2 - k_T3j_d * r2t)	# 117
	du[tags.idx["tR3"]] = (k_cat[15] * tags.c_ARS[15] / (1 + KM_ARS_aa[15] /Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[25]/tr3)) - (k_T3j_f * EFT * tR3 - k_T3j_d * r3t)	# 118
	du[tags.idx["tR4"]] = (k_cat[15] * tags.c_ARS[15] / (1 + KM_ARS_aa[15] /Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[26]/tr4)) - (k_T3j_f * EFT * tR4 - k_T3j_d * r4t)	# 119
	du[tags.idx["tR5"]] = (k_cat[15] * tags.c_ARS[15] / (1 + KM_ARS_aa[15] /Arg + KM_ARS_ATP[15]*inv_ATP + KM_ARS_tRNA[27]/tr5)) - (k_T3j_f * EFT * tR5 - k_T3j_d * r5t)	# 120
	du[tags.idx["tS1"]] = (k_cat[16] * tags.c_ARS[16] / (1 + KM_ARS_aa[16] /Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[28]/ts1)) - (k_T3j_f * EFT * tS1 - k_T3j_d * s1t)	# 121
	du[tags.idx["tS2"]] = (k_cat[16] * tags.c_ARS[16] / (1 + KM_ARS_aa[16] /Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[29]/ts2)) - (k_T3j_f * EFT * tS2 - k_T3j_d * s2t)	# 122
	du[tags.idx["tS3"]] = (k_cat[16] * tags.c_ARS[16] / (1 + KM_ARS_aa[16] /Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[30]/ts3)) - (k_T3j_f * EFT * tS3 - k_T3j_d * s3t)	# 123
	du[tags.idx["tS5"]] = (k_cat[16] * tags.c_ARS[16] / (1 + KM_ARS_aa[16] /Ser + KM_ARS_ATP[16]*inv_ATP + KM_ARS_tRNA[31]/ts5)) - (k_T3j_f * EFT * tS5 - k_T3j_d * s5t)	# 124
	du[tags.idx["tT1"]] = (k_cat[17] * tags.c_ARS[17] / (1 + KM_ARS_aa[17] /Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[32]/tt1)) - (k_T3j_f * EFT * tT1 - k_T3j_d * t1t)	# 125
	du[tags.idx["tT2"]] = (k_cat[17] * tags.c_ARS[17] / (1 + KM_ARS_aa[17] /Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[33]/tt2)) - (k_T3j_f * EFT * tT2 - k_T3j_d * t2t)	# 126
	du[tags.idx["tT4"]] = (k_cat[17] * tags.c_ARS[17] / (1 + KM_ARS_aa[17] /Thr + KM_ARS_ATP[17]*inv_ATP + KM_ARS_tRNA[34]/tt4)) - (k_T3j_f * EFT * tT4 - k_T3j_d * t4t)	# 127
	du[tags.idx["tV1"]] = (k_cat[18] * tags.c_ARS[18] / (1 + KM_ARS_aa[18] /Val + KM_ARS_ATP[18]*inv_ATP + KM_ARS_tRNA[35]/tv1)) - (k_T3j_f * EFT * tV1 - k_T3j_d * v1t)	# 128
	du[tags.idx["tV2"]] = (k_cat[18] * tags.c_ARS[18] / (1 + KM_ARS_aa[18] /Val + KM_ARS_ATP[18]*inv_ATP + KM_ARS_tRNA[36]/tv2)) - (k_T3j_f * EFT * tV2 - k_T3j_d * v2t)	# 129
	du[tags.idx["tW1"]] = (k_cat[19] * tags.c_ARS[19] / (1 + KM_ARS_aa[19] /Trp + KM_ARS_ATP[19]*inv_ATP + KM_ARS_tRNA[37]/tw1)) - (k_T3j_f * EFT * tW1 - k_T3j_d * w1t)	# 130
	du[tags.idx["tY1"]] = (k_cat[20] * tags.c_ARS[20] / (1 + KM_ARS_aa[20] /Tyr + KM_ARS_ATP[20]*inv_ATP + KM_ARS_tRNA[38]/ty1)) - (k_T3j_f * EFT * tY1 - k_T3j_d * y1t)	# 131

	# tRNA charging - concentration of ternary complexes ###############################################

	du[tags.idx["a1t"]] = k_T3j_f * EFT * tA1 - k_T3j_d * a1t	# 132
	du[tags.idx["a2t"]] = k_T3j_f * EFT * tA2 - k_T3j_d * a2t	# 133
	du[tags.idx["c1t"]] = k_T3j_f * EFT * tC1 - k_T3j_d * c1t	# 134
	du[tags.idx["d1t"]] = k_T3j_f * EFT * tD1 - k_T3j_d * d1t	# 135
	du[tags.idx["e2t"]] = k_T3j_f * EFT * tE2 - k_T3j_d * e2t	# 136 
	du[tags.idx["f1t"]] = k_T3j_f * EFT * tF1 - k_T3j_d * f1t	# 137
	du[tags.idx["g2t"]] = k_T3j_f * EFT * tG2 - k_T3j_d * g2t	# 138
	du[tags.idx["g3t"]] = k_T3j_f * EFT * tG3 - k_T3j_d * g3t	# 139
	du[tags.idx["h1t"]] = k_T3j_f * EFT * tH1 - k_T3j_d * h1t	# 140
	du[tags.idx["i1t"]] = k_T3j_f * EFT * tI1 - k_T3j_d * i1t	# 141 
	du[tags.idx["k1t"]] = k_T3j_f * EFT * tK1 - k_T3j_d * k1t	# 142 
	du[tags.idx["l1t"]] = k_T3j_f * EFT * tL1 - k_T3j_d * l1t	# 143 
	du[tags.idx["l2t"]] = k_T3j_f * EFT * tL2 - k_T3j_d * l2t	# 144 
	du[tags.idx["l3t"]] = k_T3j_f * EFT * tL3 - k_T3j_d * l3t	# 145 
	du[tags.idx["l4t"]] = k_T3j_f * EFT * tL4 - k_T3j_d * l4t	# 146 
	du[tags.idx["l5t"]] = k_T3j_f * EFT * tL5 - k_T3j_d * l5t	# 147 
	du[tags.idx["m1t"]] = k_T3j_f * EFT * tM1 - k_T3j_d * m1t	# 148 
	du[tags.idx["n1t"]] = k_T3j_f * EFT * tN1 - k_T3j_d * n1t	# 149 
	du[tags.idx["p1t"]] = k_T3j_f * EFT * tP1 - k_T3j_d * p1t	# 150 
	du[tags.idx["p2t"]] = k_T3j_f * EFT * tP2 - k_T3j_d * p2t	# 151 
	du[tags.idx["p3t"]] = k_T3j_f * EFT * tP3 - k_T3j_d * p3t	# 152 
	du[tags.idx["q1t"]] = k_T3j_f * EFT * tQ1 - k_T3j_d * q1t	# 153
	du[tags.idx["q2t"]] = k_T3j_f * EFT * tQ2 - k_T3j_d * q2t	# 154 
	du[tags.idx["r2t"]] = k_T3j_f * EFT * tR2 - k_T3j_d * r2t	# 155 
	du[tags.idx["r3t"]] = k_T3j_f * EFT * tR3 - k_T3j_d * r3t	# 156 
	du[tags.idx["r4t"]] = k_T3j_f * EFT * tR4 - k_T3j_d * r4t	# 157 
	du[tags.idx["r5t"]] = k_T3j_f * EFT * tR5 - k_T3j_d * r5t	# 158 
	du[tags.idx["s1t"]] = k_T3j_f * EFT * tS1 - k_T3j_d * s1t	# 159 
	du[tags.idx["s2t"]] = k_T3j_f * EFT * tS2 - k_T3j_d * s2t	# 160
	du[tags.idx["s3t"]] = k_T3j_f * EFT * tS3 - k_T3j_d * s3t	# 161 
	du[tags.idx["s5t"]] = k_T3j_f * EFT * tS5 - k_T3j_d * s5t	# 162 
	du[tags.idx["t1t"]] = k_T3j_f * EFT * tT1 - k_T3j_d * t1t	# 163
	du[tags.idx["t2t"]] = k_T3j_f * EFT * tT2 - k_T3j_d * t2t	# 164 
	du[tags.idx["t4t"]] = k_T3j_f * EFT * tT4 - k_T3j_d * t4t	# 165 
	du[tags.idx["v1t"]] = k_T3j_f * EFT * tV1 - k_T3j_d * v1t	# 166 
	du[tags.idx["v2t"]] = k_T3j_f * EFT * tV2 - k_T3j_d * v2t	# 167
	du[tags.idx["w1t"]] = k_T3j_f * EFT * tW1 - k_T3j_d * w1t	# 168
	du[tags.idx["y1t"]] = k_T3j_f * EFT * tY1 - k_T3j_d * y1t	# 169

	if tags.scale.V_TC_scale == 0
		du[tags.idx["ta1"]:tags.idx["ty1"]] .= 0
		du[tags.idx["tA1"]:tags.idx["tY1"]] .= 0
		du[tags.idx["a1t"]:tags.idx["y1t"]] .= 0
	end

	# Elongation factors 
	du[tags.idx["EGT"]] =  k_EFG_GTP_Ass * EFG * GTP - k_EFGGTP_Diss * EGT		# Rate of EFG_GTP # 25
	du[tags.idx["EGD"]] = -(k_EFGGDP_Diss * EGD - k_EFG_GDP_Ass * EFG * GDP)		# Rate of EFG_GDP # 26
	du[tags.idx["EFT"]] = (k_EFTsf * Est) * (EFD * GTP - (GDP * EFT/Keq_EFTu))/ (KM_GTP * EFD + KM_EFTuGDP * GTP + (k_EFTsf * Est * KM_EFTuGTP * GDP)/(k_EFTsr * Est * Keq_EFTu) + (k_EFTsf * Est * KM_GDP * EFT)/(k_EFTsr * Est * Keq_EFTu)+
			EFD * GTP + (k_EFTsf * Est * KM_EFTuGTP * EFD * GDP)/(k_EFTsr * Est * Keq_EFTu * Ki_EFTuGDP) + (k_EFTsf * Est * GDP * EFT) /(k_EFTsr * Est * Keq_EFTu) + 
			(KM_EFTuGDP * GTP * EFT)/ Ki_EFTuGTP) -
			+((k_T3j_f * EFT * tA1 - k_T3j_d * a1t),(k_T3j_f * EFT * tA2 - k_T3j_d * a2t),(k_T3j_f * EFT * tC1 - k_T3j_d * c1t),(k_T3j_f * EFT * tD1 - k_T3j_d * d1t),(k_T3j_f * EFT * tE2 - k_T3j_d * e2t),(k_T3j_f * EFT * tF1 - k_T3j_d * f1t),(k_T3j_f * EFT * tG2 - k_T3j_d * g2t),(k_T3j_f * EFT * tG3 - k_T3j_d * g3t),(k_T3j_f * EFT * tH1 - k_T3j_d * h1t),(k_T3j_f * EFT * tI1 - k_T3j_d * i1t),(k_T3j_f * EFT * tK1 - k_T3j_d * k1t),(k_T3j_f * EFT * tL1 - k_T3j_d * l1t),(k_T3j_f * EFT * tL2 - k_T3j_d * l2t),(k_T3j_f * EFT * tL3 - k_T3j_d * l3t),(k_T3j_f * EFT * tL4 - k_T3j_d * l4t),(k_T3j_f * EFT * tL5 - k_T3j_d * l5t),(k_T3j_f * EFT * tM1 - k_T3j_d * m1t),(k_T3j_f * EFT * tN1 - k_T3j_d * n1t),(k_T3j_f * EFT * tP1 - k_T3j_d * p1t),(k_T3j_f * EFT * tP2 - k_T3j_d * p2t),(k_T3j_f * EFT * tP3 - k_T3j_d * p3t),(k_T3j_f * EFT * tQ1 - k_T3j_d * q1t),(k_T3j_f * EFT * tQ2 - k_T3j_d * q2t),(k_T3j_f * EFT * tR2 - k_T3j_d * r2t),(k_T3j_f * EFT * tR3 - k_T3j_d * r3t),(k_T3j_f * EFT * tR4 - k_T3j_d * r4t),(k_T3j_f * EFT * tR5 - k_T3j_d * r5t),(k_T3j_f * EFT * tS1 - k_T3j_d * s1t),(k_T3j_f * EFT * tS2 - k_T3j_d * s2t),(k_T3j_f * EFT * tS3 - k_T3j_d * s3t),(k_T3j_f * EFT * tS5 - k_T3j_d * s5t),(k_T3j_f * EFT * tT1 - k_T3j_d * t1t),(k_T3j_f * EFT * tT2 - k_T3j_d * t2t),(k_T3j_f * EFT * tT4 - k_T3j_d * t4t),(k_T3j_f * EFT * tV1 - k_T3j_d * v1t),(k_T3j_f * EFT * tV2 - k_T3j_d * v2t),(k_T3j_f * EFT * tW1 - k_T3j_d * w1t),(k_T3j_f * EFT * tY1 - k_T3j_d * y1t)); # 27

	# steady-state of elongation factors

	T3  = zeros(Float64, tags.no_of_tRNAs);

	T3[1:end] = u[tags.idx["a1t"]:tags.idx["y1t"]]

	T3sum = sum(T3);

	# DAE system
	du[tags.idx["EFG"]] = EFG + EGT + EGD - c_EFG_t;   # EFG
	du[tags.idx["EFD"]] = EFT + EFD - Eut + T3sum;     # EFTu_GDP
	
	nothing
end