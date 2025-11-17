# like: MOD_met_kinetics.m 
# Start: 12.02.2025
# By: Cedric Meyer 
# Purpose: Contains the parameters for the metabolite module

function parameters(initial_concentrations, species_and_seqs)
    
    ## mod_met_kinetics

    # reaction kinetics: 70 S initiation complex
    k_IN_70SIC		=  0.1 * 60;# Wintermeyer 1983
    KM_mRNA			= .009;		# µM
    KM_50S			= .0119;	# µM 12*1e-3
    
    # reaction kinetics: IF2-dependent GTP hydrolysis    
    k_IN_IF2D		= .8 * 60;	# 1/min   second value in text: 1.5*60 1/min
    # k_IN_IF2D		= 1.5 * 60;	# 1/min   second value in text: 1.5*60 1/min # changed to get closer to the Niess 1.8 AA/rib/s translation rate
    
    ## kinetics from Milon et al. (2012)    
        
    # 30S IF3-binding
    k_on_30S_IF3 = 1.160 * 60; # 1/µM/min
    k_off_30SIF3 = 36 * 60; # 1/min
    
    # 30SIF3 IF1-binding
    k_on_30SIF3_IF1 = 20 * 60; # 1/µM/min
    k_off_30SIF3IF1 = 0.7 *60; # 1/min
    # 30SIF3 IF2-binding
    k_on_30SIF3_IF2 = 280 * 60; # 1/µM/min
    k_off_30SIF3IF2 = 12 * 60;  # 1/min
    # 30SIF3IF1 IF2-binding
    k_on_30SIF3IF1_IF2 = 220 * 60; # 1/µM/min
    k_off_30SIF3IF1IF2 = 15 * 60;  # 1/min
    # 30SIF3IF2 IF1-binding
    k_on_30SIF3IF2_IF1 = 12 * 60;  # 1/µM/min
    k_off_30SIF3IF2IF1 = 0.02 *60; # 1/min
    # 30SIF1IF2IF3 fMet_tRNA-binding
    k_on_30SIF123_fMet =   5 * 60;   # 1/µM/min
    k_off_30SIF123fMet = 1.5 * 60; # 1/min

    IN = txtl_types.initiation(k_IN_70SIC, KM_mRNA, KM_50S, k_IN_IF2D, k_on_30S_IF3, k_off_30SIF3, k_on_30SIF3_IF1, k_off_30SIF3IF1, k_on_30SIF3_IF2, k_off_30SIF3IF2, k_on_30SIF3IF1_IF2, k_off_30SIF3IF1IF2, k_on_30SIF3IF2_IF1, k_off_30SIF3IF2IF1, k_on_30SIF123_fMet, k_off_30SIF123fMet)
        
    ## elongation factor recycling kinetics
    
    k_EFGGDP_Diss	= 100*60;			# 1/min
    k_EFG_GDP_Ass	= 2.7e7*1e-6*60;	# 1/(µM*min) 
    k_EFG_GTP_Ass	= 1e7*1e-6*60;		# 1/(µM*min) 
    k_EFGGTP_Diss	= 400*60;			# 1/min
    
    k_EFTsf		= 30 * 60;	#1/min
    k_EFTsr		= 10 * 60;	#1/min
    Keq_EFTu	= .19;	# for Keq_EFTu Romero et al. use .19 as valje
    KM_GTP		= 50;	# µM
    KM_GDP		= 3;	# µM
    KM_EFTuGDP	= 2.5;	# µM
    KM_EFTuGTP	= 1;	# µM
    Ki_EFTuGDP	= 5.6;	# µM
    Ki_EFTuGTP	= 1;	# µM 
    c_EFG_t		= initial_concentrations.c_tot_EFG*scale.EF_scale # µM

    EF = txtl_types.EF_recycling(k_EFGGDP_Diss, k_EFG_GDP_Ass, k_EFG_GTP_Ass, k_EFGGTP_Diss, k_EFTsf, k_EFTsr, Keq_EFTu, KM_GTP, KM_GDP, KM_EFTuGDP, KM_EFTuGTP, Ki_EFTuGDP, Ki_EFTuGTP, c_EFG_t)
    
    #= tRNA charging   
    representative kinetics (low, mid, high), identified by k-means clustering: =#
    
    k_cat       = [40.2;  177; 40.2; 177; 40.2; 40.2; 97.8;  177; 97.8;  177;
                           40.2; 97.8; 40.2; 177; 97.8; 40.2; 40.2; 97.8; 97.8; 97.8]; # [1/min] 
    
    KM_ARS_tRNA = [1.13; 1.13; 1.13; 0.31; 0.31; 2.70; 0.31; 0.31; 0.31; 2.70;
                       2.70; 2.70; 2.70; 2.70; 2.70; 2.70; 2.70; 0.31; 2.70; 2.70;
                       2.70; 0.31; 0.31; 2.70; 2.70; 2.70; 2.70; 0.31; 0.31; 0.31;
                       0.31; 0.31; 0.31; 0.31; 0.31; 0.31; 0.31; 0.31];             # µM
    KM_ARS_aa   = [47.1;  7.8;  47.1;  47.1; 7.8; 47.1;  47.1; 7.8; 7.8; 7.8;
                       47.1; 47.1; 183.3; 183.3; 7.8; 47.1; 183.3; 7.8; 7.8; 7.8];  # µM
    KM_ARS_ATP  = [ 42; 230;  42; 230;  42; 42; 654;  42;  42; 230;
                       654;  42; 654;  42; 654;  42;  42; 42; 230;  42];            # µM

    tRNA = txtl_types.tRNA_charging(k_cat, KM_ARS_tRNA, KM_ARS_aa, KM_ARS_ATP)
    
    ## ternary complex formation
    
    # assumption: rate is equal for every aminoacid/aa-tRNA
    
    k_T3j_f		= 5 * 1e7 * 60 / 1e6;	# 1/(ÂµM*min) EFTuGTP + aatRNaj --> T3j   T3j= Ternary complex 
    k_T3j_d		= 1 * 60;				# 1/min T3j --> EFTuGTP + aatRNaj

    t3 = txtl_types.ternary_complex(k_T3j_f, k_T3j_d)

    q_RBS = 0
    
    p_met = txtl_types.parameters_metabolites(IN,EF,tRNA,t3,q_RBS)

    # tags.p_met = [q_RBS, k_IN_70SIC, KM_mRNA, KM_50S, k_IN_IF2D, k_on_30S_IF3, k_off_30SIF3, k_on_30SIF3_IF1, k_off_30SIF3IF1, k_on_30SIF3_IF2, k_off_30SIF3IF2, k_on_30SIF3IF1_IF2, k_off_30SIF3IF1IF2, k_on_30SIF3IF2_IF1, k_off_30SIF3IF2IF1, k_on_30SIF123_fMet, k_off_30SIF123fMet, k_EFGGDP_Diss, k_EFG_GDP_Ass, k_EFG_GTP_Ass, k_EFGGTP_Diss, k_EFTsf, k_EFTsr, Keq_EFTu, KM_GTP, KM_GDP, KM_EFTuGDP, KM_EFTuGTP, Ki_EFTuGDP, Ki_EFTuGTP, c_EFG_t, k_cat, KM_ARS_tRNA, KM_ARS_aa, KM_ARS_ATP, k_T3j_f, k_T3j_d]
   

    ## mod seq kinetics


    # Transcription kinetics

    ## mRNA production Spindler 
    mue = species_and_seqs.µ
    k_proc      = (-3.3407*mue^2 + 18.3977*mue + 29.4959) * 60;	# nuc/s --> nuc/min  # equation from Bremer&Dennis 39-55
    k_cat_Pol   = k_proc;				 # nuc/min
    C_RNP   	= initial_concentrations.c_TC[1]; 		 # average of all growth rates 6.54
    V_TCmax 	= k_cat_Pol * C_RNP;     # average of all growth rates #tags.c_TC(1); # µM/min, now directly calculated in mod_seq_rhs 
    # C_D			= 0.001;				 # nM --> 4.246e-3 µM/0.004246  /2.5e-3/one Gene --> 0.0010 µM

    ## mRNA production Nieß 
    # V_TCmax		= .0519    
    # C_D			= tags.C_D				 
    K_M_D		= 6.3e-3;				 # nM --> µM

    ## mRNA degradation 
    k_cat_E  = 11.2;             # 1/min k_cat of RNAseE
    C_RNE    = initial_concentrations.c_TC[2]     #0.86;
    # C_RNE    = 0.86     #0.86;
    # C_RNE    = 0.001     #0.86;
    V_degmax = k_cat_E * C_RNE; # average of all growth rates, now directly calculated in mod_seq_rhs 
    K_M_deg  = 5.6*10^-7 * 10^6; # µM 

    k_D = 0.289     # 1/min

    TX = txtl_types.transcription(µ, k_proc, k_cat_Pol, C_RNP, V_TCmax, K_M_D, k_cat_E, C_RNE, V_degmax, K_M_deg, k_D)


    # Translation 

    ## translation elongation
    k_TLEmaxj	= 24*60*scale.k_TLEmaxj_scale;	# 1/min Freierstroffer 1997
    KM_T3j		= .4;		# ÂµM original
    # KM_T3j		= .8;		# ÂµM
    KM_EFGGTP	= .22;		# ÂµM

    EL = txtl_types.elongation(k_TLEmaxj, KM_T3j, KM_EFGGTP)

    ## translation termination
    k_TE            = 24*60*scale.k_TLEmaxj_scale;     # 0.25*60 1/min Freierstroffer 1997
    KM_RK           = .0083;                                # µM aus Arnold
    KM_GTP          = 100;                                  # µM   # in Quelltext 20 #100
    KM_RF3GTP       = 0.076;                                # Koutmou # Marachi --> 0.02 µM

    k_RF3_diss_GTP  =  1.6516*60;                           # 1/min Shi et al. 2018
    k_RF3_diss_GDP  = 23.7903*60;                           # 1/min Shi et al. 2018
    k_RF3_diss_apo  =  1.1471*60;                           # 1/min Shi et al. 2018

    k_RF1_diss_GTP  =  0.38*60;                             # 1/min Shi et al. 2018
    k_RF1_diss_GDP  =  0.01*60;                             # 1/min Shi et al. 2018
    k_RF1_diss_apo  =  0.001*60;                            # 1/min Shi et al. 2018

    k_GTP_hydr      = 31.2;                                 # 0.52*60 1/min Shi et al. 2018
    k_RF1_ass       = 2064;                                 # 34.4*60 1/µM/min hetrick et al. 2009

    TER = txtl_types.termination(k_TE, KM_RK, KM_GTP, KM_RF3GTP, k_RF3_diss_GDP, k_RF1_diss_GTP, k_GTP_hydr, k_RF1_ass)

    # RF3 balance
    k_RF3_GDP_ass   = 1620;	                                # 27*60 1/µM/min Marraci, Rodnina 2016
    k_RF3_GTP_ass   = 72000;                                # 1.2*10^3*60 1/µM/min Marraci, Rodnina 2016
    k_RF3_GDP_diss  = 7.8;	                                # 0.13*60 1/min Marraci, Rodnina 2016
    k_RF3_GTP_diss  = 1440;	                                # 24*60 1/min Marraci, Rodnina 2016                

    RF3 =txtl_types.RF3_balance(k_RF3_GDP_ass, k_RF3_GTP_ass, k_RF3_GDP_diss, k_RF3_GTP_diss)

    ## post termination complex (PTC) recycling
    # kinetic constants from Borg et al. 2016
    A1                  = 0.041     # s
    A2                  = 1.07      # s * µM 
    A3                  = 0.126     # s  
    A4                  = 0.137     # s * µM 
    A5                  = 0.538     # s * µM ^2
    k_cat_K_M_G1        = 63*60     # 1 / µM / s
    k_cat_K_M_G2        = 8.1*60    # 1 / µM / s
    k_cat_G1            = 36*60     # 1 / s
    K_M_G1              = 0.56      # µM
    k_RRF               = 15*60     # 1 / µM / s
    q_RRF               = 58*60     # 1 / s
    K_RRF               = 3.9       # µM
    one_q_G2_k_split    = 1.1       
    k_max               = 25*60     # 1 / s

    PTC = txtl_types.PTC_recycling(A1, A2, A3, A4, A5, k_cat_K_M_G1, k_cat_K_M_G2, k_RRF, q_RRF, one_q_G2_k_split)

    # tags.p_seq = [V_degmax, k_cat_Pol, k_D, C_RNP, V_TCmax, K_M_D, C_RNE, k_cat_E, K_M_deg, k_TLEmaxj, KM_T3j, KM_EFGGTP, KM_RK, KM_RF3GTP, k_RF1_diss_GTP, k_RF3_diss_GDP, k_GTP_hydr, k_RF1_ass, k_RF3_GDP_ass, k_RF3_GTP_ass, k_RF3_GDP_diss, k_RF3_GTP_diss, A1, A2, A3, A4, A5, k_cat_K_M_G1, k_cat_K_M_G2, k_RRF, q_RRF, one_q_G2_k_split, k_TE];

    p_seq = txtl_types.parameters_sequence(TX, EL, TER, RF3, PTC)

    return p_met, p_seq
end