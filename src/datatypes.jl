# Written by: Arnav Kumar Upadhyay, Cedric Meyer
# Last updated: 05.11.2025

module txtl_types

using BioSequences

export SequenceData, scale_struct, data_struct, init_scale, init_data, init_indices

struct SequenceData
    DNA::LongDNA{2}
    RNA::LongRNA{2}
    AA1::LongAA
    Codon::Vector{LongRNA{2}}
end

struct scale
    V_TC_scale       :: Float64
    k_TLEmaxj_scale  :: Float64
    V_IN_70SIC_scale :: Float64
    V_IN_IF2D_scale  :: Float64
    Rib_scale        :: Float64
    IF_scale         :: Float64
    EF_scale         :: Float64
    tRNA_scale       :: Float64
    ARS_scale        :: Vector{Float64}
    RibRecycling     :: Float64
    V_mDeg_scale     :: Float64
end

function init_scale()
    return scale(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 .* ones(20), 0.0, 0.0)
end


mutable struct indices_dynamic_variables
		S30::Int
		S50::Int
		IF1::Int
		IF2::Int
		IF3::Int
		SC1::Int
		SC2::Int
		SC3::Int
		SC4::Int
		PIC::Int
		EFG::Int
		EFD::Int
		EGT::Int
		EGD::Int
		EFT::Int
		RF3::Int
		RGT::Int
		RGD::Int
		b1t::Int
		ta1::Int
		ta2::Int
		tc1::Int
		td1::Int
		te2::Int
		tf1::Int
		tg2::Int
		tg3::Int
		th1::Int
		ti1::Int
		tk1::Int
		t1l::Int
		t2l::Int
		t3l::Int
		t4l::Int
		t5l::Int
		tm1::Int
		tn1::Int
		tp1::Int
		tp2::Int
		tp3::Int
		tq1::Int
		tq2::Int
		tr2::Int
		tr3::Int
		tr4::Int
		tr5::Int
		ts1::Int
		ts2::Int
		ts3::Int
		ts5::Int
		tt1::Int
		tt2::Int
		tt4::Int
		tv1::Int
		tv2::Int
		tw1::Int
		ty1::Int
		tA1::Int
		tA2::Int
		tC1::Int
		tD1::Int
		tE2::Int
		tF1::Int
		tG2::Int
		tG3::Int
		tH1::Int
		tI1::Int
		tK1::Int
		tL1::Int
		tL2::Int
		tL3::Int
		tL4::Int
		tL5::Int
		tM1::Int
		tN1::Int
		tP1::Int
		tP2::Int
		tP3::Int
		tQ1::Int
		tQ2::Int
		tR2::Int
		tR3::Int
		tR4::Int
		tR5::Int
		tS1::Int
		tS2::Int
		tS3::Int
		tS5::Int
		tT1::Int
		tT2::Int
		tT4::Int
		tV1::Int
		tV2::Int
		tW1::Int
		tY1::Int
		a1t::Int
		a2t::Int
		c1t::Int
		d1t::Int
		e2t::Int
		f1t::Int
		g2t::Int
		g3t::Int
		h1t::Int
		i1t::Int
		k1t::Int
		l1t::Int
		l2t::Int
		l3t::Int
		l4t::Int
		l5t::Int
		m1t::Int
		n1t::Int
		p1t::Int
		p2t::Int
		p3t::Int
		q1t::Int
		q2t::Int
		r2t::Int
		r3t::Int
		r4t::Int
		r5t::Int
		s1t::Int
		s2t::Int
		s3t::Int
		s5t::Int
		t1t::Int
		t2t::Int
		t4t::Int
		v1t::Int
		v2t::Int
		w1t::Int
		y1t::Int
		RNA::Int
		S70::Int
		rna::Int
		PRT::Int
		PTC::Int
	end

indices = indices_dynamic_variables(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

# ---------------------------species_and_sequences! ----------------------------------------------

struct species_and_sequences
	seq_names::Vector{String}
	no_of_seqs::Int32
	µ::Float64
	avail_rib::Float64
	Tag::Vector{SubString{String}}
	Seq::Vector{String}
	Tagcell::Vector{String}
	Name::Vector{String}
	cod::Vector{String}	
	efac::Vector{Float64}
	as::Vector{String}
	tRNAs::Vector{String}	
	codAS::Vector{String}	
	t3tRNA::Vector{String}	
	AA::Vector{String}
	Tag_constant_vars::Vector{String}
	no_of_tRNAs::Int32
	no_of_variables::Int32
	cod_idx::Dict{String, Int64}
end

# ------------------------------ efficiency_factor_initialization ------------------------------------------------------

struct flags 
	positive_aminoacid::Bool
	anti_shinedalgarno::Bool
	dynamic_efficiency::Bool 
end

struct dynamic_efficiency_factors
	T_WC::Matrix{Float64}
	T_WB::Matrix{Float64}
	T_NC::Matrix{Float64}
	T_NN::Matrix{Float64}
end

# ------------------------------ calculation_of_initial_concentration ----------------------------------------------------

struct initial_concentration 
	Seq_preintegration::Vector{LongSequence{RNAAlphabet{2}}}
	constant_vars::Vector{Float64}
	c_TC::Vector{Float64}
	c_ARS::Vector{Float64}
	c_tot_EFG::Float64
end

# --------------------------------- parameters --------------------------------------------------------------------------

# for rhs_metabolites -------------------------------------

struct initiation
	k_IN_70SIC::Float64
	KM_mRNA::Float64
	KM_50S::Float64 
	k_IN_IF2D::Float64
	k_on_30S_IF3::Float64
	k_off_30SIF3::Float64
	k_on_30SIF3_IF1::Float64
	k_off_30SIF3IF1::Float64
	k_on_30SIF3_IF2::Float64
	k_off_30SIF3IF2::Float64
	k_on_30SIF3IF1_IF2::Float64
	k_off_30SIF3IF1IF2::Float64
	k_on_30SIF3IF2_IF1::Float64
	k_off_30SIF3IF2IF1::Float64
	k_on_30SIF123_fMet::Float64
	k_off_30SIF123fMet::Float64
end

struct EF_recycling
	k_EFGGDP_Diss::Float64
	k_EFG_GDP_Ass::Float64
	k_EFG_GTP_Ass::Float64
	k_EFGGTP_Diss::Float64
	k_EFTsf::Float64
	k_EFTsr::Float64
	Keq_EFTu::Float64
	KM_GTP::Float64
	KM_GDP::Float64
	KM_EFTuGDP::Float64
	KM_EFTuGTP::Float64
	Ki_EFTuGDP::Float64
	Ki_EFTuGTP::Float64
	c_EFG_t::Float64
end

struct tRNA_charging
	k_cat::Vector{Float64}
	KM_ARS_tRNA::Vector{Float64}
	KM_ARS_aa::Vector{Float64}
	KM_ARS_ATP::Vector{Float64}
end

struct ternary_complex
	k_T3j_f::Float64
	k_T3j_d::Float64
end

struct parameters_metabolites
	IN::initiation
	EF::EF_recycling
	tRNA::tRNA_charging
	t3::ternary_complex
	q_RBS::Float64
end

# for rhs_sequence ----------------------------------------

struct transcription 
	µ::Float64
	k_proc::Float64 
	k_cat_Pol::Float64
	C_RNP::Float64
	V_TCmax::Float64
	K_M_D::Float64
	k_cat_E::Float64
	C_RNE::Float64
	V_degmax::Float64
	K_M_deg::Float64
	k_D::Float64
	n_DNA_t::Float64 
	n_DNA_s::Float64
end

struct elongation 
	k_TLEmaxj::Float64
	KM_T3j::Float64
	KM_EFGGTP::Float64
end

struct termination 
	k_TE::Float64
	KM_RK::Float64
	KM_GTP::Float64
	KM_RF3GTP::Float64
	k_RF3_diss_GDP::Float64
	k_RF1_diss_GTP::Float64
	k_GTP_hydr::Float64
	k_RF1_ass::Float64
end

struct RF3_balance
	k_RF3_GDP_ass::Float64
	k_RF3_GTP_ass::Float64
	k_RF3_GDP_diss::Float64
	k_RF3_GTP_diss::Float64
end

struct PTC_recycling
	A1::Float64
	A2::Float64
	A3::Float64
	A4::Float64
	A5::Float64
	k_cat_K_M_G1::Float64
	k_cat_K_M_G2::Float64
	k_RRF::Float64
	q_RRF::Float64
	one_q_G2_k_split::Float64
end

struct parameters_sequence
	TX::transcription
	EL::elongation
	TER::termination
	RF3::RF3_balance
	PTC::PTC_recycling
end

# for dynamic_efficiencies

struct WC
	a_WC2::Float64
	a_WC3::Float64
	a_WC4::Float64
	a_WC5::Float64
	a_WC7::Float64
	a_WC8::Float64
	a_WC9::Float64
	a_WC10::Float64
	a_WC11::Float64 
end

struct WB
	a_WB2::Float64  
	a_WB3::Float64 
	a_WB4::Float64  
	a_WB5::Float64 
	a_WB7::Float64 
	a_WB8::Float64 
	a_WB9::Float64 
	a_WB10::Float64 
	a_WB11::Float64 
end

struct NC
	a_NC2::Float64
	a_NC3::Float64
	a_NC4::Float64
	a_NC5::Float64
	a_NC7::Float64
	a_NC8::Float64
	a_NC9::Float64
	a_NC10::Float64
	a_NC11::Float64
end

struct others 
	k_WC5::Float64
	a_NN2::Float64 
	a_init_selec::Float64
	k_WCrej::Float64
	k_WB5::Float64
	k_NC5::Float64
end

struct EFG_GTP_dependent
	k_WC6::Float64
	k_WCm6::Float64
	k_WC7::Float64
	k_WB6::Float64
	k_WBm6::Float64
	k_WB7::Float64
	k_NC6::Float64
	k_NCm6::Float64 
	k_NC7::Float64
end

struct parameters_dynamic_efficiency
	others::others
	WC::WC
	WB::WB
	NC::NC
	EFG_GTP::EFG_GTP_dependent
end

# ---------------------------- Preintegration ----------------------------------

struct data_struct_preintegration 
	scale::scale
    idx::indices_dynamic_variables
    species_and_sequences::species_and_sequences
    flags::flags
	dyn_eff::dynamic_efficiency_factors
    init_c::initial_concentration
    p_met::parameters_metabolites
    p_seq::parameters_sequence
	p_dyn::parameters_dynamic_efficiency
end

# ----------------------------- Integration ---------------------------------------

struct data_struct_integration 
	scale::scale
    idx::indices_dynamic_variables
    species_and_sequences::species_and_sequences
    flags::flags
	dyn_eff::dynamic_efficiency_factors
    init_c::initial_concentration
    p_met::parameters_metabolites
    p_seq::parameters_sequence
	p_dyn::parameters_dynamic_efficiency
	sequence_name::String
	Tags::Vector{SubString{String}}
	alpha_struc::Matrix{Float64}
	matrix_N::Matrix{Float64}
	tRNA_charg::Matrix{Float64}
	t3_Codon_struc::Matrix{Float64}
	efac_pos_aa::Vector{Float64}
	efac_int_sd::Vector{Float64}
	seq_names::Vector{String}
	idx_first_codon::Int32

end

end # module