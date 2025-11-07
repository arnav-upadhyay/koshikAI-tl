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

mutable struct scale_struct 
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
    return scale_struct(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 .* ones(20), 0.0, 0.0)
end

mutable struct data_struct 
    organism::String
    no_of_tRNAs::Int32
    avail_rib::Float64
    scale::scale_struct
    no_of_seqs::Int32
    growth_rate::Float64

    Tag::Vector{String}
    Seq::Vector{LongSequence{RNAAlphabet{2}}}
    seq_names::Vector{String}
    Tagcell::Vector{String}
    Name::Vector{String}

    cod::Vector{String}
    efac::Vector{Float64}
    as::Vector{String}
    tRNAs::Vector{String}
    codAS::Vector{String}

    t3tRNA::Vector{String}
    AA::Vector{String}
    alpha_struc::Matrix{Float64}
    matrix_N::Matrix{Float64}
    tRNA_charg::Matrix{Float64}

    t3_Codon_struc::Matrix{Float64}
    c_ARS::Vector{Float64}
    c_TC::Vector{Float64}
    c_tot_EFG::Float64
    efac_pos_aa::Vector{Float64}

    efac_int_sd::Vector{Float64}
    hybridEnergies::Vector{Float64}

    T_WC::Matrix{Float64}
    T_WB::Matrix{Float64}
    T_NC::Matrix{Float64}
    T_NN::Matrix{Float64}
    dynamic_efficiency::Bool

    p_met::Vector{Any}
    p_seq::Vector{Any}
    constant_vars::Vector{Float64}
	idx::Dict{String,Int}
	no_of_variables::Int32

	Tag_constant_vars::Vector{String}
    sequence_name::String

end

function init_data()
    return data_struct(
        "", 0, 0.0, init_scale(), 0, 0.0,
        String[], Vector{LongSequence{RNAAlphabet{2}}}(), String[], String[], String[],
        String[], Float64[], String[], String[], String[],
        String[], String[], zeros(0,0), zeros(0,0), zeros(0,0),
        zeros(0,0), Float64[], Float64[], 0.0, Float64[],
        Float64[], Float64[],
        zeros(0,0), zeros(0,0), zeros(0,0), zeros(0,0), true,
        Any[], Any[], Vector{Float64}(), Dict{String,Int}(), 0,
        Vector{String}(), ""
    )
end


end # module