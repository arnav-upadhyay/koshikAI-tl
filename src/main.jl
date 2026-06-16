
# ============================================================================
# ADDED: Threading support 
# ============================================================================
using Base.Threads
println("🔧 Julia started with $(Threads.nthreads()) threads")

# Packages functions and data types --------
include("packages.jl") 

# Sequences 
sequence_file	= "lacS_background.fa"
hybrid_file		= "lacS.txt"

# initialization of data containers
scale = txtl_types.init_scale()         # struct that makes it possible to influence rates by scaling them

µ = 1                    # A growth rate (0.1 to 1.8) can be chosen and the proteins of the TXTL machinery will be tuned according to that growth rate and the processivity of the RNA polymerase

## flags for efficiency factors
positive_aminoacid	= true 
anti_shinedalgarno	= false      # If no hybridization file is available put flag to false. We apply the intaRNA program's intarna function (q: aSD, t: sequence, SeedBP: 5, out: tMinE) to generate a hybridization profile along the sequence. Link to the website: https://rna.informatik.uni-freiburg.de/IntaRNA/Input.jsp;jsessionid=2AE257E91E7276283F22CD77A5A8AA03
dynamic_efficiency	= true 

flags = txtl_types.flags(positive_aminoacid, anti_shinedalgarno, dynamic_efficiency, false, false, false)

## Import information about model species and sequences 
seqs, seq_names, hyb_file_dict, species_and_seqs, idx = species_and_sequences(µ, sequence_file, hybrid_file, flags)

## Calculate efficiency factors 
efac_pos_aa_dict, efac_int_sd_dict, dynamic_efficiencies = efficiency_factor_initialization(seqs, seq_names, hyb_file_dict, species_and_seqs)

## calculation of initial concentrations
y0, initial_concentrations = calculation_of_initial_concentration(species_and_seqs, seqs, seq_names, scale, idx);

## calculate and load kinetic parameters
p_met, p_seq, p_dyn = parameters(initial_concentrations, species_and_seqs);

## Release factor balance
RF_balance!(y0, initial_concentrations)

# Preintegration
## Definition of mass matrix
conserved = [idx.EFG, idx.EFD] # EFG, EFTu_GDP - Conserved species, Algebraic equations are given by mass balances

M = Matrix{Float64}(I, length(y0), length(y0)) # Creates unity matrix in the size of y0

M[conserved, :] .= 0; # setting mass matrix to zero to  remove time derivatives from those variables. This converts these differential equations into algebraic constraints

tags_pre = txtl_types.data_struct_preintegration(scale, idx, species_and_seqs, flags, dynamic_efficiencies, initial_concentrations, p_met, p_seq, p_dyn)

## Integration of the DAE system without DNA in the system 
y_steady_state = preintegration(tags_pre, y0, M);

# Elongation factor balance 
EF_balance!(tags_pre, y0, y_steady_state)

print_ef_pool_literature_audit(y_steady_state, initial_concentrations, p_met, idx)

# Build states vector for multiple sequences 
y0_all = combined_states_vector(tags_pre, y_steady_state, seqs, seq_names)

# Prepare per-sequence tags (for use in ODE)
local_tags_list = generate_local_tags(seqs, seq_names, tags_pre, efac_pos_aa_dict, efac_int_sd_dict);

# Integration
M = Matrix{Float64}(I, length(y0_all), length(y0_all))

M[conserved, :] .= 0

tspan 	= (0.0, 100.0) ## Define time span for the simulation

p_parallel = local_tags_list

f 		= ODEFunction(ODE_parallel_distributed!, mass_matrix = M) ## Define ODEFunction to make use of mass matrix

prob 	= ODEProblem(f, y0_all, tspan, p_parallel); ## Create an ODE problem

# --- Timing the integration step ---
println("\n⏳ Integrating ODE system (threaded, $(length(seq_names)) sequences)...")

sol = solve(prob, FBDF(autodiff=false), reltol = 1e-10, abstol = 1e-6) ## Solve ODE system

timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM")
println(timestamp)
println(" Integration completed!")
