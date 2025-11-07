# Species_and_sequences! ---------------------------------------------------------------------------

"""
    seqs_import_transcripts_parallel(filename::String)

Function to read in sequences from a fasta file. Throws error if the sequence is not divisible by three.
"""


function seqs_import_transcripts(filename::String)
    records = readfasta(filename)
    n = length(records)
    Seqs = Vector{Pair{String, txtl_types.SequenceData}}(undef, n)
    
    # Track problematic sequences
    problematic_sequences = String[]
    
    for i in 1:n
        try
            record = records[i]
            header  = String(record[1])
            dna     = LongDNA{2}(String(record[2]))
            
            # Check if sequence length is divisible by 3
            if length(dna) % 3 != 0
                push!(problematic_sequences, header)
                println("⚠️ Warning: Sequence '$header' has length $(length(dna)) (not divisible by 3)")
                println("   Skipping this sequence...")
                continue
            end
            
            rna     = convert(LongRNA{2}, dna)
            aa1     = BioSequences.translate(dna)
            codons  = [LongSubSeq(rna, j:j+2) for j in 1:3:length(rna)-2]
            Seqs[i] = header => txtl_types.SequenceData(dna, rna, aa1, codons)
            
        catch e
            push!(problematic_sequences, String(records[i][1]))
            println("❌ Error processing sequence '$i': $(e)")
            println("   Skipping this sequence...")
            continue
        end
    end
    
    # Filter out any nil entries from problematic sequences
    valid_seqs = filter(x -> x !== nothing, Seqs)
    
    if isempty(valid_seqs)
        error("❌ No valid sequences found! All sequences had errors.")
    end
    
    if !isempty(problematic_sequences)
        println("\n⚠️ WARNING: $(length(problematic_sequences)) sequences were skipped due to errors:")
        for seq in problematic_sequences
            println("   - $seq")
        end
        println("   Continuing with $(length(valid_seqs)) valid sequences...")
    end
    
    SeqsDict = Dict(valid_seqs)
    sequence_names = collect(keys(SeqsDict))
    
    println("✅ Successfully loaded $(length(sequence_names)) sequences")
    return SeqsDict, sequence_names
end

"""
	indices_dictionary!(tags)

17/07/25
Cedric Meyer 
Generates dictionary that has the as key the species name and as value the index of that species.
Chosen over the indices struct due to easier handling of dynamic variable addition for multiple sequences.
"""

function indices_dictionary!(tags)
	# Initialise dictionary 
	indices = Dict{String,Int}()

	# Fill in loop 
	for (idx, species) in enumerate(tags.Tag)
		indices[species] = idx
	end

	tags.idx = indices

	nothing 
end

# ------------------Efficiency factor initialization -------------------------------------------------------------------------------

function scan_for_int_sd(tags) # ONLY FOR ONE SEQUENCE
        
	# Function   : scan_for_int_sd
	# Written by : Jan Mueller (jan.mueller@ibvt.uni-stuttgart.de), Christina
	# Giakissiklis
	# Created on : 10/11/20
	# Last Update: 06/08/21
	# Purpose    : Function to analyze nucleotide sequence for internal SD like features 

	# efac_int_sd = Dict()
	Average_hybrid_energy_per_codon = zeros(Int(size(tags.hybridEnergies,1)/3))

	size(tags.hybridEnergies,1)/3
	
	for kk = 1:Int(size(tags.hybridEnergies,1)/3)
		Average_hybrid_energy_per_codon[kk,1] = maximum(tags.hybridEnergies[3*kk-2:3*kk,1]); # mean replaced with max
	end

	temp_efac_int_sd = 1 ./(0.65 .+0.28.*(Average_hybrid_energy_per_codon*(-1))); # Li 2012
	temp_efac_int_sd[temp_efac_int_sd .>=1] .= 1;
    
    # the distance of the hybridization site is 4 codons upstream the
    # A/P site of the ribosome. Therefore, the efficiency vector has to
    # be shifted as it influences the speed of that codon n+3
    
	temp_efac_int_sd = [ones(4,1); temp_efac_int_sd[1:end-4]];
	
	efac_int_sd = temp_efac_int_sd[1:end-1,1].^1;
	# efac_int_sd["sequence"] = temp_efac_int_sd[1:end-1,1].^1;
    
    # fill reference transcriptome sequence with ones
    
	# This part is not translated yet, as it is only needed for multiple sequences - will be done later

	# # Matlab
	# for jj = 2:no_of_seqs    
    #     L = size(tags.('Seq'+string(jj)).Codon,1)-1;
    #     efac_int_sd.('int_sd'+string(jj)) = ones(L,1);
    # end
	# # Julia translation
	# L = size(tags.Seq,1)-1;
	# efac_int_sd["reference"] = ones(L,1);



    return efac_int_sd
end


"""
    scan_for_pos_aa(tags, seqs, seq_names)
Analyze amino acid sequence for stretches of positive amino acids (K, R) and negative amino acids (E, D).
Written by: Marc
"""
function scan_for_pos_aa(tags, seqs, seq_names)
    # Get amino acid sequence from the current sequence (not always the first)
    # The sequence name should be the first (and only) element in seq_names
    seq_name = seq_names[1]
    aa_seq = String(seqs[seq_name].AA1[1:end-1])
    n = length(aa_seq)
    
    # Fix: efac_pos_aa should have the same length as the amino acid sequence
    efac_pos_aa = ones(n)
    
    # Initialize arrays
    positive = zeros(Int, n)
    negative = zeros(Int, n)
    sum_positive = zeros(Int, n)
    sum_negative = zeros(Int, n)
    
    # Scan for positive (K, R) and negative (E, D) amino acids
    for jj in 1:n
        if aa_seq[jj] == 'K' || aa_seq[jj] == 'R'
            positive[jj] = 1
        elseif aa_seq[jj] == 'E' || aa_seq[jj] == 'D'
            negative[jj] = -1
        end
    end
    
    # Calculate running sums (sliding window of 30 amino acids)
    for kk in 2:n
        if kk < 31
            sum_positive[kk] = sum(positive[1:kk-1])
            sum_negative[kk] = sum(negative[1:kk-1])
        else
            sum_positive[kk] = sum(positive[kk-30:kk-1])
            sum_negative[kk] = sum(negative[kk-30:kk-1])
        end
    end
    
    # Calculate net charge and filter negative values
    sum_net = sum_positive .+ sum_negative
    sum_filtered = copy(sum_net)
    for kk in 1:n
        if sum_net[kk] < 0
            sum_filtered[kk] = 0
        end
    end
    
    # Calculate efficiency factors based on positive charge density
    # Now both arrays have the same length (n)
    efac_pos_aa .= 1 ./(0.007792 .* sum_filtered.^2 .+ 0.03117 .* sum_filtered .+ 1)
    
    return efac_pos_aa
end


function efac_import(tags, arrange,wc_temp, wb_temp,nc_temp, conc2id  )

	

	no_of_tRNAs = tags.no_of_tRNAs;
		
	## Re-arrange

	wc  = zeros(62,1);
	wb  = zeros(62,2);
	nc  = zeros(62,11);

	for ii = 1:62
		kk = arrange[ii,1];
		wc[ii,:] = wc_temp[kk,:];
		wb[ii,:] = wb_temp[kk,:];
		nc[ii,:] = nc_temp[kk,:];
	end

	## Write matrix

	T_WC = zeros(62,no_of_tRNAs);    # Watson-Crick ternary complexes
	for ii = 1:62
		if wc[ii,1] != 0
			kk = conc2id[Int(wc[ii,1])];
			T_WC[ii,kk] = 1; 
		end
	end

	T_WB = zeros(62,no_of_tRNAs);    # Wobble base ternary complexes
	for ii = 1:62
		for jj=1:2
			if wb[ii,jj] != 0
				kk = conc2id[Int(wb[ii,jj])];
				T_WB[ii,kk] = 1; 
			end
		end
	end

	T_NC = zeros(62,no_of_tRNAs);    # Near-cognate ternary complexes
	for ii = 1:62
		for jj=1:11
			if nc[ii,jj] != 0
				kk = conc2id[Int(nc[ii,jj])];
				T_NC[ii,kk] = 1; 
			end
		end
	end

	T_NN = [ones(61,no_of_tRNAs); zeros(1,no_of_tRNAs)] - (T_WC + T_WB + T_NC); # Non-cognate ternary complexes

	return T_WC, T_WB, T_NC, T_NN

end

function readmatrix(filepath, sheet, range)
	xlsx = XLSX.readxlsx(filepath) # Load the workbook
	xl_sheet = xlsx[sheet] # Access the specific sheet
	data = xl_sheet[range] # Read the data from the specified range

	return data
end


function dynEffficiencyCalculation(WC,WB,NC,NN,EFG_GTP)
    ## kinetic parameters (k)
 
    # WC                # WB                # NC                # NN
    k_WC1  = 140;       k_WB1  = 140;       k_NC1  = 140;       k_NN1  = 140;   
    k_WCm1 =  85;       k_WBm1 =  85;       k_NCm1 =  85;       k_NNm1 =  85;
    k_WC2  = 190;       k_WB2  = 190;       k_NC2  = 190;      
    k_WCm2 = 0.2;       k_WBm2 =   1;       k_NCm2 =  80;      
    k_WC3  = 260;       k_WB3  =  25;       k_NC3  = 0.4;     
    k_WC4  =  10;       k_WB4  =  10;       k_NC4  =  10;      
    k_WC5  =  20;       k_WB5  = 1.6;       k_NC5  = 0.1;      
    k_WCrej =0.1;       k_WBrej= 1.1;       k_NCrej=   6;     
    k_WC6  = 150;       k_WB6  = 150;       k_NC6  = 150;      
    k_WCm6 = 140;       k_WBm6 = 140;       k_NCm6 = 140;
    k_WC7  = 100;       k_WB7  = 100;       k_NC7  = 100;
    k_WC8  =  30;       k_WB8  =  30;       k_NC8  =  30;
    k_WC9  = 500;       k_WB9  = 500;       k_NC9  = 500;
    k_WC10 =   5;       k_WB10 =   5;       k_NC10 =   5;
    k_WC11 =  20;       k_WB11 =  20;       k_NC11 =  20;
   
    ## alphas (a)
   
    # Cognate WC (WC)
    a_WC2  = (k_WCm2+k_WC3)/(k_WC2*k_WC3);
    a_WC3  = 1/k_WC3;  
    a_WC4  = 1/k_WC4; 
    a_WC5  = 1/k_WC5; 
    a_WC6f = 1/(EFG_GTP*k_WC6);
    a_WC6b = k_WCm6/(EFG_GTP*k_WC6*k_WC7);
    a_WC7  = 1/k_WC7; 
    a_WC8  = 1/k_WC8; 
    a_WC9  = 1/k_WC9; 
    a_WC10 = 1/k_WC10; 
    a_WC11 = 1/k_WC11; 
    
    # Cognate WB (WB)
    a_WB2 = (k_WB1*(k_WB3+k_WBm2))/(k_WB3*k_WBm1+k_WBm1*k_WBm2+k_WB3*k_WB2);
    a_WB3 =         (k_WB1*k_WB2)/(k_WB3*k_WBm1+k_WBm1*k_WBm2+k_WB3*k_WB2);
    a_WB4 =   (k_WB3*k_WB1*k_WB2)/(k_WB4*(k_WB3*k_WBm1+k_WBm1*k_WBm2+k_WB3*k_WB2));
    a_WB5 =   (k_WB3*k_WB1*k_WB2)/((k_WB5+k_WBrej)*(k_WB3*k_WBm1+k_WBm1*k_WBm2+k_WB3*k_WB2));
    a_WB6f = 1/(EFG_GTP*k_WB6);
    a_WB6b = k_WBm6/(EFG_GTP*k_WB6*k_WB7);
    a_WB7  = 1/k_WB7; 
    a_WB8  = 1/k_WB8; 
    a_WB9  = 1/k_WB9; 
    a_WB10 = 1/k_WB10; 
    a_WB11 = 1/k_WB11; 
    
    # Near-cognate (NC)
    a_NC2 = (k_NC1*(k_NC3+k_NCm2))/(k_NC3*k_NCm1+k_NCm1*k_NCm2+k_NC3*k_NC2);
    a_NC3 =        (k_NC1*k_NC2) /(k_NC3*k_NCm1+k_NCm1*k_NCm2+k_NC3*k_NC2);
    a_NC4 =  (k_NC3*k_NC1*k_NC2) /(k_NC4*(k_NC3*k_NCm1+k_NCm1*k_NCm2+k_NC3*k_NC2));
    a_NC5 =  (k_NC3*k_NC1*k_NC2) /((k_NC5+k_NCrej)*(k_NC3*k_NCm1+k_NCm1*k_NCm2+k_NC3*k_NC2));
    a_NC6f = 1/(EFG_GTP*k_NC6);
    a_NC6b = k_NCm6/(EFG_GTP*k_NC6*k_NC7);    
    a_NC7  = 1/k_NC7; 
    a_NC8  = 1/k_NC8; 
    a_NC9  = 1/k_NC9; 
    a_NC10 = 1/k_NC10; 
    a_NC11 = 1/k_NC11; 
    
    # Non cognate (NN)
    a_NN2 = k_NN1/k_NNm1;
    # initial selection
    a_init_selec = (1/k_WC1) * ( (k_WCm1+k_WC2)*(k_WC3+k_WCm2)/(k_WC3*k_WC2) - (k_WCm2/k_WC3) );
    # cognate WC term
    a_WC = (1+k_WCrej/k_WC5)*(a_init_selec .+ a_WC2*WC + a_WC3*WC + a_WC4*WC) + (a_WC5+a_WC6f+a_WC6b+a_WC7+a_WC8+a_WC9+a_WC10+a_WC11)*WC;
    # competition term   
    a_comp = a_init_selec * (1+k_WCrej/k_WC5) * (a_NN2*NN + (a_NC2+a_NC3+a_NC4+a_NC5)*NC + (a_WB2+a_WB3+a_WB4+a_WB5)*WB );
    # WB incorporation
    a_WB_inc = k_WB5 * a_WB5 * (a_WB6f+a_WB6b+a_WB7+a_WB8+a_WB9+a_WB10+a_WB11) * WB;
    # mis-incorporation
    a_mis_inc = k_NC5 * a_NC5 * (a_NC6f+a_NC6b+a_NC7+a_NC8+a_NC9+a_NC10+a_NC11) * NC;
   
    ## final equation
    k_eff = ( WC +  k_WB5*a_WB5*WB*a_init_selec*(1+k_WCrej/k_WC5) + k_NC5*a_NC5*NC*a_init_selec*(1+k_WCrej/k_WC5) )./ ( a_comp + a_WC + a_init_selec*(1+k_WCrej/k_WC5).*(a_WB_inc-a_mis_inc) );

    return k_eff
end   


# ---------------------------------- Release factor balance --------------------------------------------------------------


"""
	RF_balance!(y0, tags)

15/05/2025 Cedric Meyer 
Calculates the balance between RF3 bound to GTP and GDP from total RF3 amount
"""

function RF_balance!(y0, tags)

	# Function   : RF_balance V1.0
	# Written by : Marc Buschlüter (mabuschlueter@web.de), Jan Müller (jan.mueller@ibvt.uni-stuttgart.de)
	# Created on : 08/12/2020       
	# Last Update: 23/02/2022, translation to julia 26/03/2025
	# Purpose    : Solve the RF DAE System.

	## set fixed concentrations
	c_RF3_t  = y0[16]  				   # µM - RF3 
	c_GTP    = tags.constant_vars[2]   # µM - GTP
	c_GDP    = tags.constant_vars[6]   # µM - GDP 

	## Initialize RF3    

	c_RF3 = 0; # µM

	## set kinetics

	k_on_GTP_RF3=1.2*10^9*(10^-6)*60; # GTP association rate to RF3 #Marraci, Rodnina 2015 µM^-1 s^-1
	k_on_GDP_RF3=2.7*10^7*(10^-6)*60; # GDP association rate to RF3 #Marraci, Rodnina 2015 µM^-1 s^-1


	k_RF3_diss_GTP=24*60;  # GTP dissociation rate from RF3 #Maracci, Rodnina 2016 s^-1
	k_RF3_diss_GDP=0.1*60; # GDP dissociation rate from RF3 #Maracci, Rodnina 2016 s^-1

	#K_RF3_GTP=2*10^-8*10^6;   #Affinity GTP to RF3 #Maracci, Rodnina 2016
	#K_RF3_GDP=5*10^-9*10^6;   #Affinity GDP to RF3 #Maracci, Rodnina 2016 
	## solver

	function DAE!(out, du, u, p, t)
		k_on_GTP_RF3, k_on_GDP_RF3, k_RF3_diss_GTP, k_RF3_diss_GDP, c_GTP, c_GDP, c_RF3_t = p  
		c_RF3 = u[1]

		c_RF3_GTP 	= (k_on_GTP_RF3*c_GTP*c_RF3)/k_RF3_diss_GTP;
		c_RF3_GDP 	= (k_on_GDP_RF3*c_GDP*c_RF3)/k_RF3_diss_GDP;

	    out[1] = c_RF3_t-c_RF3-c_RF3_GTP-c_RF3_GDP
	end

	p = [k_on_GTP_RF3, k_on_GDP_RF3, k_RF3_diss_GTP, k_RF3_diss_GDP, c_GTP, c_GDP, c_RF3_t]
	u0 = [0.0]
	du0 = [0.0]
	differential_vars = [false]
	tspan = (0.0, 100000.0) # ? Why so much time?

	prob = DAEProblem(DAE!, du0, u0, tspan, p, differential_vars = differential_vars, 
						maxiters = 3000, abstol = 1e-6, reltol = 1e-6, verbose = false)

	sol = solve(prob, DABDF2())

	y0[16] = c_RF3 	= sol[end][]						# RF3
	y0[17] = (k_on_GTP_RF3*c_GTP*c_RF3)/k_RF3_diss_GTP; # c_RF3_GTP
	y0[18] = (k_on_GDP_RF3*c_GDP*c_RF3)/k_RF3_diss_GDP; # c_RF3_GDP 

    return nothing

end


#------------------------ preintegration -------------------------------------------------

# Start: 26.02.2025
# By: Cedric Meyer 
# Purpose: ODE function only for preintegration (excluding sequence module)

function ODE_preintegration(du, u, p, t)

    fill!(du, 0.0)  # Reset to zero, just like initializing with zeros

    rhs_metabolites!(du, u, p)
   
    nothing 
end 

"""
	preintegration(tags, y0, M)	

15/05/2025 Cedric Meyer 
calculates equilibrium before starting with transcription
Increases stability of the DAE solver
returns y1
simulates everything except tRNA related fun
"""


function preintegration(tags, y0, M)
	## preintegration - Simulation until steady state before introduction of the disturbance of interest - the sequence

	y_steady_state = copy(y0)

	TC_scale_init = tags.scale.V_TC_scale; # array [1,1]

	tags.scale.V_TC_scale = 0.0; # array [1,1]

	# Preintegration

	# Define time span for the simulation
	tspan = (0.0, 1.0)

	# Define ODEFunction to make use of mass matrix
	f = ODEFunction(ODE_preintegration, mass_matrix = M)

	# Create an ODE problem
	prob = ODEProblem(f, y_steady_state, tspan, tags);

	# Solve ODE system
	sol = solve(prob, FBDF(autodiff=false), reltol = 3e-14);

	y_steady_state[1:end] = sol[end]

	tags.scale.V_TC_scale = TC_scale_init;

	return y_steady_state

end


# ------------------ Elongation factor balance ------------------------------------------------------


function EF_balance_exponential(E_in, y0, tags)
	# Function   : EL_factor_balance V.0
	# Written by : Jan Müller
	# Rewritten by: Cedric Meyer
	# Created on : 31/07/2019
	# Last Update: 29/04/2025, renovated and translated to julia 28/03/2024
	# Purpose    : Solve the EL DAE System with ternery complexes. Avoid negative concentrations by substitution with exponential function

	# calculation of the EL-equilibrium to get the initial values for the simulation 
	
	function elongation_factor_balance!(du, u, p, t)
		# Unpack parameters 
		c_EFG_t, vector_k, K_EFG_GTP, K_EFG_GDP, Keq_EFTu, Keq_T3, 
		c_vec_T3, c_vec_aatRNA, c_vec_tRNA, c_vec_tRNA_t, E_in = p

		# Umwandlung der logarithmierten Konzentrationen in die ursprünglichen Konzentrationen
		c_EFG			= exp(u[1]);
		c_EFTu_GDP		= exp(u[2]);
		c_vec_T3   .= exp.(u[3:40]);

		# Calculate dependent variables 
		c_EFG_GTP	= (c_EFG * E_in[120]) / K_EFG_GTP;
		c_EFG_GDP	= (c_EFG * E_in[121]) / K_EFG_GDP; 

		c_EFTu_GTP	= (Keq_EFTu * c_EFTu_GDP * E_in[120])/E_in[121];

		c_vec_tRNA .= (1 .- vector_k) .* c_vec_tRNA_t;
		c_vec_aatRNA .= c_vec_T3 ./(Keq_T3 .* c_EFTu_GTP);	

		# Residuals from DAE system 
		du[1]	  = c_EFG + c_EFG_GTP + c_EFG_GDP - c_EFG_t
		du[2]	  = c_EFTu_GTP + c_EFTu_GDP + sum(c_vec_T3) - E_in[122]
		du[3:40] .= c_vec_T3 .+ c_vec_aatRNA .+ c_vec_tRNA .- c_vec_tRNA_t

		# Return sum of squared residuals 
		return nothing
	end

	# Set parameter vector 
	c_EFG_t  = tags.c_tot_EFG*tags.scale.EF_scale; #0.8 * c_30S_t;

	vector_k = [0.45;0.55;0.80;0.75;0.99;0.99;0.99;0.99;0.65;0.90
				0.90;0.90;0.99;0.99;0.99;0.99;0.99;0.70;0.99;0.99
				0.99;0.95;0.99;0.99;0.99;0.99;0.99;0.80;0.85;0.85
				0.95;0.80;0.80;0.80;0.95;0.99;0.99;0.90]; 

	K_EFG_GTP 	=   400 / 1.0e7 * 1e6;		#2.5*1e+4*1e-6	µM
	K_EFG_GDP 	=   100 / 2.7e7 * 1e6;		#2.7*1e+5*1e-6	µM
	Keq_EFTu 	=   .19;					#factor .4 in Arnold model; .19 in Romero et al
	Keq_T3      =   50;                     #1/µM

	c_vec_T3 		= zeros(38)
	c_vec_aatRNA 	= zeros(38)
	
	c_vec_tRNA 		= y0[tags.idx["ta1"]:tags.idx["ty1"]]
	c_vec_tRNA_t 	= y0[tags.idx["ta1"]:tags.idx["ty1"]];

	p = (c_EFG_t, vector_k,K_EFG_GTP, K_EFG_GDP, Keq_EFTu, Keq_T3, c_vec_T3, c_vec_aatRNA, c_vec_tRNA, c_vec_tRNA_t, E_in)

	# Define time span for the simulation
	tspan = (0.0, 100.0)
	u0 = zeros(40)
	
	tspan = (0.0, 100.0)
	
	# Definition of mass matrix
	conserved = Int.(collect(range(1, stop=40, length=40)))
	M = Matrix{Float64}(I, length(u0), length(u0))
	M[conserved, :] .= 0;

	# Define ODEFunction to make use of mass matrix
	f = ODEFunction(elongation_factor_balance!, mass_matrix = M)

	# Create an ODE problem
	prob = ODEProblem(f, u0, tspan, p);

	# Solve ODE system
	sol = solve(prob, FBDF(autodiff=false), reltol = 1e-10, abstol = 1e-10);

	x = zeros(119);
	x[1:40]		= exp.(sol[end])
	x[41]		= (x[1] * E_in[120]) / K_EFG_GTP;
	x[42]		= (x[1] * E_in[121]) / K_EFG_GDP; 
	x[43]		= (Keq_EFTu * x[2] * E_in[120]) / E_in[121];
	x[44:81] 	= x[3:40]/(Keq_T3 * x[43]);
	x[82:119]  .= (1 .- vector_k) .* c_vec_tRNA_t;

	return x
end

"""
	EF_balance!(tags, y0)

15/05/2025 Cedric Meyer 
calculation of the EL-equilibrium to get the initial values for the simulatio
"""

function EF_balance!(tags, y0, y_steady_state)


	# Elongation factor balance 

	vec_tRNA 	= tags.Tag[tags.idx["ta1"]:tags.idx["ty1"]] # all tRNA tags 
	vec_aatRNA 	= tags.Tag[tags.idx["tA1"]:tags.idx["tY1"]] # all aa-tRNA tags 
	vec_T3 		= tags.Tag[tags.idx["a1t"]:tags.idx["y1t"]] # all ternary complexes tags 

	vect = ["EFG";"EFD"; vec_T3 ;"EGT";"EGD";"EFT"; vec_aatRNA ; vec_tRNA]     

	# first aa-tRNA then tRNA 5x3 char array containing elongation factors and GTP/GDP (;'ta1';'tA1';'a1t')
	no_bal = 3*tags.no_of_tRNAs+8;  #number of metabolites for EL_balance

	EL_FAC_in 		= zeros(no_bal,1); # ELongation FACtor    
	position_elfac 	= zeros(no_bal,1);

	position_elfac = [tags.idx[species] for species in vect]

	EL_FAC_in[1:end-3] = y0[position_elfac]

	EL_FAC_in[end-2] = tags.constant_vars[2] 	# GTP
	EL_FAC_in[end-1] = tags.constant_vars[6] 	# GDP
	EL_FAC_in[end]   = tags.constant_vars[13] 	# Eut

	EL_FAC_out = EF_balance_exponential(EL_FAC_in, y0, tags);

	y_steady_state[Int.(position_elfac[1:Int(no_bal-3)])] = EL_FAC_out[1:Int(no_bal-3)]; # writes EL_FAC in y0;

end

#--------------- local tags list -------------------------------------


"""
	generate_local_tags(seqs, seq_names, tags)

02.09.2025, Cedric Meyer, Generates local tags list needed for parallel computing in Arnavs pipeline
"""
function generate_local_tags(seqs, seq_names, tags, efac_pos_aa_dict, efac_int_sd_dict)

	local_tags_list = Vector{txtl_types.data_struct}()

	for (i, seq_name) in enumerate(seq_names)
		
		local_tags 		= deepcopy(tags)
		local_tags.Seq 	= seqs[seq_name].Codon
		local_tags.sequence_name = seq_name

		# Additional tags for model species 
		species 		= ["RNA"* string(i),"S70"* string(i),"rna"* string(i),"PRT"* string(i),"PTC"* string(i)]
		local_tags.Tag 	= [strip(s) for s in vcat(tags.Tag[1:end-5],species)] 
		
		# Arrays necessary for elongation
		alpha_struc_generation!(local_tags)
		N_generation!(local_tags)
		tRNAcharging!(local_tags)
		sum_of_t3complex!(local_tags)
		
		# Efficiency factors
		local_tags.efac_pos_aa 		= efac_pos_aa_dict[local_tags.sequence_name]
		local_tags.efac_int_sd 	    = efac_int_sd_dict[local_tags.sequence_name]

		push!(local_tags_list, local_tags)
	end

	return local_tags_list

end


function alpha_struc_generation!(tags)

	# generate matrix to calculate the consumption of aa-tRNA concentration (dependent on codon and c_aa-tRNA)

	# Function   : generation of alpha_struc
	# Written by : Maike Kuschel
	# Created on : 15/12/2014
	# Last Update: 25/03/2025 Translation to julia
	# Purpose    : generate array to calculate the t3 and tRNA concentration  

	# Array to calculate the fractional part of the elongation for codon j, in
	# which the k-th T3-complex is consumed. (Arnold p.85)
	# Alpha struc is an array which combines the t3-complexes with the 
	# correct codon from the codon sequence. 
	# (in EL.rhs 
	# - calculation of c_Codon (sum of all t3-complexes) 
	# - sorting sum of all t3 complexes in the order of the codonsequence = 
	# c_Codon_Codonsequenz
	# - dividing actual t3 complex concentration by c_Codon_Codonsequenz
	# produces array N 
	#- multiplying N with alpha struc produces an array M which contains only
	# entries if the t3 complex fits the correct codon
	# - N multiplied by V_TL_EL gives the fractional part of the elongation 
	# for codon j 

	# Generate zero matrix with the tRNAs in the rows and the codons in the columns
	alpha_struc = zeros(tags.no_of_tRNAs,length(tags.Seq)-1)

	# For every tRNA get every cognate codon
	str_cod = fill("", length(tags.codAS), 3)
	for i in 1:length(tags.codAS)
		single_codons 	= split(tags.codAS[i], ",")
		stripped_codons = map(strip, single_codons)
		str_cod[i,1:length(stripped_codons)] = stripped_codons
	end

	# Compares the dependencies of the t3 complexes to the codons and saves the
	# dependencies in array "alpha_struc"
	# if the codonsequence fits with the codon (that can be bound by a specific
	# t3-complex) the codon gets the value one
	# Array with codons of gene sequence 
	strCodon = tags.Seq[1:end-1]
	for t in 1:tags.no_of_tRNAs # tRNAas
		for j = 1:size(strCodon, 1) # codons of the sequence 
			for m = 1:size(str_cod,2) # cognate codons of tRNAs
				if str_cod[t,m]  == String(strCodon[j]) # Is the current tRNA the one that reads the current codon?
					alpha_struc[t,j] = 1.0 # When yes then put a one into the matrix
				end
			end
		end
	end

	# Extension of alpha_struc with zeros to get the same dimension as N
	# (N= matrix of loaded tRNA per Codon)
	O = zeros(size((tags.Name),1)-tags.no_of_tRNAs, length(tags.Seq)-1)
	alpha_struc = vcat(O, alpha_struc)

	tags.alpha_struc = alpha_struc

end

function N_generation!(tags)

	# Function   : generation of N
	# Written by : Maike Kuschel
	# Created on : 15/04/2015
	# Last Update: 25/03/2025 Translation to julia
	# Purpose    : generate array to calculate the tRNA (release) and T3 (consumption) concentration  

	# matrix_N : rows = length of codon sequence; columns = number of all
	# possible codons. 
	# array combines the codon sequence with the Codons --> see MOD_EL_rhs
	# (multiplication of matrix_N with C_Codon(concentrations of each codon),
	# generates only entries in matrix_N, if the codons fit with the
	# sequence)

	matrix_N = zeros(length(tags.Seq)-1, size(tags.cod,1)-1)

	for i = 1:size(tags.Seq,1)-1
		for j = 1:size(tags.cod,1)-1
			if String(tags.Seq[i]) == tags.cod[j]
				matrix_N[i,j] = 1
			end
		end
	end
	
	tags.matrix_N = matrix_N
    
end

function tRNAcharging!(tags)

	# Function   : generation V_ARS for tRNA charging
	# Written by : Maike Kuschel
	# Created on : 12/01/2015
	# Last Update: -
	# Purpose    : generate array to calculate rate V_ARS for tRNA charging

	# Some amino acids are bound by more than one tRNA molecule. The array
	# tRNA_charg is generated to combine the amino acids with the correct tRNA
	# molecule

	tRNA_charg = zeros(size(tags.AA,1),size(tags.as,1));
	for t = 1 :size(tags.AA,1)
		for j = 1:size(tags.as, 1)
				if tags.AA[t] == tags.as[j]       #compares amino acid tag with the tag of tRNA, if they fit write 1                  
					tRNA_charg[t,j] = 1;                         #(see excel file, vector tags.as contains all tRNA in form of the correct aa) 
				end
		end
	end
	
	# Extension of vector tRNA_charg with zeros --> needs to have the same
	# dimension as tempV_ARS.(tempV_ARS = matrix which matches all amino acids and tRNAs; 
	# see MOD.EL_rhs-->tRNA_charging: muliplication of tempV_ARS and tRNA_charg leads 
	# to the deletion of unsuitable couples)

	P = zeros(tags.idx["ty1"], tags.idx["ty1"]) # idx of last tRNA
	O = zeros(20,tags.idx["b1t"]) 
	temp_tRNA_charg = hcat(O,tRNA_charg)
	tRNA_charg = vcat(P,temp_tRNA_charg)
	tags.tRNA_charg = tRNA_charg
end

function sum_of_t3complex!(tags)

	# Function   : generation t3_Codon_struc
	# Written by : Maike Kuschel
	# Created on : 12/01/2015
	# Last Update: -
	# Purpose    : generate array t3_Codon_struc to calculate total codon concentration

	# Some codons are bound by more than one t3 complex. The array
	# t3_Codon_struc is generated to combine the codons with the correct t3
	# complex

	sum_t3comp_temp = fill("", length(tags.t3tRNA), 3);
	for i = 1:size(tags.t3tRNA,1)                       # t3tRNA = vector (see excel file sheet t3 comp),  which combines the ternary complexes with the correct codon												
		t3tRNA_temp = split(tags.t3tRNA[i], ',');     	# divides the cell
		sum_t3comp_temp[i,1:length(t3tRNA_temp)] = t3tRNA_temp;
	end
	# Preallocation of the array t3_Codon_struc
	t3_Codon_struc = zeros(size(sum_t3comp_temp,1),size(tags.Tagcell,1));
	for i=1:size(sum_t3comp_temp,1)
		for j=1:size(sum_t3comp_temp,2)
			for k=1:size(tags.Tagcell,1)        
				if sum_t3comp_temp[i,j] == tags.Tagcell[k]		    #compares the tags in sum_t3comp_temp with 
					t3_Codon_struc[i,k] = 1.0;                      #Tagcell(=contains tags of all metabolites 
				end                                                 #of the initial vector), writes an array 
			end                                                     #if tag of ternary complex in sum_t3comp_temp
		end                                                         #fits with Tagcell--> write 1 
	end  

	tags.t3_Codon_struc = t3_Codon_struc                                                             
end