# Species_and_sequences! ---------------------------------------------------------------------------

"""
12/11/2025 
species_and_sequences!(tags, sequence_file, hybrid_file)

loads species names, amino acid, their relations to codons, codon efficiencies, and compatibilites of codons to tRNAs
laods and converts DNA sequence to RNA and amino acids 
prepares matrices necessary for ODE 
eveything is stored in tags structure
"""

function species_and_sequences(µ, sequence_file, hybrid_file)

	sequence_path = joinpath(@__DIR__, "..", "data", "sequences")
	hybrid_path = joinpath(@__DIR__, "..", "data", "sequences","hybridization_files")

	seq_file = sequence_path * "\\" * sequence_file
	hyb_file = hybrid_path * "\\" * hybrid_file

	
	seqs, seq_names = seqs_import_transcripts(seq_file) # Function that transforms DNA Sequence into mRNA, amino acids and codons

	
	hyb_files = [hyb_file]
	hyb_file_dict = Dict{String, String}()
	for (i, seq_name) in enumerate(seq_names)
		hyb_file_dict[seq_name] = hyb_files[i]
	end
	seq_names = seq_names
	no_of_seqs = length(seq_names)

	if µ < 0.56
		avail_rib = -1.1191*tags.growth_rate^2 + 2.0385*tags.growth_rate + 0.1996; # calculate percentage of available ribosome, from Matamouros et al. 2021
	else
		avail_rib = 1.0 
	end

	# Starting conentrations 

	data_path = joinpath(@__DIR__, "..", "data", "input")
	metab 		= DataFrame(CSV.File(data_path * "\\metabolites.dat", delim=';', missingstring="")) # contains model species with 3 letter tags
	amino 		= DataFrame(CSV.File(data_path * "\\aminoacids.dat", delim=';', missingstring="")) # AS long name, 3 letters, short form for aatRNAs and tRNAs
	codon_conc 	= DataFrame(CSV.File(data_path * "\\Codons.dat", delim=';', missingstring="")) # Contains codons and the relative binding strength of the cognate tRNA to them
	as_cod 		= DataFrame(CSV.File(data_path * "\\AS_and_Codon.dat", delim=';', missingstring="")) # tRNAs and codons they recognize
	T3_comp 	= DataFrame(CSV.File(data_path * "\\T3_comp.dat", delim=';', missingstring="")) # codons and their compatible tRNAs (can be one or two tRNAs)
	const_vars  = DataFrame(CSV.File(data_path * "\\Constant_variables.dat", delim=';', missingstring="")) # Tags of constant variables  
			

	# Tag structure

	# Additional tags for model species 
	species 		= ["RNA","S70","rna","PRT","PTC"]
	Tag = [strip(String(s)) for s in vcat(metab.Tag, species)]
	
	Seq 		= seqs[seq_names[1]].Codon
	Tagcell 	= metab.Tag
	Name 		= metab.Name
	cod 		= map(strip,codon_conc.cod)
	efac 		= codon_conc.efac
	as 			= as_cod.as
	tRNAs 		= as_cod.tRNAas
	codAS 		= as_cod.codAS # contains codons related to tRNAs
	t3tRNA 		= T3_comp.t3tRNA
	AA 			= amino.triple
	Tag_constant_vars = const_vars.Tag

	no_of_tRNAs = 38 # Model can be adjusted to organisms with different tRNA number. E. coli has 38
	no_of_variables = length(Name)

	idx = create_indices_struct(Tag)

	# indices_dictionary!(tags)

	species_and_seqs = txtl_types.species_and_sequences(seq_names, no_of_seqs, µ, avail_rib, Tag, Seq, Tagcell, Name, cod, efac, as, tRNAs, codAS, t3tRNA, AA, Tag_constant_vars, no_of_tRNAs, no_of_variables)

	return seqs, seq_names, hyb_file_dict, species_and_seqs, idx
end


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

"""
	create_indices_struct(tags, indices)

23/06/2025
Cedric Meyer
Fills indices struct with indices of all dynamic variables.
"""

function create_indices_struct(Tag)

	idx = txtl_types.indices_dynamic_variables((1:lastindex(Tag)...))
	
	return idx
end


# ------------------Efficiency factor initialization -------------------------------------------------------------------------------


function efficiency_factor_initialization(seqs, seq_names, hyb_file_dict, species_and_seqs)

	# Scan for sequence features

	# scan for internal SD-like sequences & scan for positve aa clusters positive amino acids (ph 7): 
	# lysine(K), arginine(R) (histidine not considered)
	# if unflagged than than only ones are written into the arrays

	efac_int_sd_dict = Dict{String, Vector{Float64}}()

	if flags.anti_shinedalgarno == true
		for seq_name in seq_names
			rawHybridization = readdlm(hyb_file_dict[seq_name], ';')
			hybridEnergies = rawHybridization[2:end,3]
			hybridEnergies = replace(hybridEnergies, "NA" => missing)
			hybridEnergies = ifelse.(coalesce.(hybridEnergies, 0.0) .>= 0.0, 0.0, hybridEnergies);
			efac_int_sd_dict[seq_name]= scan_for_int_sd(hybridEnergies)
		end
	else
		for seq_name in seq_names
			L = size(seqs[seq_name].Codon, 1) - 1
			efac_int_sd_dict[seq_name] = ones(L)
		end
	end;

	# Positive amino acid clusters
	efac_pos_aa_dict = Dict{String, Vector{Float64}}()
	if flags.positive_aminoacid == true
		for seq_name in seq_names
			efac_pos_aa_dict[seq_name] = scan_for_pos_aa(seqs, [seq_name])
		end
	else
		for seq_name in seq_names
			efac_pos_aa_dict[seq_name] = ones(size(seqs[seq_name].Codon, 1) - 1)
		end
	end

	## dynamic efficiency factor initialization 

	T_WC, T_WB, T_NC, T_NN = efac_import(species_and_seqs)

	dynamic_efficiencies = txtl_types.dynamic_efficiency_factors(T_WC, T_WB, T_NC, T_NN)

	return efac_pos_aa_dict, efac_int_sd_dict, dynamic_efficiencies

end


function scan_for_int_sd(hybridEnergies) # ONLY FOR ONE SEQUENCE
        
	# Function   : scan_for_int_sd
	# Written by : Jan Mueller (jan.mueller@ibvt.uni-stuttgart.de), Christina
	# Giakissiklis
	# Created on : 10/11/20
	# Last Update: 06/08/21
	# Purpose    : Function to analyze nucleotide sequence for internal SD like features 

	# efac_int_sd = Dict()
	Average_hybrid_energy_per_codon = zeros(Int(size(hybridEnergies,1)/3))

	for kk = 1:Int(size(hybridEnergies,1)/3)
		Average_hybrid_energy_per_codon[kk,1] = maximum(hybridEnergies[3*kk-2:3*kk,1]); # mean replaced with max
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
    scan_for_pos_aa(seqs, seq_names)
Analyze amino acid sequence for stretches of positive amino acids (K, R) and negative amino acids (E, D).
Written by: Marc
"""
function scan_for_pos_aa(seqs, seq_names)
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


function efac_import(species_and_seqs)

	## Import
	data_path = joinpath(@__DIR__, "..", "data", "input")
	file = data_path * "\\Vieira_data_eco_revised.xlsx";

	arrange = readmatrix(file,"switch","A2:A63");
	arrange = convert(Matrix{Int64}, arrange);
	wc_temp = readmatrix(file,"wc1","A2:A63");
	wc_temp = convert(Matrix{Int64}, wc_temp)
	wb_temp = readmatrix(file,"wb1","A2:B63");
	wb_temp = convert(Matrix{Int64}, wb_temp)
	nc_temp = readmatrix(file,"nc1","A2:K63");
	nc_temp = convert(Matrix{Int64}, nc_temp)
	conc2id = readmatrix(file,"conc_to_id","D2:D44");
	conc2id = convert(Matrix{Int64}, conc2id)

	no_of_tRNAs = species_and_seqs.no_of_tRNAs
		
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

function readmatrix(filepath::String, sheet::String, range::String)::Any #---
	xlsx = XLSX.readxlsx(filepath) # Load the workbook
	xl_sheet = xlsx[sheet] # Access the specific sheet
	data = xl_sheet[range] # Read the data from the specified range

	return data
end


function dynEffficiencyCalculation(tags, EFG_GTP, u)
    
	@unpack a_WC2, a_WC3, a_WC4, a_WC5, a_WC7, a_WC8, a_WC9, a_WC10, a_WC11 = tags.p_dyn.WC

	@unpack a_WB2, a_WB3, a_WB4, a_WB5, a_WB7, a_WB8, a_WB9, a_WB10, a_WB11 = tags.p_dyn.WB 

	@unpack a_NC2 , a_NC3 , a_NC4 , a_NC5, a_NC7 , a_NC8 , a_NC9 , a_NC10, a_NC11 = tags.p_dyn.NC 

	@unpack k_WC5, a_NN2, a_init_selec, k_WCrej, k_WB5, k_NC5 = tags.p_dyn.others

	@unpack k_WC6, k_WCm6, k_WC7, k_WB6, k_WBm6, k_WB7, k_NC6, k_NCm6, k_NC7 = tags.p_dyn.EFG_GTP

	# GTP dependent 
	a_WC6f = 1/(EFG_GTP*k_WC6);
    a_WC6b = k_WCm6/(EFG_GTP*k_WC6*k_WC7);

    a_WB6f = 1/(EFG_GTP*k_WB6);
    a_WB6b = k_WBm6/(EFG_GTP*k_WB6*k_WB7);

    a_NC6f = 1/(EFG_GTP*k_NC6);
    a_NC6b = k_NCm6/(EFG_GTP*k_NC6*k_NC7); 
	
	c_T3 = u[tags.idx.a1t : tags.idx.y1t];

	WC = tags.dyn_eff.T_WC * c_T3;
	WB = tags.dyn_eff.T_WB * c_T3;
	NC = tags.dyn_eff.T_NC * c_T3;
	NN = tags.dyn_eff.T_NN * c_T3;

    # cognate WC term
    a_WC = (1+k_WCrej/k_WC5)*(a_init_selec .+ a_WC2.*WC .+ a_WC3.*WC .+ a_WC4.*WC) + (a_WC5+a_WC6f+a_WC6b+a_WC7+a_WC8+a_WC9+a_WC10+a_WC11).*WC;
    # competition term   
    a_comp = a_init_selec * (1+k_WCrej/k_WC5) * (a_NN2.*NN .+ (a_NC2+a_NC3+a_NC4+a_NC5).*NC .+ (a_WB2+a_WB3+a_WB4+a_WB5).*WB );
    # WB incorporation
    a_WB_inc = k_WB5 * a_WB5 * (a_WB6f+a_WB6b+a_WB7+a_WB8+a_WB9+a_WB10+a_WB11) .* WB;
    # mis-incorporation
    a_mis_inc = k_NC5 * a_NC5 * (a_NC6f+a_NC6b+a_NC7+a_NC8+a_NC9+a_NC10+a_NC11) .* NC;
   
    ## final equation
	scale = a_init_selec*(1+k_WCrej/k_WC5)
    k_eff = ( WC .+  k_WB5*a_WB5.*WB.*scale + k_NC5*a_NC5.*NC.*scale )./ ( a_comp + a_WC + scale.*(a_WB_inc-a_mis_inc) );

    return k_eff
end   

# -------------------------------- Calculation of initial concentrations --------------------

"""
	calculation_of_initial_concentration!(tags)

14/05/2025 Cedric Meyer
calculates species concentrations for the predefined growth rate 
outputs y0, mutates tags
"""

function calculation_of_initial_concentration(species_and_seqs, seqs, seq_names, scale, idx)

	# Function   : calculation_of_initial_concentration
	# Written by : Jan Mueller (jan.mueller@ibvt.uni-stuttgart.de)
	# Created on : 18/08/2022
	# Last Update: 23/04/2025 - Restructuring of the code to make it more readable
	# Purpose    : script to calculate concentration from fit equations

	data_path = joinpath(@__DIR__, "..", "data", "input")


	# Set tags.Seq for preintegration (using first sequence)
	Seq_preintegration = seqs[seq_names[1]].Codon

	## setting of initial concentration vector y0
	y0 = zeros(size(species_and_seqs.Tag,1)); # One position for every state variable 
	codon_0 = zeros(size(species_and_seqs.Seq,1)-1) # one position for each codon
	y0 = vcat(y0,codon_0); # combines both vectors
	constant_vars = zeros(37) # Nucleotides, Eut, Est, RFt, RRF, amino acids, plasmid

	GR = species_and_seqs.µ

	# aminoacyl-tRNA synthetases
	CoeffvalsARS = readmatrix(data_path*"\\FitCoefficients_eco_V9.xlsx", "aaRS", "B2:D21")
	ARSconcCalc  = CoeffvalsARS[:,3] .+(CoeffvalsARS[:,1].*GR)./(CoeffvalsARS[:,2].+GR)
	c_ARS 		 = ARSconcCalc;

	# Translation machinery
	CoeffvalsMet 		= readmatrix(data_path*"\\FitCoefficients_eco_V9.xlsx", "Translation_machinery", "B2:D14")
	METconcCalc = CoeffvalsMet[:,3]+(CoeffvalsMet[:,1].*GR)./(CoeffvalsMet[:,2].+GR);
	sum_IF      = sum(METconcCalc[3:5])
	total       = 1+0.4103+0.4219
	IF1_fix_ratio = sum_IF * 0.4103/total
	IF2_fix_ratio = sum_IF * 0.4219/total
	IF3_fix_ratio = sum_IF * 1/total
	y0[idx.IF1]    = IF1_fix_ratio;
	y0[idx.IF2]    = IF2_fix_ratio;
	y0[idx.IF3]    = IF3_fix_ratio;
	S50 = METconcCalc[2]
	S50_058      = CoeffvalsMet[2,3] +(CoeffvalsMet[2,1].*0.58)./(CoeffvalsMet[2,2]+0.58);
	y0[idx.S30:idx.S50] = METconcCalc[1:2];  # 30S, 50S, IF1, IF2 & IF3
	y0[idx.EFG]    = METconcCalc[6];    # EFG
	y0[idx.EFD]    = 3;                 # EFTu_GDP has to be greater than zero for initial solving of the DAE system
	constant_vars[13:15] = METconcCalc[7:9];  # 
	y0[idx.RF3]    = METconcCalc[12];
	constant_vars[16]    = METconcCalc[13];

	# tRNAs
	CoeffvalsReltRNA 	= readmatrix(data_path*"\\FitCoefficients_eco_V9.xlsx", "reltRNA", "B2:D40")
	tRNAconcCalc = (CoeffvalsReltRNA[:,1].*exp.(-CoeffvalsReltRNA[:,2].*GR).+CoeffvalsReltRNA[:,3]).*S50 
	y0[idx.b1t:idx.ty1] = tRNAconcCalc[1:39];

	# Nucleotides 
	NucConc058 = [9630; 4870; 8290; 2730; 555;  676; 1790;  100; 281; 23.7;  100;  360];
	NucconcCalc  = NucConc058*S50/S50_058;
	constant_vars[1:12]  = NucconcCalc[1:12]; # nucleotides 

	# Amino acids
	AAConc058  = [2550;   39; 4230; 96000; 18.2; 351; 67.6; 151.5;  405; 151.5; 145;  511;  385;  3810;  569;  68;  179;  4020; 12.1;  28.9];
	AAconcCalc   = AAConc058 *S50/S50_058;        
	constant_vars[17:36] = AAconcCalc[1:20];

	c_TC  = [METconcCalc[10]; METconcCalc[11]];

	y0[idx.EFG]    = y0[idx.EFG]   .*scale.EF_scale; #EFG
	constant_vars[13:14] = constant_vars[13:14].*scale.EF_scale; #EFTu_total & EFTs_total

	# Plasmid 
	C_D			= 0.001;				 # nM --> 4.246e-3 µM/0.004246  /2.5e-3/one Gene --> 0.0010 µM
	constant_vars[37] = C_D


	## read out total concentrations for mass balances
	c_tot_EFG = y0[idx.EFG]

	initial_concentrations = txtl_types.initial_concentration(Seq_preintegration, constant_vars, c_TC, c_ARS, c_tot_EFG)

	return y0, initial_concentrations
end

# ---------------------------------- Release factor balance --------------------------------------------------------------


"""
	RF_balance!(y0, tags)

15/05/2025 Cedric Meyer 
Calculates the balance between RF3 bound to GTP and GDP from total RF3 amount
"""

function RF_balance!(y0, initial_concentrations)

	# Function   : RF_balance V1.0
	# Written by : Marc Buschlüter (mabuschlueter@web.de), Jan Müller (jan.mueller@ibvt.uni-stuttgart.de)
	# Created on : 08/12/2020       
	# Last Update: 23/02/2022, translation to julia 26/03/2025
	# Purpose    : Solve the RF DAE System.

	## set fixed concentrations
	c_RF3_t  = y0[16]  				   # µM - RF3 
	c_GTP    = initial_concentrations.constant_vars[2]   # µM - GTP
	c_GDP    = initial_concentrations.constant_vars[6]   # µM - GDP 

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


function preintegration(tags, y0, M, scale)
	## preintegration - Simulation until steady state before introduction of the disturbance of interest - the sequence

	y_steady_state = copy(y0)

	TC_scale_init = scale.V_TC_scale; # array [1,1]

	scale.V_TC_scale = 0.0; # array [1,1]

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

	scale.V_TC_scale = TC_scale_init;

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
	c_EFG_t  = tags.init_c.c_tot_EFG*tags.scale.EF_scale; #0.8 * c_30S_t;

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
	
	c_vec_tRNA 		= y0[tags.idx.ta1:tags.idx.ty1]
	c_vec_tRNA_t 	= y0[tags.idx.ta1:tags.idx.ty1];

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

	vec_tRNA 	= tags.species_and_sequences.Tag[tags.idx.ta1:tags.idx.ty1] # all tRNA tags 
	vec_aatRNA 	= tags.species_and_sequences.Tag[tags.idx.tA1:tags.idx.tY1] # all aa-tRNA tags 
	vec_T3 		= tags.species_and_sequences.Tag[tags.idx.a1t:tags.idx.y1t] # all ternary complexes tags 

	vect = ["EFG";"EFD"; vec_T3 ;"EGT";"EGD";"EFT"; vec_aatRNA ; vec_tRNA]     

	# first aa-tRNA then tRNA 5x3 char array containing elongation factors and GTP/GDP (;'ta1';'tA1';'a1t')
	no_bal = 3*tags.species_and_sequences.no_of_tRNAs+8;  #number of metabolites for EL_balance

	EL_FAC_in 		= zeros(no_bal,1); # ELongation FACtor    
	position_elfac 	= zeros(no_bal,1);

	position_elfac = [getfield(tags.idx, Symbol(f)) for f in vect]

	EL_FAC_in[1:end-3] = y0[position_elfac]

	EL_FAC_in[end-2] = tags.init_c.constant_vars[2] 	# GTP
	EL_FAC_in[end-1] = tags.init_c.constant_vars[6] 	# GDP
	EL_FAC_in[end]   = tags.init_c.constant_vars[13] 	# Eut

	EL_FAC_out = EF_balance_exponential(EL_FAC_in, y0, tags);

	y_steady_state[Int.(position_elfac[1:Int(no_bal-3)])] = EL_FAC_out[1:Int(no_bal-3)]; # writes EL_FAC in y0;

end

#--------------- local tags list -------------------------------------


"""
	generate_local_tags(seqs, seq_names, tags)

02.09.2025, Cedric Meyer, Generates local tags list needed for parallel computing in Arnavs pipeline
"""
function generate_local_tags(seqs, seq_names, tags, efac_pos_aa_dict, efac_int_sd_dict)

	local_tags_list = Vector{txtl_types.data_struct_integration}()

	for (i, seq_name) in enumerate(seq_names)

		
		sequence 		= seqs[seq_name].Codon

		# Additional tags for model species 
		species = ["RNA"* string(i),"S70"* string(i),"rna"* string(i),"PRT"* string(i),"PTC"* string(i)]
		Tags 	= [strip(s) for s in vcat(tags.species_and_sequences.Tag[1:end-5],species)] 
		
		# Arrays necessary for elongation
		alpha_struc 	= alpha_struc_generation(tags)
		matrix_N 		= N_generation(tags)
		tRNA_charg 		= tRNAcharging(tags)
		t3_Codon_struc 	= sum_of_t3complex(tags)
		
		# Efficiency factors
		efac_pos_aa 		= efac_pos_aa_dict[seq_name]
		efac_int_sd 	    = efac_int_sd_dict[seq_name]

		species_and_seqs	= txtl_types.species_and_sequences(tags.species_and_sequences.seq_names, tags.species_and_sequences.no_of_seqs, tags.species_and_sequences.µ, 
		tags.species_and_sequences.avail_rib, tags.species_and_sequences.Tag, sequence, tags.species_and_sequences.Tagcell, tags.species_and_sequences.Name, tags.species_and_sequences.cod, 
		tags.species_and_sequences.efac, tags.species_and_sequences.as, tags.species_and_sequences.tRNAs, tags.species_and_sequences.codAS, tags.species_and_sequences.t3tRNA, tags.species_and_sequences.AA, 
		tags.species_and_sequences.Tag_constant_vars, tags.species_and_sequences.no_of_tRNAs, tags.species_and_sequences.no_of_variables)

		idx_first_codon = length(Tags)+1
		
		local_tags 		= txtl_types.data_struct_integration(tags.scale, tags.idx, species_and_seqs, tags.flags, tags.dyn_eff, 
																tags.init_c, tags.p_met, tags.p_seq, tags.p_dyn, seq_name, Tags, alpha_struc, 
																matrix_N, tRNA_charg, t3_Codon_struc, efac_pos_aa, efac_int_sd, seq_names, idx_first_codon)

		push!(local_tags_list, local_tags)
	end

	return local_tags_list

end


function alpha_struc_generation(tags)

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
	alpha_struc = zeros(tags.species_and_sequences.no_of_tRNAs,length(tags.species_and_sequences.Seq)-1)

	# For every tRNA get every cognate codon
	str_cod = fill("", length(tags.species_and_sequences.codAS), 3)
	for i in 1:length(tags.species_and_sequences.codAS)
		single_codons 	= split(tags.species_and_sequences.codAS[i], ",")
		stripped_codons = [strip(String(c)) for c in single_codons] #changed -11.11.25-arnav
		str_cod[i,1:length(stripped_codons)] = stripped_codons
	end

	# Compares the dependencies of the t3 complexes to the codons and saves the
	# dependencies in array "alpha_struc"
	# if the codonsequence fits with the codon (that can be bound by a specific
	# t3-complex) the codon gets the value one
	# Array with codons of gene sequence 
	strCodon = tags.species_and_sequences.Seq[1:end-1]
	for t in 1:tags.species_and_sequences.no_of_tRNAs # tRNAas
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
	O = zeros(size((tags.species_and_sequences.Name),1)-tags.species_and_sequences.no_of_tRNAs, length(tags.species_and_sequences.Seq)-1)
	alpha_struc = vcat(O, alpha_struc)

	return alpha_struc

end

function N_generation(tags)

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

	matrix_N = zeros(length(tags.species_and_sequences.Seq)-1, size(tags.species_and_sequences.cod,1)-1)

	for i = 1:size(tags.species_and_sequences.Seq,1)-1
		for j = 1:size(tags.species_and_sequences.cod,1)-1
			if String(tags.species_and_sequences.Seq[i]) == tags.species_and_sequences.cod[j]
				matrix_N[i,j] = 1
			end
		end
	end
	
	return matrix_N
    
end

function tRNAcharging(tags)

	# Function   : generation V_ARS for tRNA charging
	# Written by : Maike Kuschel
	# Created on : 12/01/2015
	# Last Update: -
	# Purpose    : generate array to calculate rate V_ARS for tRNA charging

	# Some amino acids are bound by more than one tRNA molecule. The array
	# tRNA_charg is generated to combine the amino acids with the correct tRNA
	# molecule

	tRNA_charg = zeros(size(tags.species_and_sequences.AA,1),size(tags.species_and_sequences.as,1));
	for t = 1 :size(tags.species_and_sequences.AA,1)
		for j = 1:size(tags.species_and_sequences.as, 1)
				if tags.species_and_sequences.AA[t] == tags.species_and_sequences.as[j]       #compares amino acid tag with the tag of tRNA, if they fit write 1                  
					tRNA_charg[t,j] = 1;          #(see excel file, vector tags.as contains all tRNA in form of the correct aa) 
				end
		end
	end
	
	# Extension of vector tRNA_charg with zeros --> needs to have the same
	# dimension as tempV_ARS.(tempV_ARS = matrix which matches all amino acids and tRNAs; 
	# see MOD.EL_rhs-->tRNA_charging: muliplication of tempV_ARS and tRNA_charg leads 
	# to the deletion of unsuitable couples)

	P = zeros(tags.idx.ty1, tags.idx.ty1) # idx of last tRNA
	O = zeros(20,tags.idx.b1t) 
	temp_tRNA_charg = hcat(O,tRNA_charg)
	tRNA_charg = vcat(P,temp_tRNA_charg)
	return tRNA_charg
end

function sum_of_t3complex(tags)

	# Function   : generation t3_Codon_struc
	# Written by : Maike Kuschel
	# Created on : 12/01/2015
	# Last Update: -
	# Purpose    : generate array t3_Codon_struc to calculate total codon concentration

	# Some codons are bound by more than one t3 complex. The array
	# t3_Codon_struc is generated to combine the codons with the correct t3
	# complex

	sum_t3comp_temp = fill("", length(tags.species_and_sequences.t3tRNA), 3);
	for i = 1:size(tags.species_and_sequences.t3tRNA,1)                       # t3tRNA = vector (see excel file sheet t3 comp),  which combines the ternary complexes with the correct codon												
		t3tRNA_temp = split(tags.species_and_sequences.t3tRNA[i], ',');     	# divides the cell
		sum_t3comp_temp[i,1:length(t3tRNA_temp)] = t3tRNA_temp;
	end
	# Preallocation of the array t3_Codon_struc
	t3_Codon_struc = zeros(size(sum_t3comp_temp,1),size(tags.species_and_sequences.Tagcell,1));
	for i=1:size(sum_t3comp_temp,1)
		for j=1:size(sum_t3comp_temp,2)
			for k=1:size(tags.species_and_sequences.Tagcell,1)        
				if sum_t3comp_temp[i,j] == tags.species_and_sequences.Tagcell[k]		    #compares the tags in sum_t3comp_temp with 
					t3_Codon_struc[i,k] = 1.0;                      #Tagcell(=contains tags of all metabolites 
				end                                                 #of the initial vector), writes an array 
			end                                                     #if tag of ternary complex in sum_t3comp_temp
		end                                                         #fits with Tagcell--> write 1 
	end  

	return t3_Codon_struc                                                             
end

"""
	combined_states_vector(tags, y_steady_state, seqs, seq_names)

12.11.2025
Function to build a state vector that includes metabolites, outputvariables for sequ. 1, codons sequ. 1, output 2, sequence 2 ...
"""

function combined_states_vector(tags, y_steady_state, seqs, seq_names)

	# 1. Build combined state vector for this set
	y0_metab_output = y_steady_state[1:length(tags.species_and_sequences.Tag)]                                        # All species without the codons
	y0_metab = y0_metab_output[1:end-5]                                                         # Species of the metabolite pool
	y0_output = y0_metab_output[end - 4: end]                                                   # Output species  "RNA" "S70"  "rna"  "PRT"  "PTC"
	ribosome_positions = [zeros(length(seqs[seq_name].Codon) - 1) for seq_name in seq_names]    # All codons per sequence 
	output_sequences = reduce(vcat, ([y0_output; ribosome_positions[i]...] for i in eachindex(ribosome_positions))) # One set of output variables per sequence 
	y0_all = vcat(y0_metab, output_sequences);                                                   # concatenate everything: metabolites pool, output 1, sequence 1, output 2, sequence 2, output 3, sequence 3 ...

	return y0_all
end

# ------------------------------- Integration -----------------------------------------

# 3. Define the ODE function for solving (distributed, focused on derivatives)
function ODE_parallel_distributed!(du, u, p, t)

    local_tags_list = p				# Extraction of local_tags_list from p
    tags            = local_tags_list[1]

    n_metab = length(tags.species_and_sequences.Tag) -5 		# number of metabolites in the shared pool that is used by all sequences 

    @views y_metab =  u[1:n_metab]				# Save all pool species from u in a vector 
    idx = n_metab						# Index of last pool species
   
    y_codons_lengths = Int[length(tags.species_and_sequences.Seq) + 4 for seq_name in tags.seq_names]	# List with sequence lengths - 1 stop codon + 5 observable variables

    y_codons = Vector{SubArray{Float64,1,typeof(u),Tuple{UnitRange{Int64}},true}}(undef, length(y_codons_lengths))
    idx = n_metab
    @views for (i,L) in enumerate(y_codons_lengths)
        y_codons[i] = u[idx+1:idx+L]  # ✅ View statt Kopie
        idx += L
    end

    
    n_seq = tags.species_and_sequences.no_of_seqs 		# save number of sequences
    du_metab_all = [zeros(n_metab) for _ in 1:n_seq]	# List to store a du vector for the metabolite pool for every sequence 
    du_codons = Vector{Vector{Float64}}(undef, n_seq)	# List to store a du vector for the codons for every sequence 

    # Distributed parallelization of per-sequence rate calculations
    results = map(1:n_seq) do i						    # Mapping der Schleife auf die worker
        local_tags = local_tags_list[i]		# Extraction of local tags
        u_seq = vcat(y_metab, y_codons[i])				# Sequenzspezifischer state vector
        du_seq = zeros(length(u_seq))					# Sequenzspezifischer rates vector
        ODE!(du_seq, u_seq, local_tags, t)	            # Aufruf des ODE systems
        return du_seq									# Rückgabe von du, abgespeichert in results
    end

    @views for i in 1:n_seq								# Loop over number of sequences
        du_metab_all[i] .= results[i][1:n_metab]	    # Abspeichern aller rates der Metabolite 
        du_codons[i] = results[i][n_metab+1:end]		# Abspeichern aller rates für die Ribosomen
    end

    du_metab = zeros(n_metab)						    # Summation der n_metab du
    for i in 1:n_seq
        du_metab .+= du_metab_all[i]
    end

    @views du[1:n_metab] .= du_metab                    # Zusammensetzung des gesamten rates vectors
    idx = n_metab
    for i in 1:n_seq
        len = length(du_codons[i])
        @views du[idx+1:idx+len] .= du_codons[i]
        idx += len
    end			

    return nothing
end