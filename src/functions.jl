# Species_and_sequences! ---------------------------------------------------------------------------

"""
12/11/2025 
species_and_sequences!(tags, sequence_file, hybrid_file)

loads species names, amino acid, their relations to codons, codon efficiencies, and compatibilites of codons to tRNAs
laods and converts DNA sequence to RNA and amino acids 
prepares matrices necessary for ODE 
eveything is stored in tags structure
"""

function species_and_sequences(genes, µ, sequence_file, hybrid_file, flags, seqs_to_use)

	
	hybrid_path = joinpath(@__DIR__, "..", "data", "sequences","hybridization_files")
	hyb_file = hybrid_path * "\\" * hybrid_file

	
	seqs_all, seq_names_all = seqs_import_transcripts(sequence_file, seqs_to_use) # Function that transforms DNA Sequence into mRNA, amino acids and codons


	# !!!!!! Shine dalgarno needs to be adapted to competition mode !!!!!!!!!!!!!!!!
	hyb_file_dict = Dict{String, String}()
	# if flags.anti_shinedalgarno == true
	# 	hyb_files = [hyb_file]
	# 	for (i, seq_name) in enumerate(seq_names)
	# 		hyb_file_dict[seq_name] = hyb_files[i]
	# 	end
	# else
	# 		hyb_file_dict["error"] = "No shine dalgarno sequence interaction"
	# end

	# ------------ for competition mode
	seq_names = genes		
	# filter seqs_all
	seqs = Dict(k => seqs_all[k] for k in seq_names if haskey(seqs_all, k));
	no_of_seqs = length(seq_names)

	if µ < 0.56
		avail_rib = -1.1191*tags.growth_rate^2 + 2.0385*tags.growth_rate + 0.1996; # calculate percentage of available ribosome, from Matamouros et al. 2021
	else
		avail_rib = 1.0 
	end

	# Starting conentrations 

	data_path = joinpath(@__DIR__, "..", "data", "input")
	metab_path = joinpath(data_path, "Metabolites.dat")
	amino_path = joinpath(data_path, "aminoacids.dat")
	codons_path = joinpath(data_path, "Codons.dat")
	AC_Cod_path = joinpath(data_path, "AS_and_Codon.dat")
	T3_path = joinpath(data_path, "T3_comp.dat")
	Const_path = joinpath(data_path, "Constant_variables.dat")

	metab 		= DataFrame(CSV.File(metab_path, delim=';', missingstring="")) # contains model species with 3 letter tags
	amino 		= DataFrame(CSV.File(amino_path, delim=';', missingstring="")) # AS long name, 3 letters, short form for aatRNAs and tRNAs
	codon_conc 	= DataFrame(CSV.File(codons_path, delim=';', missingstring="")) # Contains codons and the relative binding strength of the cognate tRNA to them
	as_cod 		= DataFrame(CSV.File(AC_Cod_path, delim=';', missingstring="")) # tRNAs and codons they recognize
	T3_comp 	= DataFrame(CSV.File(T3_path, delim=';', missingstring="")) # codons and their compatible tRNAs (can be one or two tRNAs)
	const_vars  = DataFrame(CSV.File(Const_path, delim=';', missingstring="")) # Tags of constant variables  
			

	# Tag structure

	# Additional tags for model species 
	species 		= ["S70","rna","PRT","PTC"]
	Tag = [strip(String(s)) for s in vcat(metab.Tag, species)]
	
	Seq 		= String.(seqs[seq_names[1]].Codon)
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

	# Codon lookup table for V_TL_EL 

	cod_index = Dict{String,Int}()
	for i in eachindex(cod)
		cod_index[cod[i]] = i
	end

	species_and_seqs = txtl_types.species_and_sequences(seq_names, no_of_seqs, µ, avail_rib, Tag, Seq, Tagcell, Name, cod, efac, as, 
														tRNAs, codAS, t3tRNA, AA, Tag_constant_vars, no_of_tRNAs, no_of_variables, cod_index)

	return seqs, seq_names, hyb_file_dict, species_and_seqs, idx
end


"""
    seqs_import_transcripts(filename::String)

Function to read in sequences from a fasta file. Throws error if the sequence is not divisible by three.
"""


function seqs_import_transcripts(sequence_file::String, seqs_to_use)
	
	sequence_path = joinpath(@__DIR__, "..", "data", "sequences")
	# filename = sequence_path * "\\" * sequence_file
	filename = joinpath(sequence_path, sequence_file)

    records = readfasta(filename)
    n = length(records)
	no_of_seqs = seqs_to_use[2] - seqs_to_use[1] + 1
    Seqs = Vector{Pair{String, txtl_types.SequenceData}}(undef, no_of_seqs)

    
    # Track problematic sequences
    problematic_sequences = String[]
    
    for i in 1:n
		if i >= seqs_to_use[1] && i <= seqs_to_use[2]
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
    
    # println("✅ Successfully loaded $(length(sequence_names)) sequences")
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


function efficiency_factor_initialization(seqs, seq_names, hyb_file_dict, species_and_seqs, flags)

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
	file = joinpath(data_path, "Vieira_data_eco_revised.xlsx")
	# file = data_path * "\\Vieira_data_eco_revised.xlsx";

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
	
	@views c_T3 = u[tags.idx.a1t : tags.idx.y1t];

	WC = tags.dyn_eff.T_WC * c_T3;
	WB = tags.dyn_eff.T_WB * c_T3;
	NC = tags.dyn_eff.T_NC * c_T3;
	NN = tags.dyn_eff.T_NN * c_T3;

    # Preallocate arrays 
    a_WC = zeros(62)
    a_comp = zeros(62)
    a_WB_inc = zeros(62)
    k_eff = zeros(62)

    # cognate WC term
    @inbounds for i in eachindex(a_WC)
    a_WC[i] =
        (1 + k_WCrej/k_WC5) * (
            a_init_selec +
            (a_WC2 + a_WC3 + a_WC4) * WC[i]
        ) +
        (a_WC5 + a_WC6f + a_WC6b + a_WC7 + a_WC8 + a_WC9 + a_WC10 + a_WC11) * WC[i]
    end 

    # # competition term   
    scale = a_init_selec * (1 + k_WCrej/k_WC5)
    sum_NC = a_NC2 + a_NC3 + a_NC4 + a_NC5
    sum_WB = a_WB2 + a_WB3 + a_WB4 + a_WB5

    @inbounds for i in eachindex(a_comp)
        a_comp[i] = scale * (
            a_NN2 * NN[i] +
            sum_NC * NC[i] +
            sum_WB * WB[i]
        )
    end

    # # WB incorporation
    coeff = k_WB5 * a_WB5 * (a_WB6f + a_WB6b + a_WB7 + a_WB8 + a_WB9 + a_WB10 + a_WB11)

    @inbounds for i in eachindex(a_WB_inc)
        a_WB_inc[i] = coeff * WB[i]
    end

    # mis-incorporation
    a_mis_inc = k_NC5 * a_NC5 * (a_NC6f+a_NC6b+a_NC7+a_NC8+a_NC9+a_NC10+a_NC11) .* NC;
   
    # ## final equation
    @inbounds for i in eachindex(k_eff)
        num = WC[i] +
            k_WB5 * a_WB5 * WB[i] * scale +
            k_NC5 * a_NC5 * NC[i] * scale

        den = a_comp[i] +
            a_WC[i] +
            scale * (a_WB_inc[i] - a_mis_inc[i])

        k_eff[i] = num / den
    end

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
	CoeffvalsARS = readmatrix(data_path*"/FitCoefficients_eco_V9.xlsx", "aaRS", "B2:D21")
	ARSconcCalc  = CoeffvalsARS[:,3] .+(CoeffvalsARS[:,1].*GR)./(CoeffvalsARS[:,2].+GR)
	c_ARS 		 = ARSconcCalc;

	# Translation machinery
	CoeffvalsMet 		= readmatrix(data_path*"/FitCoefficients_eco_V9.xlsx", "Translation_machinery", "B2:D14")
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
	CoeffvalsReltRNA 	= readmatrix(data_path*"/FitCoefficients_eco_V9.xlsx", "reltRNA", "B2:D40")
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
<

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

# long_met_rhs_01
# like: MOD_met_rhs.m 
# Start: 13.03.2025
# By: Cedric Meyer 
# Purpose: ODE system of the metabolite module, long form without the loops 

# function ODE_met_long(du, u, p::data_struct, t)

# Can mod_met_rhs and mod_seq_rhs be combined so that du can be directly saved into du?

function rhs_metabolites_preintegration!(du, u, p)


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
	S30, S50, IF1, IF2, IF3, SC1, SC2, SC3, SC4, PIC, EFG, EFD, EGT, EGD, EFT, RF3, RGT, RGD, b1t, ta1, ta2, tc1, td1, te2, tf1, tg2, tg3, th1, ti1, tk1, t1l, t2l, t3l, t4l, t5l, tm1, tn1, tp1, tp2, tp3, tq1, tq2, tr2, tr3, tr4, tr5, ts1, ts2, ts3, ts5, tt1, tt2, tt4, tv1, tv2, tw1, ty1, tA1, tA2, tC1, tD1, tE2, tF1, tG2, tG3, tH1, tI1, tK1, tL1, tL2, tL3, tL4, tL5, tM1, tN1, tP1, tP2, tP3, tQ1, tQ2, tR2, tR3, tR4, tR5, tS1, tS2, tS3, tS5, tT1, tT2, tT4, tV1, tV2, tW1, tY1, a1t, a2t, c1t, d1t, e2t, f1t, g2t, g3t, h1t, i1t, k1t, l1t, l2t, l3t, l4t, l5t, m1t, n1t, p1t, p2t, p3t, q1t, q2t, r2t, r3t, r4t, r5t, s1t, s2t, s3t, s5t, t1t, t2t, t4t, v1t, v2t, w1t, y1t, RNA, S70, rna, PRT, PTC, RIB1 = u 

	RNA = tags.RNA 

	# Define RIB1
	idx_first_codon = length(tags.species_and_sequences.Tag)+1
	@views c_RIB	 = u[idx_first_codon:end];

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


	du[tags.idx.S30] = - k_on_30S_IF3*S30*IF3 + k_off_30SIF3*SC1 # 13
	du[tags.idx.S50] = - V_IN_70SIC # 14
	du[tags.idx.IF1] = - k_on_30SIF3_IF1*SC1*IF1    + k_off_30SIF3IF1*SC2 - k_on_30SIF3IF2_IF1*SC3*IF1 + k_off_30SIF3IF2IF1*SC4 + V_IN_IF2D # 15
	du[tags.idx.IF2] = - k_on_30SIF3_IF2*SC1*IF2    + k_off_30SIF3IF2*SC3 - k_on_30SIF3IF1_IF2*SC2*IF2 + k_off_30SIF3IF1IF2*SC4 + V_IN_IF2D# 16
	du[tags.idx.IF3] = - k_on_30S_IF3*S30*IF3 + k_off_30SIF3*SC1 + V_IN_IF2D # 17
	du[tags.idx.SC1] = + k_on_30S_IF3*S30*IF3 - k_off_30SIF3*SC1 + k_off_30SIF3IF1*SC2  - k_on_30SIF3_IF1*SC1*IF1 + k_off_30SIF3IF2*SC3    - k_on_30SIF3_IF2*SC1*IF2 #18
	du[tags.idx.SC2] = + k_on_30SIF3_IF1*SC1*IF1 - k_off_30SIF3IF1*SC2 + k_off_30SIF3IF1IF2*SC4    - k_on_30SIF3IF1_IF2*SC2*IF2 # 19
	du[tags.idx.SC3] = + k_on_30SIF3_IF2*SC1*IF2 - k_off_30SIF3IF2*SC3 + k_off_30SIF3IF2IF1*SC4    - k_on_30SIF3IF2_IF1*SC3*IF1 # 20
	du[tags.idx.SC4] = + k_off_30SIF123fMet*PIC       - k_on_30SIF123_fMet*SC4*b1t + k_on_30SIF3IF1_IF2*SC2*IF2 - k_off_30SIF3IF1IF2*SC4 + k_on_30SIF3IF2_IF1*SC3*IF1 - k_off_30SIF3IF2IF1*SC4 # 21
	du[tags.idx.PIC] =  + k_on_30SIF123_fMet*SC4*b1t - k_off_30SIF123fMet*PIC - (((k_IN_70SIC * PIC)/ ( 1+(KM_mRNA/RNA+(KM_50S/S50)) )) .* q_RBS) # 22 30S PIC
	du[tags.idx.b1t] = + k_off_30SIF123fMet*PIC - k_on_30SIF123_fMet*SC4*b1t + V_IN_IF2D # 35, db1t
	# du[tags.idx.RNA] = - V_IN_70SIC # 170, dRNA
	du[tags.idx.S70] = V_IN_70SIC - V_IN_IF2D # 171, dS70
	du[tags.idx.rna] = V_IN_70SIC # 172, drna

	du[idx_first_codon] = V_IN_IF2D # 175, dRIB

    du[tags.idx.ta1:tags.idx.ty1] .= 0
    du[tags.idx.tA1:tags.idx.tY1] .= 0
    du[tags.idx.a1t:tags.idx.y1t] .= 0


	# Elongation factors 
	du[tags.idx.EGT] =  k_EFG_GTP_Ass * EFG * GTP - k_EFGGTP_Diss * EGT		# Rate of EFG_GTP # 25
	du[tags.idx.EGD] = -(k_EFGGDP_Diss * EGD - k_EFG_GDP_Ass * EFG * GDP)		# Rate of EFG_GDP # 26
	du[tags.idx.EFT] = (k_EFTsf * Est) * (EFD * GTP - (GDP * EFT/Keq_EFTu))/ (KM_GTP * EFD + KM_EFTuGDP * GTP + (k_EFTsf * Est * KM_EFTuGTP * GDP)/(k_EFTsr * Est * Keq_EFTu) + (k_EFTsf * Est * KM_GDP * EFT)/(k_EFTsr * Est * Keq_EFTu)+
			EFD * GTP + (k_EFTsf * Est * KM_EFTuGTP * EFD * GDP)/(k_EFTsr * Est * Keq_EFTu * Ki_EFTuGDP) + (k_EFTsf * Est * GDP * EFT) /(k_EFTsr * Est * Keq_EFTu) + 
			(KM_EFTuGDP * GTP * EFT)/ Ki_EFTuGTP) -
			+((k_T3j_f * EFT * tA1 - k_T3j_d * a1t),(k_T3j_f * EFT * tA2 - k_T3j_d * a2t),(k_T3j_f * EFT * tC1 - k_T3j_d * c1t),(k_T3j_f * EFT * tD1 - k_T3j_d * d1t),(k_T3j_f * EFT * tE2 - k_T3j_d * e2t),(k_T3j_f * EFT * tF1 - k_T3j_d * f1t),(k_T3j_f * EFT * tG2 - k_T3j_d * g2t),(k_T3j_f * EFT * tG3 - k_T3j_d * g3t),(k_T3j_f * EFT * tH1 - k_T3j_d * h1t),(k_T3j_f * EFT * tI1 - k_T3j_d * i1t),(k_T3j_f * EFT * tK1 - k_T3j_d * k1t),(k_T3j_f * EFT * tL1 - k_T3j_d * l1t),(k_T3j_f * EFT * tL2 - k_T3j_d * l2t),(k_T3j_f * EFT * tL3 - k_T3j_d * l3t),(k_T3j_f * EFT * tL4 - k_T3j_d * l4t),(k_T3j_f * EFT * tL5 - k_T3j_d * l5t),(k_T3j_f * EFT * tM1 - k_T3j_d * m1t),(k_T3j_f * EFT * tN1 - k_T3j_d * n1t),(k_T3j_f * EFT * tP1 - k_T3j_d * p1t),(k_T3j_f * EFT * tP2 - k_T3j_d * p2t),(k_T3j_f * EFT * tP3 - k_T3j_d * p3t),(k_T3j_f * EFT * tQ1 - k_T3j_d * q1t),(k_T3j_f * EFT * tQ2 - k_T3j_d * q2t),(k_T3j_f * EFT * tR2 - k_T3j_d * r2t),(k_T3j_f * EFT * tR3 - k_T3j_d * r3t),(k_T3j_f * EFT * tR4 - k_T3j_d * r4t),(k_T3j_f * EFT * tR5 - k_T3j_d * r5t),(k_T3j_f * EFT * tS1 - k_T3j_d * s1t),(k_T3j_f * EFT * tS2 - k_T3j_d * s2t),(k_T3j_f * EFT * tS3 - k_T3j_d * s3t),(k_T3j_f * EFT * tS5 - k_T3j_d * s5t),(k_T3j_f * EFT * tT1 - k_T3j_d * t1t),(k_T3j_f * EFT * tT2 - k_T3j_d * t2t),(k_T3j_f * EFT * tT4 - k_T3j_d * t4t),(k_T3j_f * EFT * tV1 - k_T3j_d * v1t),(k_T3j_f * EFT * tV2 - k_T3j_d * v2t),(k_T3j_f * EFT * tW1 - k_T3j_d * w1t),(k_T3j_f * EFT * tY1 - k_T3j_d * y1t)); # 27

	# steady-state of elongation factors

	T3  = zeros(Float64, tags.species_and_sequences.no_of_tRNAs);

	T3[1:end] = u[tags.idx.a1t:tags.idx.y1t]

	T3sum = sum(T3);

	# DAE system
	du[tags.idx.EFG] = EFG + EGT + EGD - c_EFG_t;   # EFG
	du[tags.idx.EFD] = EFT + EFD - Eut + T3sum;     # EFTu_GDP
	
	nothing
end





# Start: 26.02.2025
# By: Cedric Meyer 
# Purpose: ODE function only for preintegration (excluding sequence module)

function ODE_preintegration(du, u, p, t)

    fill!(du, 0.0)  # Reset to zero, just like initializing with zeros

    rhs_metabolites_preintegration!(du, u, p)
   
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

	# CRITICAL FIX: Clamp physically non-negative species to prevent negative values
	# Negative S70 (free 70S ribosome pool) was causing multi-sequence simulations to fail
	# When S70 goes negative in preintegration, it propagates to all sequences in y0_all
	# This ensures the steady-state vector is physically valid before being used
	y_steady_state[tags.idx.S70] = max(y_steady_state[tags.idx.S70], 0.0)
	
	# Also clamp rna (transcribed RNA) to ensure non-negative
	# RNA/rna are output species that should also be non-negative
	y_steady_state[tags.idx.rna] = max(y_steady_state[tags.idx.rna], 0.0)


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
	generate_p(local_tags_list)

09.12.2025, generates a struct that has local tags list and initialized vectors for ODE_parallel_distributed
"""

function generate_p(local_tags_list)

    # Initialiization of vectors and lists
    y_codons_lengths    = Vector{Int}(undef, length(local_tags_list))   # Initialization
    @inbounds for i in eachindex(local_tags_list)                                     # List with sequence lengths - 1 stop codon + 4 observable variables
        y_codons_lengths[i] = length(local_tags_list[i].species_and_sequences.Seq) + 3
    end	

	offsets = cumsum(vcat(0, y_codons_lengths[1:end-1]))

    n_metab = length(local_tags_list[1].species_and_sequences.Tag) -4     # number of metabolites in the shared pool that is used by all sequences 
    n_seq   = local_tags_list[1].species_and_sequences.no_of_seqs 		                # save number of sequences

    du_metab_all    = [zeros(n_metab) for _ in 1:n_seq]             # Initialization - List to store a du vector for the metabolite pool for every sequence 
    du_codons       = Vector{Vector{Float64}}(undef, n_seq)	        # Initialization - List to store a du vector for the codons for every sequence     
    y_codons        = Vector{Vector{Float64}}(undef, n_seq)			# Initialization
    results         = Vector{Vector{Float64}}(undef, n_seq)         # Initialization - List to store all sequence specific metabolite du and codon du

    du_seq          = Vector{Vector{Float64}}(undef, n_seq)	        # Initialization - List to store a du vector for the codons for every sequence     
    for (i,y_length) in enumerate(y_codons_lengths)
        L = y_length+n_metab
        du_seq[i] = zeros(L)
    end   

	du_metab = zeros(n_metab)

    vectors = txtl_types.vectors(y_codons_lengths, offsets, n_metab, n_seq, du_metab_all, du_codons, y_codons, results, du_seq, du_metab)

    p = txtl_types.parameters(local_tags_list, vectors)

    return p 
end

"""
	generate_local_tags(seqs, seq_names, tags)

02.09.2025, Cedric Meyer, Generates local tags list needed for parallel computing in Arnavs pipeline
"""
function generate_local_tags(seqs, seq_names, tags, efac_pos_aa_dict, efac_int_sd_dict, mRNA)

	local_tags_list = Vector{txtl_types.data_struct_integration}()

	for (i, seq_name) in enumerate(seq_names)

		
		sequence 		= String.(seqs[seq_name].Codon)

		# Additional tags for model species 
		species = ["S70"* string(i),"rna"* string(i),"PRT"* string(i),"PTC"* string(i)]
		Tags 	= [strip(s) for s in vcat(tags.species_and_sequences.Tag[1:end-4],species)] 
		
		# Arrays necessary for elongation
		alpha_one, alpha_two	= alpha_struc_generation(tags, sequence)
		matrix_N 				= N_generation(tags, sequence)
		tRNA_charg 				= tRNAcharging(tags)
		t3_Codon_struc 			= sum_of_t3complex(tags)

		# Efficiency factors
		efac_pos_aa 		= efac_pos_aa_dict[seq_name]
		efac_int_sd 	    = efac_int_sd_dict[seq_name]

		species_and_seqs	= txtl_types.species_and_sequences(tags.species_and_sequences.seq_names, tags.species_and_sequences.no_of_seqs, tags.species_and_sequences.µ, 
		tags.species_and_sequences.avail_rib, tags.species_and_sequences.Tag, sequence, tags.species_and_sequences.Tagcell, tags.species_and_sequences.Name, tags.species_and_sequences.cod, 
		tags.species_and_sequences.efac, tags.species_and_sequences.as, tags.species_and_sequences.tRNAs, tags.species_and_sequences.codAS, tags.species_and_sequences.t3tRNA, tags.species_and_sequences.AA, 
		tags.species_and_sequences.Tag_constant_vars, tags.species_and_sequences.no_of_tRNAs, tags.species_and_sequences.no_of_variables, tags.species_and_sequences.cod_idx)

		idx_first_codon = length(Tags)+1

		mRNA_concentration = filter(row -> row.gene == seq_name, mRNA).mRNA[]

		# numbers and arrays for rhs_sequence
		num_species 	= length(tags.species_and_sequences.Name)
		len_sequence 	= length(sequence) # Codon length with stop codon 
		rate_RIB        = zeros(Float64, len_sequence-1)
		q_R_elo     	= ones(len_sequence-1)
		M 				= zeros(num_species, len_sequence-1)
		mtRNA 			= zeros(size(M,1),2);

		N_one 	= zeros(len_sequence-1)
		N_two 	= zeros(len_sequence-1)
		K_T3 	= zeros(num_species)
		K_tRNA 	= zeros(num_species) 

		local_tags 		= txtl_types.data_struct_integration(
			tags.scale, tags.idx, species_and_seqs, tags.flags, tags.dyn_eff,
			tags.init_c, tags.p_met, tags.p_seq, tags.p_dyn, seq_name,
			Tags, alpha_one, alpha_two, matrix_N, tRNA_charg, t3_Codon_struc, 
			efac_pos_aa, efac_int_sd, seq_names, idx_first_codon, mRNA_concentration,
			num_species, len_sequence, rate_RIB, q_R_elo, M, 
			mtRNA, N_one, N_two, K_T3, K_tRNA)

		push!(local_tags_list, local_tags)
	end

	return local_tags_list

end


function alpha_struc_generation(tags, sequence)

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
	alpha_struc = zeros(tags.species_and_sequences.no_of_tRNAs,length(sequence)-1)

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
	strCodon = sequence[1:end-1]
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
	O = zeros(size((tags.species_and_sequences.Name),1)-tags.species_and_sequences.no_of_tRNAs, length(sequence)-1)
	alpha_struc = vcat(O, alpha_struc)

	alpha_one = zeros(Int, length(sequence)-1)
	alpha_two = zeros(Int, length(sequence)-1)
	alpha_three = zeros(Int, length(sequence)-1)
	num_col = 0

	for column in eachcol(alpha_struc)
		num_col += 1
		num_element = 0.0

		for (index, element) in enumerate(column)
			if element !=0 && num_element == 0.0
				alpha_one[num_col] = index
				num_element += 1.0

			elseif element !=0 && num_element == 1.0
				alpha_two[num_col] = index
				num_element += 1.0

			elseif element !=0 && num_element == 2.0
				alpha_three[num_col] = index
				num_element += 1.0
			end
		end
	end

	return alpha_one, alpha_two

end

function N_generation(tags, sequence)

	# Function   : generation of N
	# Written by : Maike Kuschel
	# Created on : 15/04/2015
	# Last Update: 25/03/2025 Translation to julia, 10.12.2025 now it gives a vector with indices instead of a sparse matrix
	# Purpose    : generate array to calculate the tRNA (release) and T3 (consumption) concentration  

	# matrix_N : rows = length of codon sequence; columns = number of all
	# possible codons. 
	# array combines the codon sequence with the Codons --> see MOD_EL_rhs
	# (multiplication of matrix_N with C_Codon(concentrations of each codon),
	# generates only entries in matrix_N, if the codons fit with the
	# sequence)

    sequence = sequence[1:end-1]

    matrix_N = Vector{Int32}(undef, length(sequence))

    for i in eachindex(sequence) 

        for j = 1:size(tags.species_and_sequences.cod,1)-1

            if String(sequence[i]) == tags.species_and_sequences.cod[j]
                matrix_N[i] = j
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
	y0_metab_output = y_steady_state[1:length(tags.species_and_sequences.Tag)]                  # All species without the codons
	y0_metab = y0_metab_output[1:end-4]                                                         # Species of the metabolite pool
	y0_output = y0_metab_output[end - 3: end]                                                   # Output species  "RNA" "S70"  "rna"  "PRT"  "PTC"
	ribosome_positions = [zeros(length(seqs[seq_name].Codon) - 1) for seq_name in seq_names]    # All codons per sequence 
	output_sequences = reduce(vcat, ([y0_output; ribosome_positions[i]...] for i in eachindex(ribosome_positions))) # One set of output variables per sequence 
	y0_all = vcat(y0_metab, output_sequences);                                                   # concatenate everything: metabolites pool, output 1, sequence 1, output 2, sequence 2, output 3, sequence 3 ...
	
	return y0_all
end

"""
	include_mRNA_concentrations!(y0_all, tags_pre, seq_names, seqs, mRNA)

26.11.2025
Puts mRNA concentrations (From Lu, Vogel et al. 2007) into states vector.
"""

function include_mRNA_concentrations!(y0_all, tags_pre, seq_names, seqs, mRNA)
    # Fill y0_all with mRNA concentrations

    RNA_positions = calculate_RNA_position(tags_pre, seq_names, seqs)

    for (i, seq_name) in enumerate(seq_names) 
        mRNA_concentration = subset(mRNA, :gene => ByRow(==(seq_name)))
        y0_all[RNA_positions[i]] = mRNA_concentration.mRNA[1]
    end

    return nothing
end

# ------------------------------- Integration -----------------------------------------

function ODE!(du, u, p, t)

    fill!(du, 0.0)  # Reset to zero, just like initializing with zeros

    rhs_metabolites!(du, u, p)
    
    rhs_sequence!(du, u, p)

    nothing 
end 


# 3. Define the ODE function for solving (THREADED VERSION - only change: map -> Threads.@threads)
"""
Explanation of function:

Aim: Enable parallel computing of du for multiple sequences.

1. Create sequence specific states vector (u) 
- Generate numbers to access certain parts of u
- Generate parts of u that need to be assembled per specific to the sequence. This is the metabolite pool and the codons of a specific sequence.

2. Compute rates vector (du)
- Put together metabolite pool and codon vector, initialize du
- Give those sequence specific vectors to the ODE! function 
- Save sequence specific result in list 

3. Construction of one du for all sequences 
- Save all metabolite and codon du vectors in lists
- Summation of du metabolites 
- Combining du metabolites with all du of codons 
"""

function ODE_parallel_distributed!(du, u, p, t)


    # Vectors and numbers 
    @unpack y_codons_lengths, offsets, n_metab, n_seq, du_metab_all, du_codons, y_codons, results, du_seq, du_metab = p.vectors

    # Generate parts of u that need to be assembled specific to the sequence   
    @views u_metab =  u[1:n_metab]				                            # Save all pool species from u in a vector 

    # thread-local accumulator for du[1:start]
    #du_head_threads = [zero(u_metab) for _ in 1:Threads.nthreads()] ---- wrong for servers! - Arnav
	# Jupyter/REPL often runs on an extra "interactive" thread whose id is
    # Threads.nthreads() + 1; allocate one extra slot to avoid BoundsError.
    du_head_threads = [zero(u_metab) for _ in 1:(Threads.nthreads() + 1)]

    Threads.@threads for i in 1:n_seq				                        # Threading instead of map, # ONLY CHANGE: map -> Threads.@threads for  
        tid = Threads.threadid()

        local_tags = p.local_tags_list[i]		                            # Extraction of local tags
        idx = n_metab + offsets[i]                                          # Create index to access u on the sequence specific indexes
        len = y_codons_lengths[i]            
        
        # thread-safe local views
        u_codons = @view u[idx+1 : idx+len]

        # assemble without allocation
        u_seq = similar(u, n_metab + len)
        @views u_seq[1:n_metab] .= u_metab
        @views u_seq[n_metab+1:end] .= u_codons

        du_seq_local = du_seq[i]					                        # Initialization of sequence specific rates vector

        ODE!(du_seq_local, u_seq, local_tags, t)	                        # Give those sequence specific vectors to the ODE! function  

        # thread-local accumulation
        @views du_head_threads[tid] .+= du_seq_local[1:n_metab]               # Save du results from sequence specific ODEs to common rates vector

        @views du[idx+1:idx+len] .= du_seq_local[n_metab+1:end]
    end

    # Add up du for metabolites for all sequences
    du[1:n_metab] .= 0
    for acc in du_head_threads
        @views du[1:n_metab] .+= acc
    end

    return nothing 
end

"""
	matrix_free_KT3(tags, V_TL_EL, c_codon_seq, C_vect)

12.01.2026 Cedric Meyer, Calculates rates vectors K_tRNA and K_T3
"""

function matrix_free_KT3!(tags, V_TL_EL, c_codon_seq, C_vect)
    
    ### Division of tRNA concentration by codon concentration # ---------------------------------------------------------------------

	fill!(tags.N_one, 0.0)
	fill!(tags.N_two, 0.0)
	fill!(tags.K_T3, 0.0)
	fill!(tags.K_tRNA, 0.0)

    # N_three = zeros(tags.len_sequence-1)

    counter = 0
    @inbounds for (index, c_codon) in zip(tags.alpha_one, c_codon_seq)
        counter += 1
        tags.N_one[counter] = C_vect[index] / c_codon
    end

    counter = 0
    @inbounds for (index, c_codon) in zip(tags.alpha_two, c_codon_seq)
        if index != 0.0
            counter += 1
            tags.N_two[counter] = C_vect[index] / c_codon
        elseif index ==0.0
            counter += 1
            tags.N_two[counter] = 0.0
        end
    end

     # ---------------------------------------------------------------------
    @inbounds for (codon_position, tRNA_position) in enumerate(tags.alpha_one) 
        tags.K_T3[tRNA_position] += tags.N_one[codon_position]*V_TL_EL[codon_position]
    end

    @inbounds for (codon_position, tRNA_position) in enumerate(tags.alpha_two)
        if tRNA_position != 0 
            tags.K_T3[tRNA_position] +=tags.N_two[codon_position]*V_TL_EL[codon_position]
        end
    end

    # ---------------------------------------------------------------------

    V_TL_EL_subset  = @view V_TL_EL[3:end] 
    N_one_subset    = @view tags.N_one[2:end-1]
    N_two_subset    = @view tags.N_two[2:end-1]
    alpha_one_subset= @view tags.alpha_one[2:end-1]
    alpha_two_subset= @view tags.alpha_two[2:end-1]

    @inbounds for (codon_position, tRNA_position) in enumerate(alpha_one_subset) 
        tags.K_tRNA[tRNA_position] += N_one_subset[codon_position]*V_TL_EL_subset[codon_position]
    end

    @inbounds for (codon_position, tRNA_position) in enumerate(alpha_two_subset)
        if tRNA_position != 0 
            tags.K_tRNA[tRNA_position] += N_two_subset[codon_position]*V_TL_EL_subset[codon_position]
        end
    end

    return nothing

end 

# ---------------------------- figures -----------------------------------------

"""
	sequence_dataframes(sol, local_tags_list)

Creates a dataframe with shared state variables and sequence specific ribosome concentrations. One dataframe per sequence. Gives back a list with a dataframe per sequence.
"""

function sequence_dataframes(sol, local_tags_list)
    
    df_list = []
    seq_names = local_tags_list[1].seq_names
    n_metab = length(local_tags_list[1].Tags)
	idx = n_metab - 5
    
	for (seq_idx, seq_name) in enumerate(seq_names)
        local_tags = local_tags_list[seq_idx]

		# Prepare data frame 
        df = DataFrame()
        df.time = sol.t

		# Metabolites
        for i in 1:n_metab 
            var_name = i <= length(local_tags.species_and_sequences.Tag) ? local_tags.species_and_sequences.Tag[i] : "State_$i"
            df[!, var_name] = [u[i] for u in sol.u]
        end

		# Output variables 

		out_start 	= n_metab - 4
		out_end		= n_metab 

		L = 0

		for i in out_start:out_end
			var_name = i <= length(local_tags.species_and_sequences.Tag) ? local_tags.species_and_sequences.Tag[i] : "State_$i"
            df[!, var_name] = [u[i + L] for u in sol.u]
		end

		# Ribosome complexes
        L = length(local_tags.species_and_sequences.Seq)-1
        ribo_start = out_end + 1
        ribo_end = out_end + L

        for i in 1:L
            df[!, "Ribosome_$(i)"] = [u[ribo_start-1+i] for u in sol.u]
        end

        push!(df_list, df)

        idx += L + 5

    end

    return df_list
end

"""
	calculate_V_TL_EL(tags, y)

Calculates V_TL_EL just like in rhs_sequence. Needed for plotting the translation elongation rate profile.
"""


function calculate_V_TL_EL(tags, y)

    # Calculate q_R_elo 
    # necessary parameters 
    len_sequence = length((tags.species_and_sequences.Seq))
    gene_sizes = length((tags.species_and_sequences.Seq))-1
    num_species = length(tags.species_and_sequences.Name)

    @unpack µ, k_proc, k_cat_Pol, C_RNP, V_TCmax, K_M_D, k_cat_E, C_RNE, V_degmax, K_M_deg, k_D, n_DNA_t, n_DNA_s = tags.p_seq.TX
    @unpack k_TLEmaxj, KM_T3j, KM_EFGGTP = tags.p_seq.EL
    @unpack k_TE, KM_RK, KM_GTP, KM_RF3GTP, k_RF3_diss_GDP, k_RF1_diss_GTP, k_GTP_hydr, k_RF1_ass = tags.p_seq.TER
    @unpack k_RF3_GDP_ass, k_RF3_GTP_ass, k_RF3_GDP_diss, k_RF3_GTP_diss = tags.p_seq.RF3
    @unpack A1, A2, A3, A4, A5, k_cat_K_M_G1, k_cat_K_M_G2, k_RRF, q_RRF, one_q_G2_k_split = tags.p_seq.PTC 

    @unpack k_IN_70SIC, KM_mRNA, KM_50S, k_IN_IF2D, k_on_30S_IF3, k_off_30SIF3, k_on_30SIF3_IF1, k_off_30SIF3IF1, k_on_30SIF3_IF2, k_off_30SIF3IF2, k_on_30SIF3IF1_IF2, k_off_30SIF3IF1IF2, k_on_30SIF3IF2_IF1, k_off_30SIF3IF2IF1, k_on_30SIF123_fMet, k_off_30SIF123fMet = tags.p_met.IN
    @unpack k_EFGGDP_Diss, k_EFG_GDP_Ass, k_EFG_GTP_Ass, k_EFGGTP_Diss, k_EFTsf, k_EFTsr, Keq_EFTu, KM_GTP, KM_GDP, KM_EFTuGDP, KM_EFTuGTP, Ki_EFTuGDP, Ki_EFTuGTP, c_EFG_t = tags.p_met.EF
    @unpack k_cat, KM_ARS_tRNA, KM_ARS_aa, KM_ARS_ATP = tags.p_met.tRNA
    @unpack k_T3j_f, k_T3j_d = tags.p_met.t3 
    q_RBS = tags.p_met.q_RBS

    # Unpack dynamic variables 
    S30, S50, IF1, IF2, IF3, SC1, SC2, SC3, SC4, PIC, EFG, EFD, EGT, EGD, EFT, RF3, RGT, RGD, b1t, ta1, ta2, tc1, td1, te2, tf1, tg2, tg3, th1, ti1, tk1, t1l, t2l, t3l, t4l, t5l, tm1, tn1, tp1, tp2, tp3, 
    tq1, tq2, tr2, tr3, tr4, tr5, ts1, ts2, ts3, ts5, tt1, tt2, tt4, tv1, tv2, tw1, ty1, tA1, tA2, tC1, tD1, tE2, tF1, tG2, tG3, tH1, tI1, tK1, tL1, tL2, tL3, tL4, tL5, tM1, tN1, tP1, tP2, tP3, tQ1, tQ2, 
    tR2, tR3, tR4, tR5, tS1, tS2, tS3, tS5, tT1, tT2, tT4, tV1, tV2, tW1, tY1, a1t, a2t, c1t, d1t, e2t, f1t, g2t, g3t, h1t, i1t, k1t, l1t, l2t, l3t, l4t, l5t, m1t, n1t, p1t, p2t, p3t, q1t, q2t, r2t, r3t, 
    r4t, r5t, s1t, s2t, s3t, s5t, t1t, t2t, t4t, v1t, v2t, w1t, y1t, RNA, S70, rna, PRT, PTC = y



    # calculation of q_R_elo 
    c_RIB = y[tags.idx_first_codon:end]
    # rna_free = y[tags.idx.RNA]
    rna_free = tags.RNA
    total_RNA = rna_free + sum(c_RIB[1:12])
    s70_ini = y[tags.idx.S70]
    PTC = y[tags.idx.PTC]

    q_R_elo     = ones(gene_sizes);

	q_R_elo[1] = 1 - (sum(c_RIB[1:12]) +s70_ini)/(rna_free+sum(c_RIB[1:12]))
        
    @views for i in 2:(len_sequence - 13)
        q_R_elo[i] 	= 1 - c_RIB[i+12] / (total_RNA-sum(c_RIB[i:i+11]));
    end 


    @views for i in (len_sequence-13) : (len_sequence-1)
            q_R_elo[i] = 1 - (c_RIB[end] + PTC) / (total_RNA - sum(c_RIB[i:end-1]))
    end

    # dynamic efficiency factor 
    if tags.flags.dynamic_efficiency == true
        efac = dynEffficiencyCalculation(tags, EGT, y);
    end

    # V_EL
    C_vect = reshape(y[1:num_species], 1, num_species)  # Directly assign
    c_Codon = tags.t3_Codon_struc * C_vect'

    V_TL_EL	= zeros(len_sequence-1);

    EGT_inv = inv(EGT)


    for i in 1:len_sequence-2

        cod_next = tags.species_and_sequences.Seq[i+1]
        cod_curr = tags.species_and_sequences.Seq[i]

        j = tags.species_and_sequences.cod_idx[cod_next]    # O(1) Lookup
        k = tags.species_and_sequences.cod_idx[cod_curr]    # O(1) Lookup

        cc = c_Codon[j]

        eff = efac[k] * min(tags.efac_pos_aa[i], tags.efac_int_sd[i])
        denom = 1 + KM_T3j/cc + KM_EFGGTP*EGT_inv
        V_TL_EL[i+1] = q_R_elo[i] * k_TLEmaxj * c_RIB[i] * eff / denom

    end

    return V_TL_EL

end

"""
	V_TL_EL_smooth(tags, V_TL_EL, c_RIB)

V_TL_EL divided by ribosome concentration and smoothed with sliding window. Minutes are transformed to seconds.
"""


function V_TL_EL_smooth(tags, V_TL_EL, c_RIB)
	# parameters and arrays 
	len_sequence 	= length(tags.species_and_sequences.Seq) # Codon length with stop codon
	L 				= len_sequence - 2
	V_ELO 			= zeros(L)
	V_EL_smooth  	= zeros(L)
	h1				= zeros(L)

	# V_TL_EL per second
	V_ELO = (V_TL_EL[2:end]./60)./(c_RIB[1:end-1]);
	V_EL_smooth[1:3] 	= [mean(V_ELO[1:i+3])   for i in 1:3]
	V_EL_smooth[4:L-3] 	= [mean(V_ELO[i-3:i+3]) for i in 4:L-3]
	V_EL_smooth[L-2:L] 	= [mean(V_ELO[i-2:L])   for i in L-2:L];

	# Plot
	h1 .= mean(V_ELO)

	return V_EL_smooth, h1

end

"""
	calculate_translation_rate(df)

Calculates translation rate.
"""

function translation_rate_m(df)
    x1 = df.time[end]
    x2 = df.time[end-1]
    y1 = df.PRT[end]
    y2 = df.PRT[end-1]

    # Calculate increase 
    m = (y2 - y1) / (x2 - x1)

    return m 

end

"""
	protein_production(sol, tags)

16/05/2025 Cedric Meyer 
Calculates the translation rate at timepoint 100
"""

function translation_rate_ODE(sol, p)
	tags = p.local_tags_list[1]
	# Zustand zum Zeitpunkt t
	t_query = 100.0
	u_t = sol(t_query)

	# Raten zum Zeitpunkt t berechnen
	du_t = similar(u_t)  # Speicher für Ableitungen
	ODE_parallel_distributed!(du_t, u_t, p, t_query)

	return Float64(du_t[tags.idx.PRT])
end

function save_results(sol, local_tags_list, save_data)
	
	if save_data == true
		# Save results 
		timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM")
		filename = "output/sol_" * timestamp * ".jld2"

		@save filename sol local_tags_list

		println("Solution saved at $filename.")
	end
end

function save_translation_rates(translation_rates)
	
	# Save results 
	timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM")
	filename = "../data/output/translation_rates_" * timestamp * ".jld2"

	@save filename translation_rates

	println("Solution saved at $filename.")

end

# --------------------------- mRNA competition -----------------------------------------

"""
	calculate_RNA_position(tags, seq_names, seqs)
	
04.09.2025. Cedric Meyer. Calculates the position of RNA in y0_all states vector. Needed to set mRNA concentrations from file.
"""

function calculate_RNA_position(tags, seq_names, seqs)

	n_metab 	= tags.idx.PTC
	n_ribo 		= Vector{Int}()

	for seq_name in seq_names
		push!(n_ribo, length(seqs[seq_name].Codon)-1)
	end

	n_RNA 	= Vector{Int}()
	idx 	= n_metab -4

	for sequence_length in n_ribo 
		idx += 1
		push!(n_RNA, idx)
		idx += 4
		idx += sequence_length
	end

	return n_RNA
end

"""
	function to_DataFrame(mat)

11/07/2025 Cedric Meyer 
Conversion of a matrix to a DataFrame with header
"""

function to_DataFrame(mat)
	# Erste Zeile als Spaltennamen extrahieren
	header = Symbol.(mat[1, :])  # in Symbole umwandeln (empfohlen für DataFrame)

	# Restliche Zeilen sind die Daten
	data = mat[2:end, :]

	# DataFrame erzeugen
	df = DataFrame( data , header)
	return df 
end

# ---------------------- Model validation ----------------------------------

"""
	call_model(sequence_file, heterologous_genes, num_background_genes)

01.12.2025; Model wrapped in a function to call it with one line.
"""

function call_model(sequence_file, heterologous_genes, num_background_genes)

	include("parameters.jl")
	
	hybrid_file		    = "lacS.txt"
	seqs_to_use 	    = (1,2401)                                 # Range of sequences in Fasta file. Takes all sequences from sequence a to sequence b

	# initialization of data containers
	scale = txtl_types.init_scale()         # struct that makes it possible to influence rates by scaling them

	µ = 1                    # A growth rate (0.1 to 1.8) can be chosen and the proteins of the TXTL machinery will be tuned according to that growth rate and the processivity of the RNA polymerase

	## flags for efficiency factors
	positive_aminoacid	= true 
	anti_shinedalgarno	= false      # If no hybridization file is available put flag to false. We apply the intaRNA program's intarna function (q: aSD, t: sequence, SeedBP: 5, out: tMinE) to generate a hybridization profile along the sequence. Link to the website: https://rna.informatik.uni-freiburg.de/IntaRNA/Input.jsp;jsessionid=2AE257E91E7276283F22CD77A5A8AA03
	dynamic_efficiency	= true 
	save_data           = false

	flags = txtl_types.flags(positive_aminoacid, anti_shinedalgarno, dynamic_efficiency)

	# Data import --------
	mRNA 		= to_DataFrame(readmatrix("../data/input/mRNA_concentrations.xlsx", "mRNA_concentration", "A1:B2402"))
	mRNA.mRNA 	= Float64.(mRNA.mRNA)

	background_genes = String.(first(mRNA.gene, num_background_genes))
	genes = vcat(heterologous_genes, background_genes)
	println("Number of sequences: $(length(genes))") 

	mRNA = subset(mRNA, :gene => ByRow(in(genes)))

	## Import information about model sp ecies and sequences 
	seqs, seq_names, hyb_file_dict, species_and_seqs, idx = species_and_sequences(genes, µ, sequence_file, hybrid_file, flags, seqs_to_use)

	## Calculate efficiency factors 
	efac_pos_aa_dict, efac_int_sd_dict, dynamic_efficiencies = efficiency_factor_initialization(seqs, seq_names, hyb_file_dict, species_and_seqs, flags)

	## calculation of initial concentrations
	y0, initial_concentrations = calculation_of_initial_concentration(species_and_seqs, seqs, seq_names, scale, idx);

	## calculate and load kinetic parameters
	p_met, p_seq, p_dyn = parameters(initial_concentrations, species_and_seqs, scale, µ);

	## Release factor balance
	RF_balance!(y0, initial_concentrations)

	# Preintegration
	## Definition of mass matrix
	conserved = [idx.EFG, idx.EFD] # EFG, EFTu_GDP - Conserved species, Algebraic equations are given by mass balances

	M = Matrix{Float64}(I, length(y0), length(y0)) # Creates unity matrix in the size of y0

	M[conserved, :] .= 0; # setting mass matrix to zero to  remove time derivatives from those variables. This converts these differential equations into algebraic constraints


	tags_pre = txtl_types.data_struct_preintegration(scale, idx, species_and_seqs, flags, dynamic_efficiencies, initial_concentrations, p_met, p_seq, p_dyn, mRNA.mRNA[1])

	## Integration of the DAE system without DNA in the system 
	y_steady_state = preintegration(tags_pre, y0, M);


	# Elongation factor balance 
	EF_balance!(tags_pre, y0, y_steady_state)


	# Build states vector for multiple sequences 
	y0_all = combined_states_vector(tags_pre, y_steady_state, seqs, seq_names)

	# Prepare per-sequence tags (for use in ODE)
	local_tags_list = generate_local_tags(seqs, seq_names, tags_pre, efac_pos_aa_dict, efac_int_sd_dict, mRNA);


	# Integration
	# PERFORMANCE FIX: Use Diagonal mass matrix instead of dense matrix
	# For 100 sequences with ~22k states, dense matrix would be ~4GB and cause O(n³) scaling
	n_all = length(y0_all)
	mass_diag_all = ones(Float64, n_all)
	mass_diag_all[conserved] .= 0.0
	M = Diagonal(mass_diag_all)  # Much more memory-efficient: O(n) instead of O(n²)

	tspan 	= (0.0, 100.0) ## Define time span for the simulation

	p       = generate_p(local_tags_list)

	f 		= ODEFunction(ODE_parallel_distributed!, mass_matrix = M) ## Define ODEFunction to make use of mass matrix

	prob 	= ODEProblem(f, y0_all, tspan, p); ## Create an ODE problem

	# --- Timing the integration step ---
	# println("\n⏳ Integrating ODE system (distributed, test set)...")

	sol     = solve(prob, FBDF(autodiff=false), reltol = 1e-8, abstol = 1e-6); ## Solve ODE system

	println("Integration completed!")

	translation_rate = translation_rate_ODE(sol, p)


	return translation_rate
end



"""
	calculate_absolute_protein_µgmL(c_x0, µ, t_process, weight_cell, volume_cell, translation_rate, mass_protein)

30.10.2025, Function calculates the protein concentration at the end of the cultivation. Given translation rates and all other parameters specified in function.
"""
function calculate_absolute_protein_µgmL(c_x0, µ, t_process, weight_cell, volume_cell, translation_rate, mass_protein) # give mass_protein in Dalton, t in hours, everything else in SI units 
	
	qp 				= calculate_qp(translation_rate, weight_cell, volume_cell) 	# µmol/gDW/h
	c_p 			= qp *c_x0 / µ *(exp(µ*t_process) -1) 						# µM
	mass_protein 	= mass_protein/1000 										# g/µmol
	c_gL 			= c_p * mass_protein										# g/L -> mg/mL
	c_µgmL			= c_gL * 1000												# µg/mL
	
	return c_µgmL
end


"""
	calculate_qp(translation_rate)

Cedric Meyer, 01/08/2025
Transforms the translation rate that is calculated in the model into biomass specific production rate qp
"""

function calculate_qp(translation_rate, weight_cell, volume_cell)
	cells_per_gram 		= 1/weight_cell						# 1/g
	volume_per_gram		= cells_per_gram*volume_cell		# L_cell/gDW 
	v_t 				= translation_rate*60				# µmol/L_cell/h
	qp  				= v_t * volume_per_gram				# µmol/gDW/h
	return qp 
end