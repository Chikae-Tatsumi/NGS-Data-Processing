# Necessary to do in one directry (put OTU(ASV) table and fasta file at Database folder)

library(Tax4Fun2)
library(seqinr)
setwd("~/R/Database/Tax4Fun2")
dir.create("Tax4Fun2")

#Step 2: Generate your own reference datasets
#1. Extracting SSU seqeunces (16S rRNA and 18S rRNA)
extractSSU(genome_file = "OneProkaryoticGenome.fasta", file_extension = "fasta", path_to_reference_data ="Tax4Fun2_ReferenceData_v2")
#2. Assigning functions to prokayotic genomes
assignFunction(genome_file = "OneProkaryoticGenome.fasta", file_extension = "fasta", path_to_reference_data = "Tax4Fun2_ReferenceData_v2", num_of_threads = 1, fast = TRUE)
#3. Generate the reference data
generateUserData(path_to_reference_data ="Tax4Fun2_ReferenceData_v2", path_to_user_data = ".", name_of_user_data = "User_Ref0", SSU_file_extension = "_16SrRNA.ffn", KEGG_file_extension = "_funPro.txt")

#Step 3: Making functional predictions
#1. Making functional predictions using the default reference data only
#1. Run the reference blast
runRefBlast(path_to_otus ="12_otus_97.fasta" , path_to_reference_data ="Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "Tax4Fun2", database_mode = "Ref99NR", use_force = T, num_threads = 6)
# 2) Predicting functional profiles
# Remove the first row's # and the second row's # of "17_otu_97_table_taxonomy.txt" --> "17_otu_97_table_taxonomy_tax4fun.txt"
makeFunctionalPrediction(path_to_otu_table = "17_otu_97_table_taxonomy_tax4fun.txt", path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "Tax4Fun2", database_mode = "Ref99NR", normalize_by_copy_number = TRUE, min_identity_to_reference = 0.97, normalize_pathways = FALSE)
# note. normalize_pathways = FALSE will affiliate the rel. abundance of each KO to each pathway it belongs to. By setting it to true, the rel. abundance is equally distributed to all pathways it was assigned to.)

#Step 4: Calculating (multi-)functional redundancy indices (experimental)
calculateFunctionalRedundancy(path_to_otu_table = "17_otu_97_table_taxonomy_tax4fun.txt", path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "Tax4Fun2", database_mode = "Ref99NR", min_identity_to_reference = 0.97)

# Don't forget to move the output file (Tax4Fun2) for your project file
