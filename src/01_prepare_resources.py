# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 10:44:29 2024

@author: cartama

    The script preparares resources used in OncoVI. 
    
"""

# --------------------------------------------------------------------------- #

# Libraries
import pandas as pd
import os  
import json

# --------------------------------------------------------------------------- #

# Base directory
base_dir = "V:/gruppen/AG_Bioinfo/members/Carta/PhD_Project/Guidelines_implementation/Preparation_for_submission"

# Directory for original resources (i.e., resources not modified)
orig_resdir = os.path.join(base_dir, "02_original_resources")

# Directory for final resources 
final_resdir = os.path.join(base_dir, "03_resources")


# --------------------------------------------------------------------------- #
# COSMIC RESOURCES ---------------------------------------------------------- # 

# Read in the cancer gene census (CGC) csv file in a dataframe 
# (The dataframe was downloaded as .csv format from here:
#  https://cancer.sanger.ac.uk/census . Accessed 1 May 2024.)

cgc_df = pd.read_csv(os.path.join(orig_resdir,"Census_allWed_May_1_17_27_56_2024.csv"),
                     sep=',')

# Select genes with "Tier" == 1

cgc_tier1 = cgc_df.loc[cgc_df['Tier'] == 1] 
# 581 genes from CGC Tier 1 

# Select genes from Tier 1 with "Role in Cancer" containing TSGs 
# (i.e. tumor supressor genes)

tsg_tier1 = cgc_tier1[cgc_tier1['Role in Cancer'].str.contains('TSG', na=False)] 
# 273 TSGs

# Saving the dataframe as csv file 
tsg_tier1.to_csv(os.path.join(final_resdir,"tsg_tier1.csv"), sep='\t', index=False)


# --------------------------------------------------------------------------- #
# OncoKB RESOURCES ---------------------------------------------------------- # 

# Read in the OncoKB Cancer Genes tsv file in a dataframe 
# (The dataframe was downloaded as .tsv format from here:
#  https://www.oncokb.org/cancer-genes . Accessed 29 April 2024.)

cancer_genes_oncokb_df = pd.read_csv(os.path.join(orig_resdir,"cancerGeneList.tsv"),
                                     sep='\t')

# Select genes with "Is Tumor Suppressor Gene" equal to "Yes"

tsg_oncokb = cancer_genes_oncokb_df[cancer_genes_oncokb_df['Is Tumor Suppressor Gene'] == "Yes"] 
# 364 TSGs

# Saving the dataframe as csv file 
tsg_oncokb.to_csv(os.path.join(final_resdir,"tsg_oncokb.csv"), sep='\t', index=False)

# --------------------------------------------------------------------------- #

## Define bona fide tumor supressor genes for OVS1 criterion
## (i.e., genes that have TSG role and Tier 1 in the COSMIC 
##  cancer gene census or in OncoKB Cancer Genes)

# Convert tsg_tier1 to list
tsg_tier1_list = tsg_tier1["Gene Symbol"].to_list()

# Convert tsg_oncokb to list
tsg_oncokb_list = tsg_oncokb["Hugo Symbol"].to_list()

# Merge the two list and extract unique genes 
bona_fide_tsg = list(set(tsg_tier1_list) | set(tsg_oncokb_list))
# 467 bona fide tumor supressor genes

# Saving the dataframe as txt file
with open(os.path.join(final_resdir,"bona_fide_tsg.txt"), 'w') as file:
    
    # Join the list elements into a single string with a newline character
    data_to_write = '\n'.join(bona_fide_tsg)
    
    # Write the data to the file
    file.write(data_to_write)
    
# --------------------------------------------------------------------------- #
# COSMIC RESOURCES ---------------------------------------------------------- # 

# Select genes with "Role in Cancer" containing oncogene
oncogenes = cgc_df[cgc_df['Role in Cancer'].str.contains('oncogene', na=False)]
# 320 oncogenes

# Saving the dataframe as csv file 
oncogenes.to_csv(os.path.join(final_resdir,"oncogenes_cgc.csv"), sep='\t', index=False)

# --------------------------------------------------------------------------- #
# OncoKB RESOURCES ---------------------------------------------------------- # 

# Select genes with "Is Oncogene" equal to "Yes"

og_oncokb = cancer_genes_oncokb_df[cancer_genes_oncokb_df['Is Oncogene'] == "Yes"] 
# 412 OGs

# Saving the dataframe as csv file 
og_oncokb.to_csv(os.path.join(final_resdir,"og_oncokb.csv"), sep='\t', index=False)

# --------------------------------------------------------------------------- #

## Define oncogenes (i.e., genes that have oncogene role 
## in the COSMIC cancer gene census or in OncoKB Cancer Genes)

# Convert oncogenes to list
oncogenes_list = oncogenes["Gene Symbol"].to_list()

# Convert og_oncokb to list
og_oncokb_list = og_oncokb["Hugo Symbol"].to_list()

# Merge the two list and extract unique genes 
ogs_list = list(set(oncogenes_list) | set(og_oncokb_list))
# 553 oncogenes

# Saving the dataframe as txt file
with open(os.path.join(final_resdir,"ogs_list.txt"), 'w') as file:
    
    # Join the list elements into a single string with a newline character
    data_to_write = '\n'.join(ogs_list)
    
    # Write the data to the file
    file.write(data_to_write)
    
# --------------------------------------------------------------------------- #
# COSMIC RESOURCES ---------------------------------------------------------- # 

# Select genes with "Role in Cancer" containing TSG
tsg_cgc = cgc_df[cgc_df['Role in Cancer'].str.contains('TSG', na=False)]
# 320 TSGs

# Saving the dataframe as csv file 
tsg_cgc.to_csv(os.path.join(final_resdir,"tsg_cgc.csv"), sep='\t', index=False)

# --------------------------------------------------------------------------- #

## Define tumor supressor genes (i.e., genes that have TSG role 
## in the COSMIC cancer gene census or in OncoKB Cancer Genes)

# Convert tsg_cgc to list
tsg_cgc_list = tsg_cgc["Gene Symbol"].to_list()

# Merge the two list and extract unique genes 
tsg_list = list(set(tsg_cgc_list) | set(tsg_oncokb_list))
# 498 TSGs

# Saving the dataframe as txt file
with open(os.path.join(final_resdir,"tsg_list.txt"), 'w') as file:
    
    # Join the list elements into a single string with a newline character
    data_to_write = '\n'.join(tsg_list)
    
    # Write the data to the file
    file.write(data_to_write)
    

# --------------------------------------------------------------------------- #
# COSMIC RESOURCES ---------------------------------------------------------- # 

# Read in a reduction of the COSMIC Census Genes Mutations txt file in a dataframe 
# (The dataframe was originally downloaded as .tar format from here:
#  https://cancer.sanger.ac.uk/cosmic/download/cosmic/v99/mutantcensus .
#  Accessed 09 April 2024.)

# The original file contains "All coding mutations in genes listed in the Cancer
# Gene Census.

# The archive was extracted and the correspondent 
# Cosmic_MutantCensus_v99_GRCh38.tsv.gz obtained. The file was then unzipped:
# Cosmic_MutantCensus_v99_GRCh38.tsv. 

# All these files can be found here:
# W:\pathologie\bioinfo-archive\Giulia_PhD\Guidelines_resources\COSMIC.

# The reduced version of the .tsv file was obtained by filtering columns:
# - GENE_SYMBOL
# - MUTATION_CDS 
# - MUTATION_AA
# - HGVSG

cosmic_df_dir = os.path.join(orig_resdir, "Cosmic_MutantCensus_v99_GRCh38_red.txt")

# Open the file as text file
file_cosmic = open(cosmic_df_dir, "r")

# Each line is saved in a list
lines = file_cosmic.readlines()

# Pandas dataframe 
cosmic_df = pd.DataFrame(lines, columns = ['Data'])

# Split column into multiple named columns
cosmic_df[['GENE_SYMBOL', 'MUTATION_CDS', 'MUTATION_AA', 'HGVSG']] = cosmic_df['Data'].str.split(' ', expand=True)

# Remove the original column
cosmic_df = cosmic_df.drop('Data', axis=1)

# Remove string from last column
cosmic_df['HGVSG'] = cosmic_df['HGVSG'].str.replace('\n', '')

# Remove the header
cosmic_df = cosmic_df.drop(0)

# Print the header
cosmic_df.head()

#  GENE_SYMBOL MUTATION_CDS MUTATION_AA             HGVSG
#1       EP300    c.5350C>T    p.Q1784*  22:g.41177061C>T
#2       EP300    c.3566C>T    p.A1189V  22:g.41158476C>T
#3       EP300    c.4145G>A    p.G1382D  22:g.41168840G>A
#4       EP300    c.3256A>G    p.I1086V  22:g.41155108A>G
#5       EP300    c.4154G>T    p.C1385F  22:g.41168849G>T

# Create a dictionary in which the key is "HGVSG", the value is the
# number of samples in COSMIC with that key
cosmic_hgvsg_dict = {}
        
cosmic_df_hgvsg = cosmic_df

# Remove rows with empty HGVSG 
cosmic_df_hgvsg_valid = cosmic_df_hgvsg[cosmic_df_hgvsg['HGVSG'] != '']

# Create the dictionary
cosmic_hgvsg_dict = cosmic_df_hgvsg_valid.groupby(['HGVSG'])['MUTATION_AA'].apply(list).to_dict()

# Save the dictionary to file 
with open(os.path.join(final_resdir,"cosmic_hgvsg_dictionary.txt"), 'w') as fp:
    json.dump(cosmic_hgvsg_dict, fp)

# --------------------------------------------------------------------------- #

# Read in a reduction of the COSMIC Cancer Mutation Census tsv file in a dataframe 
# (The dataframe was originally downloaded as .tar format from here:
#  https://cancer.sanger.ac.uk/cosmic/download/cancer-mutation-census/v99/alldata-cmc .
#  Accessed 07 May 2024.)

# The original file contains "all coding somatic mutations collected by COSMIC"

# The archive was extracted and the correspondent 
# CancerMutationCensus_AllData_v99_GRCh38.tsv.gz obtained. 
# The file was then unzipped:
# CancerMutationCensus_AllData_v99_GRCh38.tsv. 

# All these files can be found here:
# W:\pathologie\bioinfo-archive\Giulia_PhD\Guidelines_resources\COSMIC.

# The reduced version of the .tsv file was obtained by filtering columns:
# - GENE_NAME
# - Mutation CDS 
# - Mutation AA
# - AA_MUT_START
# - Mutation genome position GRCh38
# - COSMIC_SAMPLE_MUTATED

cosmic_all_data_dir = os.path.join(orig_resdir, "CancerMutationCensus_All_v99_GRCh38_red.tsv")

# Open the file as text file
file_cosmic_all_data = open(cosmic_all_data_dir, "r")

# Each line is saved in a list
lines_all_data = file_cosmic_all_data.readlines()

# Pandas dataframe 
cosmic_all_data = pd.DataFrame(lines_all_data, columns = ['Data'])

# Split column into multiple named columns
cosmic_all_data[['GENE', 'MUT_CDS', 'MUT_AA','AA_REF','Mut_genom_pos_GRCh38', 'COSMIC_SAMPLE_MUTATED']] = cosmic_all_data['Data'].str.split('\t', expand=True)

# Remove the original column
cosmic_all_data = cosmic_all_data.drop('Data', axis=1)

# Remove string from the last column
cosmic_all_data['COSMIC_SAMPLE_MUTATED'] = cosmic_all_data['COSMIC_SAMPLE_MUTATED'].str.replace('\n', '')

# Remove the header
cosmic_all_data = cosmic_all_data.drop(0)

# Print the header
cosmic_all_data.head()

#   GENE    MUT_CDS  ... Mut_genom_pos_GRCh38 COSMIC_SAMPLE_MUTATED
#1  PODN   c.888C>A  ...  1:53077690-53077690                     1
#2  PODN   c.631C>T  ...  1:53075877-53075877                     2
#3  PODN   c.936G>A  ...  1:53077738-53077738                     1
#4  PODN  c.1137del  ...  1:53078503-53078503                     1
#5  PODN  c.1234G>A  ...  1:53078600-53078600                     1

# 1) How many variants with "Mut_genom_pos_GRCh38" == None?
# In case of None values instead of NaN replace() cannot be used to replace
# None with empty strings

mask = cosmic_all_data.applymap(lambda x: x is None)
cols = cosmic_all_data.columns[(mask).any()]
for col in cosmic_all_data[cols]:
    cosmic_all_data.loc[mask[col], col] = ''
    
# Variants with valid genomic position
start_valid = cosmic_all_data[cosmic_all_data['Mut_genom_pos_GRCh38'] != ''] 
# 5,170,495

# Variants with non valid genomic position
start_not_valid = cosmic_all_data[cosmic_all_data['Mut_genom_pos_GRCh38'] == ''] 
# 55,319

# Variants with description on the protein level
mutAA_valid = cosmic_all_data[cosmic_all_data['MUT_AA'] != ''] 
# all variants have a mut AA

# Variants with description on the cDNA level
mutCDS_valid = cosmic_all_data[cosmic_all_data['MUT_CDS'] != ''] 
# all variants have a mut CDS

# Variants with a reference AA
refAA_valid = cosmic_all_data[cosmic_all_data['AA_REF'] != ''] 
# all variants have a ref AA

# Add a new column to the refAA_valid dataframe
refAA_valid = refAA_valid.assign(MUT_AA_mod = None)

# Change nomenclature of synonymous variants (whose AA change they cause is the
# same as the reference AA)
for index, row in refAA_valid.iterrows():
    
    AA_change = refAA_valid.loc[index, "MUT_AA"]
    
    if "=" in AA_change:
        
        AA_change_mod = AA_change.replace("=", AA_change[2])
        refAA_valid.loc[index, "MUT_AA_mod"] = AA_change_mod
        
    else:
        
        refAA_valid.loc[index, "MUT_AA_mod"] = AA_change
        
# The dictionary is populated with variants with a genomic position GRCh38 
cosmic_all_dict = {}
        
for index, row in refAA_valid.iterrows():
    
    gene = refAA_valid.loc[index, "GENE"]
    refAA = refAA_valid.loc[index, "AA_REF"]
    variant = refAA_valid.loc[index, "MUT_AA_mod"]
    variant_c = refAA_valid.loc[index, "MUT_CDS"]
    samples = refAA_valid.loc[index, "COSMIC_SAMPLE_MUTATED"]
    
    
    if gene not in cosmic_all_dict.keys():
        
        cosmic_all_dict[gene] = {}
        
        if refAA not in cosmic_all_dict[gene].keys():
            
            cosmic_all_dict[gene][refAA] = {}
            
            if variant not in cosmic_all_dict[gene][refAA].keys():
                
                cosmic_all_dict[gene][refAA][variant] = {}
                
                if variant_c not in cosmic_all_dict[gene][refAA][variant].keys():
                    
                    cosmic_all_dict[gene][refAA][variant][variant_c] = ""
                    cosmic_all_dict[gene][refAA][variant][variant_c] = int(samples)
            
            else:
                if variant_c not in cosmic_all_dict[gene][refAA][variant].keys():
                    
                    cosmic_all_dict[gene][refAA][variant][variant_c] = ""
                    cosmic_all_dict[gene][refAA][variant][variant_c] = int(samples)
    
        else:
            if variant not in cosmic_all_dict[gene][refAA].keys():
                
                cosmic_all_dict[gene][refAA][variant] = {}
                
                if variant_c not in cosmic_all_dict[gene][refAA][variant].keys():
                    
                    cosmic_all_dict[gene][refAA][variant][variant_c] = ""
                    cosmic_all_dict[gene][refAA][variant][variant_c] = int(samples)
            
            else:
                if variant_c not in cosmic_all_dict[gene][refAA][variant].keys():
                    
                    cosmic_all_dict[gene][refAA][variant][variant_c] = ""
                    cosmic_all_dict[gene][refAA][variant][variant_c] = int(samples)
    else:
        
        if refAA not in cosmic_all_dict[gene].keys():
            
            cosmic_all_dict[gene][refAA] = {}
            
            if variant not in cosmic_all_dict[gene][refAA].keys():
                
                cosmic_all_dict[gene][refAA][variant] = {}
                
                if variant_c not in cosmic_all_dict[gene][refAA][variant].keys():
                    
                    cosmic_all_dict[gene][refAA][variant][variant_c] = ""
                    cosmic_all_dict[gene][refAA][variant][variant_c] = int(samples)
            
            else:
                
                if variant_c not in cosmic_all_dict[gene][refAA][variant].keys():
                    
                    cosmic_all_dict[gene][refAA][variant][variant_c] = ""
                    cosmic_all_dict[gene][refAA][variant][variant_c] = int(samples)
        else:
            
            if variant not in cosmic_all_dict[gene][refAA].keys():
                
                cosmic_all_dict[gene][refAA][variant] = {}
                
                if variant_c not in cosmic_all_dict[gene][refAA][variant].keys():
                    
                    cosmic_all_dict[gene][refAA][variant][variant_c] = ""
                    cosmic_all_dict[gene][refAA][variant][variant_c] = int(samples)
                
            else:
                
                if variant_c not in cosmic_all_dict[gene][refAA][variant].keys():
                    
                    cosmic_all_dict[gene][refAA][variant][variant_c] = ""
                    cosmic_all_dict[gene][refAA][variant][variant_c] = int(samples)

# Save the dictionary to file 
with open(os.path.join(final_resdir,"cosmic_all_dictionary.txt"), 'w') as fp:
    json.dump(cosmic_all_dict, fp)
    

# --------------------------------------------------------------------------- #
# MutSpliceDB RESOURCE ------------------------------------------------------ # 

# Read in the MutSpliceDB csv file in a dataframe 
# (The dataframe was downloaded as .csv format from here:
#  https://brb.nci.nih.gov/cgi-bin/splicing/splicing_main.cgi . 
#  Accessed 06 February 2024.)

mut_splice_db = pd.read_csv(os.path.join(orig_resdir, "MutSpliceDB_BRP_2024-02-06.csv"),
                            sep=',')

# Rename columns
mut_splice_db = mut_splice_db.rename(columns={'Gene Symbol': 'Symbol',
                                              'Entrez Gene ID': 'Entrez_Gene_ID',
                                              'Allele registry ID': 'Allele_registry_ID',
                                              'Splicing effect': 'Splicing_effect',
                                              'RNA-seq evidence': 'RNA-seq_evidence'})

mut_splice_dict = {}

for index, row in mut_splice_db.iterrows():
    
    gene = mut_splice_db.loc[index, "Symbol"]
    mutation = mut_splice_db.loc[index, "Mutation"]
    transcript, dna_change = mutation.split(":")
    transcript = transcript.split(".")[0]
    
    if gene not in mut_splice_dict.keys():
        
        mut_splice_dict[gene] = {}
        
        if transcript not in mut_splice_dict[gene].keys():
            
            mut_splice_dict[gene][transcript] = []
            mut_splice_dict[gene][transcript].append(dna_change)
            
        else:
            
            mut_splice_dict[gene][transcript].append(dna_change)
        
    else:
        if transcript not in mut_splice_dict[gene].keys():
            
            mut_splice_dict[gene][transcript] = []
            mut_splice_dict[gene][transcript].append(dna_change)
            
        else:
            
            mut_splice_dict[gene][transcript].append(dna_change)

# Saving the dictionary as txt file
with open(os.path.join(final_resdir,"mut_splice_dict.txt"), 'w') as fp:
     json.dump(mut_splice_dict, fp)


# --------------------------------------------------------------------------- #
# Amino acid conversion RESOURCE -------------------------------------------- # 

# Read in a table for amino acid conversion from three letter convention to 
# one letter convention.

amino_conversion = pd.read_csv(os.path.join(orig_resdir, "amino_conversion.txt"), 
                               sep='\t', header=0)

# Convert dataframe to dictionary
amino_dict = amino_conversion.set_index('Three_letters').One_letter.to_dict()

# Saving the dictionary as txt file
with open(os.path.join(final_resdir,"amino_dict.txt"), 'w') as fp:
    json.dump(amino_dict, fp)

# --------------------------------------------------------------------------- #
# cancerhotspots.org RESOURCE -------------------------------------------- # 

# Read in the cancerhotspots csv file in a dataframe 
# (The dataframe was downloaded as .txt format from here:
#  https://www.cancerhotspots.org/#/home . 
#  Accessed 26 September 2023.)

cancerhotspots = pd.read_csv(os.path.join(orig_resdir, "hotspots.txt"), sep='\t')

# Extract single residues
single_residue = cancerhotspots[cancerhotspots['Type'] == 'single residue']
# 1110

# Dictionary for single residues
single_residue_dict = {}

for index, row in single_residue.iterrows():
    
    gene = single_residue.loc[index, "Gene"]
    residue = single_residue.loc[index, "Residue"]
    variants = single_residue.loc[index, "Variants"]
    
    if gene not in single_residue_dict.keys():
        
        single_residue_dict[gene] = {}
        
        if residue not in single_residue_dict[gene].keys():
            
            single_residue_dict[gene][residue] = {}
            
            variant_residue_list = variants.split('|')
            
            di = {}
            for var in variant_residue_list:
                key, value = var.split(':')
                
                # if key == "*":
                #     key = "Ter"

                if key not in di.keys():
                    di[key] = ''
                    di[key] = value
                    
                else:
                    di[key] = value
            
            single_residue_dict[gene][residue] = di
            
        
    else:
        if residue not in single_residue_dict[gene].keys():
            
            single_residue_dict[gene][residue] = {}
            
            variant_residue_list = variants.split('|')
            
            di = {}
            for var in variant_residue_list:
                key, value = var.split(':')
                
                # if key == "*":
                #     key = "Ter"
                
                if key not in di.keys():
                    di[key] = ''
                    di[key] = value
                    
                else:
                    di[key] = value
            
            single_residue_dict[gene][residue] = di
            
# Saving the dictionary as txt file
with open(os.path.join(final_resdir,"single_residue_dict.txt"), 'w') as fp:
    json.dump(single_residue_dict, fp)     
            
# Extract in-frame indel
in_frame_indel = cancerhotspots[cancerhotspots['Type'] == 'in-frame indel']
# 55

# Dictionary for in-frame indel 
in_frame_indel_dict = {}

for index, row in in_frame_indel.iterrows():
    
    gene = in_frame_indel.loc[index, "Gene"]
    residue = in_frame_indel.loc[index, "Residue"]
    variants = in_frame_indel.loc[index, "Variants"]
    
    if gene not in in_frame_indel_dict.keys():
        
        in_frame_indel_dict[gene] = {}
        
        if residue not in in_frame_indel_dict[gene].keys():
            
            in_frame_indel_dict[gene][residue] = {}
            
            variant_residue_list = variants.split('|')
            
            dii = {}
            for var in variant_residue_list:
                key, value = var.split(':')
                
                if key not in dii.keys():
                    dii[key] = ''
                    dii[key] = value
                    
                else:
                    dii[key] = value
            
            in_frame_indel_dict[gene][residue] = dii
            
        
    else:
        if residue not in in_frame_indel_dict[gene].keys():
            
            in_frame_indel_dict[gene][residue] = {}
            
            variant_residue_list = variants.split('|')
            
            dii = {}
            for var in variant_residue_list:
                key, value = var.split(':')
                
                if key not in dii.keys():
                    dii[key] = ''
                    dii[key] = value
                    
                else:
                    dii[key] = value
            
            in_frame_indel_dict[gene][residue] = dii
            

# Saving the dictionary as txt file
with open(os.path.join(final_resdir,"inframe_indel_dict.txt"), 'w') as fp:
    json.dump(in_frame_indel_dict, fp)     

# --------------------------------------------------------------------------- #
# Cancer Genome Interpreter RESOURCE ---------------------------------------- # 

# Read in the Cancer Genome Interpreter tsv file in a dataframe 
# (The dataframe was downloaded as .txt format from here:
#  https://www.cancergenomeinterpreter.org/mutations . 
#  Accessed 01 February 2024.)

CGI_db = pd.read_csv(os.path.join(orig_resdir,"catalog_of_validated_oncogenic_mutations20240201.tsv"), 
                     sep='\t')

# Extract rows with more entities in the gdna column, separated by 
# the "__" character
CGI_db_red = CGI_db[CGI_db['gdna'].str.contains('__')]

# Remove CGI_db_red columns from CGI_db which will be processed afterwards
CGI_db_excluded = pd.merge(CGI_db, CGI_db_red, how='outer', indicator=True)
CGI_db_excluded = CGI_db_excluded.loc[CGI_db_excluded["_merge"] == "left_only"].drop("_merge", axis=1)

# Process rows in CGI_db_red, assign to a new df and concat this new df with 
# CGI_db_excluded, from which those rows were previously extracted.

# Create a dictionary for variants in CGI_db_red
cgi_db_red_dict = {}
count = 0

for index, row in CGI_db_red.iterrows():
    
    gdna = CGI_db_red.loc[index, "gdna"]
    gdna_list = gdna.split("__")
    
    # put in multiple lines the rows with more entities in the gdna column
    for element in gdna_list:
        
        count = 1 + count
        
        cgi_db_red_dict[count] = {'gene' : CGI_db_red.loc[index, "gene"],
                           'gdna' : element,
                           'protein': CGI_db_red.loc[index, "protein"],
                           'transcript': CGI_db_red.loc[index, "transcript"],
                           'info': CGI_db_red.loc[index, "info"],
                           'context': CGI_db_red.loc[index, "context"],
                           'cancer_acronym': CGI_db_red.loc[index, "cancer_acronym"],
                           'source': CGI_db_red.loc[index, "source"],
                           'reference': CGI_db_red.loc[index, "reference"]}

# Create dataframe from dictionary    
processed_CGI_db = pd.DataFrame.from_dict(cgi_db_red_dict,  orient='index',
                                          columns = ["gene", "gdna", "protein",
                                                     "transcript", "info", 
                                                     "context","cancer_acronym",
                                                     "source", "reference"])

# Concatenate CGI_db_excluded and processed_CGI_db
CGI_final = pd.concat([CGI_db_excluded, processed_CGI_db], ignore_index=True)

# Convert the dataframe of somatic variants to dictionary 
cgi_dict = {}

for index, row in CGI_final.iterrows():
    
    prot_change = CGI_final.loc[index, "protein"]
    source = CGI_final.loc[index, "source"]
    
    if prot_change == "." and source == "ClinVar":
        
        reference = CGI_final.loc[index, "reference"]
        
        if ">" in reference:
            
            if "__" in reference:
                
                reference = reference.split("__")[0]
                
                string, transcript, cdna = reference.split(":")
                cdna = cdna.split(" ")[0]
            
                gene = transcript[:-1]
                gene = gene.split("(")[1]
    
                cdna_plus_ref, alt_value = cdna.split(">")
                
                ref_value = cdna_plus_ref[-1]
                
                key = gene + ":" + cdna 
                
                if key not in cgi_dict:
                    
                    cgi_dict[key] = []
                    cgi_dict[key].append(alt_value)
                        
                else:
                    cgi_dict[key].append(alt_value)
                
            else:
                
                string, transcript, cdna = reference.split(":")
                cdna = cdna.split(" ")[0]
            
                gene = transcript[:-1]
                gene = gene.split("(")[1]
    
                cdna_plus_ref, alt_value = cdna.split(">")
                
                ref_value = cdna_plus_ref[-1]
                
                key = gene + ":" + cdna 
                
                if key not in cgi_dict:
  
                    cgi_dict[key] = []
                    cgi_dict[key].append(alt_value)
                        
                else:
                    cgi_dict[key].append(alt_value)
        
    else:
        
        key = CGI_final.loc[index, "gene"] + ":" + CGI_final.loc[index, "protein"]
        # Define reference and alternative alleles 
        ref_value = CGI_final.loc[index,'gdna'].strip()[-3]
        alt_value = CGI_final.loc[index,'gdna'].strip()[-1]
    
        if key not in cgi_dict:
            
            cgi_dict[key] = []
            cgi_dict[key].append(alt_value)
            
        else:
            cgi_dict[key].append(alt_value)
        
# Saving the dictionary as txt file
with open(os.path.join(final_resdir,"cgi_dictionary.txt"), 'w') as fp:
     json.dump(cgi_dict, fp)