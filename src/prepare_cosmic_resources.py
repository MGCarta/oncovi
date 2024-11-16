# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 10:44:29 2024

@author: cartama

    The script preparares COSMIC resources used in OncoVI. 
    
"""

# --------------------------------------------------------------------------- #

# Libraries
import pandas as pd
import os  
import json

# --------------------------------------------------------------------------- #

# Base directory
base_dir = "/path/to/working/directory"

# Directory for original resources (i.e., resources not modified)
orig_resdir = os.path.join(base_dir, "01_original_resources")

# Directory for final resources 
final_resdir = os.path.join(base_dir, "02_resources")


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

# Save the file to directory
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
    
# Variants with a reference AA
refAA_valid = cosmic_all_data[cosmic_all_data['AA_REF'] != ''] 
# all variants have a ref AA

# Add a new column to the refAA_valid dataframe
refAA_valid = refAA_valid.assign(MUT_AA_mod = None)

# Change nomenclature of synonymous variants
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

# Save the file to directory
with open(os.path.join(final_resdir,"cosmic_all_dictionary.txt"), 'w') as fp:
    json.dump(cosmic_all_dict, fp)
    