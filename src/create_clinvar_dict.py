# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:22:38 2024

@author: cartama

"""

import pandas as pd
import os
import json

clinvar_df_dir = "/path/to/variant_summary_red.txt"

######## ClinVar DATABASE ##########

# Open the file as text file
file_clinvar = open(clinvar_df_dir, "r")

# Each line is saved in a list
lines = file_clinvar.readlines()

# Remove the first line, i.e., the header of the file
lines.pop(0)

# Create the dictionary to save the lines of the file 
clinvar_dict = {}

clinvar_df = pd.DataFrame({'DATA':lines})
clinvar_df_split = clinvar_df['DATA'].str.split("|", expand=True)

header = ["GeneSymbol",
          "ClinicalSignificance",
          "Chromosome",
          "Start",
          "VariationID",
          "ReferenceAlleleVCF",
          "AlternateAlleleVCF",
          "ReviewStatus",
          "NumberSubmitters"]

# add the header
clinvar_df_w_header = clinvar_df_split.rename(columns=lambda x: header[x])

# create the dictionary 
for index, row in clinvar_df_w_header.iterrows():
    
    chromosome = clinvar_df_w_header.loc[index, "Chromosome"]
    start = clinvar_df_w_header.loc[index, "Start"]
    refAA = clinvar_df_w_header.loc[index, "ReferenceAlleleVCF"]
    altAA = clinvar_df_w_header.loc[index, "AlternateAlleleVCF"]
    
    clinvar_data = [clinvar_df_w_header.loc[index, "VariationID"],
                    clinvar_df_w_header.loc[index, "ClinicalSignificance"],
                    clinvar_df_w_header.loc[index, "ReviewStatus"],
                    clinvar_df_w_header.loc[index, "NumberSubmitters"]]
    
    
    if chromosome not in clinvar_dict.keys():
        
        clinvar_dict[chromosome] = {}
        
        if start not in clinvar_dict[chromosome].keys():
            
            clinvar_dict[chromosome][start] = {}
            
            if refAA not in clinvar_dict[chromosome][start].keys():
                
                clinvar_dict[chromosome][start][refAA] = {}
                
                if altAA not in clinvar_dict[chromosome][start][refAA].keys():
                    
                    clinvar_dict[chromosome][start][refAA][altAA] = ""
                    clinvar_dict[chromosome][start][refAA][altAA] = clinvar_data
            
            else:
                if altAA not in clinvar_dict[chromosome][start][refAA].keys():
                    
                    clinvar_dict[chromosome][start][refAA][altAA] = ""
                    clinvar_dict[chromosome][start][refAA][altAA] = clinvar_data
    
        else:
            if refAA not in clinvar_dict[chromosome][start].keys():
                
                clinvar_dict[chromosome][start][refAA] = {}
                
                if altAA not in clinvar_dict[chromosome][start][refAA].keys():
                    
                    clinvar_dict[chromosome][start][refAA][altAA] = ""
                    clinvar_dict[chromosome][start][refAA][altAA] = clinvar_data
            
            else:
                if altAA not in clinvar_dict[chromosome][start][refAA].keys():
                    
                    clinvar_dict[chromosome][start][refAA][altAA] = ""
                    clinvar_dict[chromosome][start][refAA][altAA] = clinvar_data
    else:
        
        if start not in clinvar_dict[chromosome].keys():
            
            clinvar_dict[chromosome][start] = {}
            
            if refAA not in clinvar_dict[chromosome][start].keys():
                
                clinvar_dict[chromosome][start][refAA] = {}
                
                if altAA not in clinvar_dict[chromosome][start][refAA].keys():
                    
                    clinvar_dict[chromosome][start][refAA][altAA] = ""
                    clinvar_dict[chromosome][start][refAA][altAA] = clinvar_data
            
            else:
                
                if altAA not in clinvar_dict[chromosome][start][refAA].keys():
                    
                    clinvar_dict[chromosome][start][refAA][altAA] = ""
                    clinvar_dict[chromosome][start][refAA][altAA] = clinvar_data
        else:
            
            if refAA not in clinvar_dict[chromosome][start].keys():
                
                clinvar_dict[chromosome][start][refAA] = {}
                
                if altAA not in clinvar_dict[chromosome][start][refAA].keys():
                    
                    clinvar_dict[chromosome][start][refAA][altAA] = ""
                    clinvar_dict[chromosome][start][refAA][altAA] = clinvar_data
                
            else:
                
                if altAA not in clinvar_dict[chromosome][start][refAA].keys():
                    
                    clinvar_dict[chromosome][start][refAA][altAA] = ""
                    clinvar_dict[chromosome][start][refAA][altAA] = clinvar_data
                    
output_dir = "/path/to/output"

# Save the dictionary to file 
with open(os.path.join(output_dir, "clinvar_all_dictionary.txt"), 'w') as fp:
    json.dump(clinvar_dict, fp)