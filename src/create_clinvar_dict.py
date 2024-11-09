# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:22:38 2024

@author: cartama
"""

import pandas as pd
import os
import json

clinvar_df_dir = "W:/pathologie/bioinfo-archive/VEP/ClinVar/variant_summary_red.txt"

######## ClinVar DATABASE ##########

# Open the file as text file
file_clinvar = open(clinvar_df_dir, "r")

# Each line is saved in a list
lines = file_clinvar.readlines()

# Remove the first line, which is the header of the file
lines.pop(0)

# Create the dictionary in which I want to save the lines of the file 
clinvar_dict = {}

clinvar_df = pd.DataFrame({'DATA':lines})
clinvar_df2 = clinvar_df['DATA'].str.split("|", expand=True)

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
clinvar_df2_w_header = clinvar_df2.rename(columns=lambda x: header[x])

# create the dictionary 
for index, row in clinvar_df2_w_header.iterrows():
    
    chromosome = clinvar_df2_w_header.loc[index, "Chromosome"]
    start = clinvar_df2_w_header.loc[index, "Start"]
    refAA = clinvar_df2_w_header.loc[index, "ReferenceAlleleVCF"]
    altAA = clinvar_df2_w_header.loc[index, "AlternateAlleleVCF"]
    
    clinvar_data = [clinvar_df2_w_header.loc[index, "VariationID"],
                    clinvar_df2_w_header.loc[index, "ClinicalSignificance"],
                    clinvar_df2_w_header.loc[index, "ReviewStatus"],
                    clinvar_df2_w_header.loc[index, "NumberSubmitters"]]
    
    
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
                    
datadir = "V:/gruppen/AG_Bioinfo/members/Carta/TSO500_HRD/data"

# Save the dictionary to file 
with open(os.path.join(datadir, "clinvar_all_dictionary.txt"), 'w') as fp:
    json.dump(clinvar_dict, fp)