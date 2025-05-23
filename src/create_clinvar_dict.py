# -*- coding: utf-8 -*-
"""
Created on Wed May 21 12:22:38 2025

    "The script creates the dictionary clinvar_all_dictionary.txt
    
    from the variant_summary_GRCh38_red.txt file."

@author: cartama

"""

# Import libraries 
import argparse
import pandas as pd
import os
import json
from collections import defaultdict

if __name__ == '__main__':
    
    # ArgumentParser object to analyse multiple strings from command-line 
    parser = argparse.ArgumentParser(prog='Creation of clinvar_all_dictionary.txt dictionary', 
                                     description="Creates the clinvar_all_dictionary.txt dictionary starting from variant_summary_GRCh38_red.txt file")
    
    # Parser arguments    
    
    # Input file variant_summary_GRCh38_red.txt file
    parser.add_argument(
        '-i','--input',
        type = str,
        required = True,
        help = 'Path to the folder containing the variant_summary_GRCh38_red.txt file.'
        )   
    
    # Output folder
    parser.add_argument(
        '-r','--resources',
        type = str,
        required = True,
        help = 'Path to the folder where to save the clinvar_all_dictionary.txt dictionary.'
        )
    
    global args
    args = parser.parse_args()

# all unique variants present in ClinVar 
clinvar_df_dir = args.input

# Create directory for final dictionary
final_datadir = args.resources

# If the final directory doesn´t exist create it
os.makedirs(final_datadir, exist_ok=True)

######## Read in ClinVar DATABASE ###############
print("######## Read in ClinVar DATABASE ###############")

# Upload the ClinVar dictionary
with open(clinvar_df_dir, "r") as f:
    lines = f.readlines()

# Parsing in DataFrame
clinvar_df = pd.DataFrame({'DATA': lines})
clinvar_df2 = clinvar_df['DATA'].str.split("|", expand=True)

header = ["GeneSymbol",
          "ClinicalSignificance", 
          "Chromosome", 
          "ReviewStatus", 
          "NumberSubmitters",
          "VariationID", 
          "PositionVCF", 
          "ReferenceAlleleVCF", 
          "AlternateAlleleVCF",
          "SomaticClinicalImpact", 
          "ReviewStatusClinicalImpact", 
          "Oncogenicity", 
          "ReviewStatusOncogenicity"]

clinvar_df2.columns = header

# Dinamic and nested structure with defaultdict
def nested_dict():
    return defaultdict(nested_dict)

clinvar_dict = nested_dict()

# Optimized iteration
for row in clinvar_df2.itertuples(index=False):
    
    var_id = row.VariationID
    chromosome = row.Chromosome
    pos = row.PositionVCF
    refAA = row.ReferenceAlleleVCF
    altAA = row.AlternateAlleleVCF

    clinical_significance = row.ClinicalSignificance
    review_status = row.ReviewStatus
    oncogenicity = row.Oncogenicity
    onco_review_status = row.ReviewStatusOncogenicity
    somatic_clinical_impact = row.SomaticClinicalImpact
    somatic_review_status = row.ReviewStatusClinicalImpact

    clinvar_data = [
        var_id,
        clinical_significance,
        review_status,
        oncogenicity,
        onco_review_status,
        somatic_clinical_impact,
        somatic_review_status
    ]

    clinvar_dict[chromosome][pos][refAA][altAA] = clinvar_data

######## Save the final ClinVar output ##########
print("######## Save the final ClinVar output ##########")

# Final conversion to dictionary
clinvar_dict = json.loads(json.dumps(clinvar_dict))

# Save the dictionary
with open(os.path.join(final_datadir, "clinvar_all_dictionary.json"), "w") as f_out:
    json.dump(clinvar_dict, f_out)
