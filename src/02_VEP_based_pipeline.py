# -*- coding: utf-8 -*-
"""
Created on Sat Jun 31 09:27:46 2025

@author: cartama

    Script to annotate variants with the Variant Effect Predictor (VEP) from
    
    Ensembl.

"""

## LIBRARIES

import argparse
from pathlib import Path

from colorama import init, Fore

# Automatically reset the style to normal after each print
init(autoreset=True)

import pandas as pd
import numpy as np
import os
import json
import subprocess

## ------------------------------------------------------------------------ ##

## FUNCTIONS

pd.DataFrame.iteritems = pd.DataFrame.items

# The function replaces in an empty df the data from the VEP output 
def dataframe_fill (records_list):
    
    dataframe = pd.DataFrame(records_list, columns=res_row_split[:7])
    
    return dataframe


# The function calls VEP
def call_vep_shell_script (input_file_dir, origin_path, output_file, vep_dir):
    
    session = subprocess.call(['bash', 'vep.sh', input_file_dir, origin_path, output_file, vep_dir])
    
    return session

# The function checks the presence of the string "CSQ" 
def slicer(my_str,sub):
    
   index=my_str.find(sub)
   
   if index !=-1 :
       
         return my_str[index:] 
     
   else :
       
         raise Exception('Sub string not found!')
         
def make_hyperlink(url, value):

    return '=HYPERLINK("%s", "%s")' % (url.format(value), value)

#---------------------------------------------------------------------------##
#     MAIN
#---------------------------------------------------------------------------##

if __name__ == '__main__':
    
    
    # ArgumentParser object to analyse strings given from command-line 
    # by the user 
    parser = argparse.ArgumentParser(prog='VEP-based Variant Annotation Pipeline', 
                                     description="Annotates input file with VEP by Ensembl")
    
    # Input txt/VCF file
    parser.add_argument(
        '-i','--input',
        type = str,
        required = True,
        help = 'Full path to input txt/VCF file'
        )   
    
    global args
    args = parser.parse_args()
    
## ------------------------------------------------------------------------ ##

### EXAMPLE DATA DIRECTORY

print("")
print(f"{Fore.MAGENTA}=" * 100)
print(f"{Fore.MAGENTA}Starting VEP-based Variant Annotation Pipeline")
print(f"{Fore.MAGENTA}=" * 100)
print("")

# Input data with absolute path
input_data_dir = Path(args.input).resolve()

# Parent dir of input data
full_path_input = os.path.dirname(Path(args.input).resolve())
print(f"Your test data is located here: {full_path_input}")
print("")

## ------------------------------------------------------------------------ ##

### VEP DIRECTORY

# Re-build path of .vep from input data full path  
vep_dir = os.path.join(os.path.dirname(full_path_input), ".vep")
                       
# Check on the presence of VEP API
if os.path.exists(vep_dir):
    
    print(f"Your VEP is located here: {vep_dir}")
    print("")

else:
    print(f"Your VEP is not located here as expected: {vep_dir}")
    print("")
    
#---------------------------------------------------------------------------##

### CLINVAR RESOURCE

print(f"{Fore.CYAN}Upload ClinVar resource")
print(f"{Fore.CYAN}-" * 100)
print("")

# Re-build path of resources
resources_data_dir = os.path.join(os.path.dirname(full_path_input), "resources")

# Open the clinvar database as dictionary
with open(os.path.join(resources_data_dir, "clinvar_all_dictionary.json"), 'r') as fp:
    clinvar_all_dict = json.load(fp)

#---------------------------------------------------------------------------##

### PREPARE DATA FOR ANNOTATION WITH VEP 

# origin_path: parent dir of the input data file 
# input_data: name of the input data file 
origin_path = full_path_input

input_data = os.path.split(input_data_dir)[1]

print(f"Your test data: {input_data}")
print("")

# Input data name 
input_data_name = input_data.split(".txt")[0]

print(f"Name of your test data: {input_data_name}")
print("")

# Path to VEP output 
output_file = os.path.join(input_data_name + "_annotated.vcf")

print(f"Name of your annotated output: {output_file}")
print("")

#---------------------------------------------------------------------------##

### RUN VEP 

print(f"{Fore.MAGENTA}=" * 100)
print(f"{Fore.MAGENTA}Creating VEP output")
print(f"{Fore.MAGENTA}=" * 100)
print("")

# Call the VEP function
vep_output = call_vep_shell_script(input_data_dir, origin_path, output_file, vep_dir)

#---------------------------------------------------------------------------##

### POST-PROCESS VEP OUTPUT

print("")
print(f"{Fore.MAGENTA}=" * 100)
print(f"{Fore.MAGENTA}Post-processing VEP output")
print(f"{Fore.MAGENTA}=" * 100)
print("")

# Annotated VEP output
annotated_output = os.path.join(origin_path, "vep_" + output_file)
    
# Read in the VEP output in a list and remove the header
output = []
with open(annotated_output, 'r') as f:
    lines = f.readlines()
    criteria_1 = '##'
    output = [line for line in lines if criteria_1 not in line]
    
    # Isolate the description of the INFO field 
    criteria_2 = '##INFO=<ID=CSQ'
    info_header = next(line for line in lines if criteria_2 in line)
    
# Parse the output file for data preprocessing
records_list = [] 
dict_of_predictions = {}

for index, row in enumerate(output):

    if row.startswith('#'):

        res_row = row.replace('#','')

        res_row_split = res_row.split()
        
    else:
        
        row_split_data = row.split()
        
        # Store the annotated consequences in a dictionary
        dict_of_predictions[index] = row_split_data[7]
        records_list.append(row_split_data[:7])


global annotated_df
    
# Create dataframe from data with the function dataframe_fill
annotated_df = dataframe_fill(records_list)

annotated_df.index = np.arange(1, len(annotated_df)+1)           
    
#---------------------------------------------------------------------------##

### HANDLING INFO FIELD

# Create new dataframe only for the INFO field     
df_of_info = pd.DataFrame({'INFO':dict_of_predictions})

# Predictions ordered in a dataframe (predictions with the same index are referred to the same variant) 
df_info_ordered = (df_of_info["INFO"].str.split(",", expand=True).stack().to_frame("words").reset_index(1, drop=True))

df_info_ordered["count"] = df_info_ordered.groupby(level=0).cumcount() // 1

# Predictions referred to the same variant are now grouped in the same row 
df_info_ordered = df_info_ordered.rename_axis("idx").groupby(["idx", "count"])["words"].agg(" ".join).unstack(1)

# Now, handle the INFO header from the "Allele" string 
header = slicer(info_header, "Allele")

# Remove last three elements form header
header = header[:-3]
header_to_list = header.split('|')

#---------------------------------------------------------------------------##

### HANDLING CANONICAL TRANSCRIPTS

# Add a flag that indicates the type of transcript in the first 
# position in the VEP output. Allowed values are: Canonical and VEP first
annotated_df["Transcript #1"] = None

# Put the prediction associated to the canonical transcript in the first position,
# if present.      
for index_row, row in df_info_ordered.iterrows():
    
    row_true = row.str.contains('YES')
    
    if row_true.any():
                
        annotated_df.loc[index_row, "Transcript #1"] = "Canonical"
        
        # Count how many predictions wew annotated for the specific variant 
        len_row = list(range(len(row)))
        
        # Get the prediction/s for which CANONICAL is YES
        index_to_shift = row_true[row_true == True].index
        
        # Check if more than 1 prediction associated with canonical transcripts 
        len_trues = len(index_to_shift)
        
        # If more than 1 prediction
        if len_trues > 1:

            for true_idx in index_to_shift:
                
                # Split the predictions by pipe and create a data frame
                row_splitted_list = row[true_idx].split("|")
                row_splitted_df = pd.DataFrame([row_splitted_list]) 
                row_splitted_df_wheader = row_splitted_df.rename(columns=lambda x: header_to_list[x])

                # If the "HGVSc" field is empty, then refuse the prediction and
                # accept the following
                if row_splitted_df_wheader.loc[0,"HGVSc"] == '':
                    
                    # Pass to the next iteration 
                    continue
                
                else:

                    # Move the prediction associated with the canonical transcript 
                    # and without empty HGVSc in the first position
                    new_idx = [i for i in len_row if i!=true_idx]
                    
                    row_modified = row.loc[[true_idx]+new_idx]
        
                    row_modified_list = row_modified.to_list()
        
                    df_info_ordered.loc[index_row, :] = row_modified_list
                    
                    break
                    
        else:
            
            # Move the prediction associated with the canonical transcript 
            # and without empty HGVSc in the first position
            new_idx = [i for i in len_row if i!=index_to_shift[0]]
            
            row_modified = row.loc[[index_to_shift[0]]+new_idx]

            row_modified_list = row_modified.to_list()

            df_info_ordered.loc[index_row, :] = row_modified_list
            
    else:
        
        annotated_df.loc[index_row,"Transcript #1"] = "VEP first"
        continue 

#---------------------------------------------------------------------------##

### COUNT TRANSCRIPTS PER VARIANT

# Count and store in the dataframe the number of consequences for each variant
annotated_df['count_cons'] = (~df_info_ordered.isna()).sum(1)
    
# Columns of interest in the header 
info_header_red = ["Consequence","SYMBOL","MANE_SELECT","EXON","INTRON",
                   "HGVSc","HGVSp","HGVSg","Amino_acids","Existing_variation",
                   "CANONICAL","AF","gnomADe_AF","gnomADe_AFR_AF",
                   "gnomADe_EAS_AF","gnomADe_NFE_AF","gnomADe_AMR_AF",
                   "gnomADe_SAS_AF","gnomADg_AF","gnomADg_AFR_AF",
                   "gnomADg_EAS_AF","gnomADg_NFE_AF","gnomADg_AMR_AF",
                   "gnomADg_SAS_AF","PUBMED",
                   "phyloP100way_vertebrate_rankscore",
                   "phastCons100way_vertebrate_rankscore",
                   "SpliceAI_cutoff"]


# Dataframe with a variant in each line 
for column_name, item in df_info_ordered.iteritems():
    
    dataframe_n = df_info_ordered.iloc[:,column_name]
    df_info_splitted = dataframe_n.str.rsplit("|", expand=True)
    df_info_splitted = df_info_splitted.rename(columns=lambda x: header_to_list[x])
    
    # Keep only the reduced header 
    df_info_filtered = df_info_splitted[info_header_red]
    annotated_df = pd.concat([annotated_df, df_info_filtered], axis=1)

# Replace Nan values with empty string
final_annotated_df = annotated_df.fillna('')    

# Path to the final output 
file_output_path = os.path.join(origin_path,str(input_data_name) + "_vep.xlsx")

print(f"{Fore.CYAN}Saving the final output:")
print(f"{Fore.CYAN}{file_output_path}")
print(f"{Fore.CYAN}-" * 100)
print("")

# Save the final output to file
final_annotated_df.to_excel(file_output_path, index=False)
  
#---------------------------------------------------------------------------##

### VARIANT CLASSIFICATION PIPELINE

print(f"{Fore.MAGENTA}=" * 100)
print(f"{Fore.MAGENTA}Variant Classification Pipeline")
print(f"{Fore.MAGENTA}=" * 100)
print("")

## Input to the variant classification pipeline
vep_output = pd.read_excel(file_output_path, na_values = 'NA', keep_default_na=False)

# Subset the final output to the columns of interest 
vep_output_subset = vep_output[['SYMBOL','HGVSc','HGVSp','HGVSg','MANE_SELECT',
                                'EXON','CHROM','POS','REF', 'ALT','Amino_acids',
                                'Existing_variation','Consequence',
                                'Transcript #1', 'gnomADe_AF','gnomADe_AFR_AF',
                                'gnomADe_EAS_AF','gnomADe_NFE_AF',
                                'gnomADe_AMR_AF','gnomADe_SAS_AF','gnomADg_AF',
                                'gnomADg_AFR_AF','gnomADg_EAS_AF',
                                'gnomADg_NFE_AF','gnomADg_AMR_AF',
                                'gnomADg_SAS_AF',
                                'phyloP100way_vertebrate_rankscore',
                                'phastCons100way_vertebrate_rankscore',
                                'SpliceAI_cutoff']]                 

#---------------------------------------------------------------------------##

### CLINVAR DATABASE QUERY

print(f"{Fore.CYAN}ClinVar database query")
print(f"{Fore.CYAN}-" * 100)
print("")

# List of ClinVar data 
list_all_clinvar = []

list_indexes_all_clinvar = []

# Iterate over the input dataframe
for index, row in vep_output_subset.iterrows():
    
    try:
        list_all_clinvar.append(clinvar_all_dict[str(vep_output_subset.loc[index, "CHROM"])][str(vep_output_subset.loc[index, "POS"])][str(vep_output_subset.loc[index, "REF"])][str(vep_output_subset.loc[index, "ALT"])])
        
        list_indexes_all_clinvar.append(index)
        
    except:
        pass
    
#  Convert the data collected from ClinVar to dataframe
columns = ["ClinVar_VariationID",
           "ClinVar_germline", 
           "ClinVar_germline_ReviewStatus",
           "ClinVar_Oncogenicity",
           "ClinVar_Oncogenicity_ReviewStatus",
           "ClinVar_Somatic_Clinical_Impact",
           "ClinVar_Somatic_ReviewStatus"]

clinvar_df = pd.DataFrame(list_all_clinvar, columns = columns)

# Change indexes
new_labels = np.array(list_indexes_all_clinvar)
clinvar_df = clinvar_df.set_index(new_labels)

# Merge the ClinVar dataframe to the original one
vep_output_subset = pd.concat([vep_output_subset, clinvar_df], axis=1)

# base clinvar url 
clinvar_url = "http://www.ncbi.nlm.nih.gov/clinvar/?term={}"

# Fill the Clinvar_url column of the dataframe with the URLs masked by the ClinVar_VariationID column names
vep_output_subset['Clinvar_url'] = vep_output_subset['ClinVar_VariationID'].apply(lambda x: make_hyperlink(clinvar_url, x))

# Remove the temporary column that contained the names that replaced the URLs 
vep_output_subset = vep_output_subset.drop('ClinVar_VariationID', axis=1)

#---------------------------------------------------------------------------##

### SAVING THE CLASSIFICATION OUTPUT

# Save the final output
file_predicted_output_path = os.path.join(origin_path, str(input_data_name) + "_vep_classification.xlsx")        
vep_output_subset.to_excel(file_predicted_output_path, index = False)

print(f"{Fore.CYAN}Saving the classification output:")
print(f"{Fore.CYAN}{file_predicted_output_path}")
print(f"{Fore.CYAN}-" * 100)
print("")

#---------------------------------------------------------------------------##

### RUN OncoVI

print("")
print(f"{Fore.MAGENTA}Running OncoVI on the VEP annotated output")
print(f"{Fore.MAGENTA}=" * 100)
print("")

# Resources for OncoVI
resources_oncovi = resources_data_dir

# Path to python script for OncoVI evaluation of the oncogenicity guidelines
main = os.path.join(os.getcwd(), "03_OncoVI_SOP.py")

# Command to run OncoVI script
cmd = " python "+main+" -i "+file_predicted_output_path+" -r "+resources_oncovi+" -d "+input_data_name
os.system(cmd)

#---------------------------------------------------------------------------##

### ANALYSIS COMPLETED

print(f"{Fore.MAGENTA}=" * 100)
print(f"{Fore.MAGENTA}Analysis completed")
print(f"{Fore.MAGENTA}=" * 100)
print("")