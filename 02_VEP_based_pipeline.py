# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 09:27:46 2024

@author: cartama

    This script annotates the variants of the SOP data set relying on the 
    
    Variant Effect Predictor (VEP) of Ensembl and collects the classification
    
    of pathogenicity via the ClinVar API. The output is a tabular data 
    
    containing the inmformation required by OncoVI for the evaluation of the 
    
    oncogenicity guidelines.

"""

## ------------------------------------------------------------------------ ##

# Directory in which VEP is located 
vep_dir = "/research/storage-normal/vep111/vep/.vep"

## ------------------------------------------------------------------------ ##

# Libraries
import pandas as pd
import numpy as np
import argparse
import os
import json
import subprocess
import xml.etree.cElementTree as ET

pd.DataFrame.iteritems = pd.DataFrame.items

# This function replaces in an empty df the data of the VEP output 
def dataframe_fill (records_list):
    
    dataframe = pd.DataFrame(records_list, columns=res_row_split[:7])
    
    return dataframe


# This function calls the VEP
def call_vep_shell_script (file_dir, file, dest_path, vep_dir):
    
    session = subprocess.call(['bash', 'vep.sh', file_dir, file, dest_path, vep_dir])
    
    return session

# This function checks the presence of the string "CSQ" 
def slicer(my_str,sub):
    
   index=my_str.find(sub)
   
   if index !=-1 :
       
         return my_str[index:] 
     
   else :
       
         raise Exception('Sub string not found!')
         

# This function makes an URL clickable 
def make_clickable(val):
    
    return '<a href="{}">{}</a>'.format(val, val)

# This function created an hyperlink
def make_hyperlink(url, value):

    return '=HYPERLINK("%s", "%s")' % (url.format(value), value)


#---------------------------------------------------------------------------##
#     MAIN
#---------------------------------------------------------------------------##

if __name__ == '__main__':
    
    
    # ArgumentParser object to analyse strings given from command-line 
    # by the user 
    parser = argparse.ArgumentParser(prog='In-house VEP-based pipeline',
                                     description="Annotates variants with the Variant Effect Predictor (VEP) by Ensembl v. 111.0")
    
    # Input txt file
    parser.add_argument(
        '-i','--input',
        type = str,
        required = True,
        help = 'Full path to input txt file'
        )   
    
    global args
    args = parser.parse_args()

# Check on the txt file given as input
input_data_dir = args.input

if vep_dir:
    
    pass

else:
    
    raise TypeError("Path to vep/.vep.")

removal_list_of_consequences = "/care/storage-normal/vep/.vep/resources/resources_2023_01/cons_to_remove_2023_01.txt"

# Consequences to filter out
removal_list_of_cons = []

with open (removal_list_of_consequences, 'r') as f:
    
    lines = f.readlines()
    removal_list_of_cons = [line.strip() for line in lines if line not in removal_list_of_cons]

    
#---------------------------------------------------------------------------##

print("\n")
print('******************    VEP-based Variant Annotation Pipeline    ******************')
print('===============================================================================\n')
    
#---------------------------------------------------------------------------##
        
print('*****************    Uploading ClinVar database    ******************')
print('=====================================================================\n')

######## ClinVar DATABASE ##########

clinvar_data_dir = os.path.join(vep_dir, "ClinVar")

# Open the clinvar database as dictionary
with open(os.path.join(clinvar_data_dir, "clinvar_all_dictionary.txt"), 'r') as fp:
    clinvar_all_dict = json.load(fp)

#---------------------------------------------------------------------------##

print("")
print(args)
print("\n")
print("Input txt file:", args.input)
print("\n")    

# origin_path: parent dir of the input data file 
# input_data: name of the input data file 
origin_path, input_data = os.path.split(input_data_dir)

# Input data name 
input_data_name = input_data.split(".txt")[0]

# Path to VEP output 
dest_path_output = os.path.join(origin_path, input_data_name + "_annotated.vcf")


print('**************************    VEP Output Computing    ***************************')
print('===============================================================================\n')
print("")

# Call the VEP function
vep_output = call_vep_shell_script(input_data_dir, input_data, dest_path_output, vep_dir)
    
print("\n")
print('************************    Post-Processing VEP Output    ***********************')
print('===============================================================================\n')
print("")

# Annotated VEP output
annotated_output = dest_path_output
    
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
    
# Create dataframe from data calling the function dataframe_fill
annotated_df = dataframe_fill(records_list)

annotated_df.index = np.arange(1, len(annotated_df)+1)           
    
########## INFO FIELD  ########################################################

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

########## CANONICAL TRANSCRIPTS ##############################################

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

########## COUNT TRANSCRIPTS PER VARIANT ######################################   

# Count and store in the dataframe the number of consequences for each variant
annotated_df['count_cons'] = (~df_info_ordered.isna()).sum(1)
    
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
    
    # Keep only the header of interest 
    df_info_filtered = df_info_splitted[info_header_red]
    annotated_df = pd.concat([annotated_df, df_info_filtered], axis=1)

# Replace Nan values with empty string
final_annotated_df = annotated_df.fillna('')    

# Path to the final output 
file_output_path = os.path.join(origin_path,str(input_data_name)+"_vep.xlsx")

print("Saving the final output:", file_output_path)
print("")

# Save the final output to file
final_annotated_df.to_excel(file_output_path, index=False)
  

print('********************    Variant Classification Pipeline    ********************')
print('=============================================================================\n')
print("\n")  

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

print('********************        ClinVar database query         ********************')
print('=============================================================================\n')
print("\n")  

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
           "NumberSubmitters"]

clinvar_df = pd.DataFrame(list_all_clinvar, columns = columns)

# Change indexes
new_labels = np.array(list_indexes_all_clinvar)
clinvar_df = clinvar_df.set_index(new_labels)

# Merge the dataframe of clinvar data to the original one
vep_output_subset = pd.concat([vep_output_subset, clinvar_df], axis=1)

# base clinvar url 
clinvar_url = "http://www.ncbi.nlm.nih.gov/clinvar/?term={}"

# Fill the Clinvar_url column of the dataframe with the URLs masked by the ClinVar_VariationID column names
vep_output_subset['Clinvar_url'] = vep_output_subset['ClinVar_VariationID'].apply(lambda x: make_hyperlink(clinvar_url, x))

# Remove the temporary column that contained the names that replaced the URLs 
vep_output_subset = vep_output_subset.drop('ClinVar_VariationID', axis=1)

                                                            
# dbSNP URLs to display data 
# *****************************************************************************
#
# The script investigates the "Existing_variation" from the VEP output.
#
# If a "Reference SNP ID" (i.e., "rs" ID) is available, then it is provided.
#       
# *****************************************************************************

# New column for the "dbSNP identifier"
vep_output_subset = vep_output_subset.assign(dbSNP_key = None)


print('********************          dbSNP database query         ********************')
print('=============================================================================\n')
print("\n")

# Create a df containing the "Reference SNP ID"
Existing_variation_dict = vep_output.Existing_variation.to_dict()

# Base dbSNP url 
dbSNP_url = "http://www.ncbi.nlm.nih.gov/snp/?term={}"

for index, row in vep_output_subset.iterrows():
    
    vep_output_subset.loc[index, 'dbSNP_key'] = next((i for i in Existing_variation_dict[index].split("&") if i.startswith('rs')), None)
            
# URLs masked by the dbSNP_key column names
vep_output_subset['dbSNP_url'] = vep_output_subset['dbSNP_key'].apply(lambda x: make_hyperlink(dbSNP_url, x) if x != None else x)


print('********************     Saving the classification output      ****************')
print('=============================================================================\n')
print("\n")

file_predicted_output_path = os.path.join(origin_path, str(input_data_name)+"_prediction.xlsx")        
vep_output_subset.to_excel(file_predicted_output_path, index = False)

print('****************************    Analysis Completed    ***************************')
print('===============================================================================\n')