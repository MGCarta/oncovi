# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:27:14 2025


    The script applies the oncogenicity guidelines implemented in OncoVI to the
    
    93 somatic variants of the Standard Operating Procedure (SOP) data set.    

"""

## ------------------------------------------------------------------------ ##

## Libraries
import json
import re
import os
import pandas as pd

## ------------------------------------------------------------------------ ##

## Functions

# Computes grantham score
def compute_grantham_score(ref_amino_acid, alt_amino_acid, original_alt_amino):
    
    # ref_amino_acid is the reference AA of the variant evaluated and 
    # of the variant we know to be oncogenic (because in CGI)
    
    # alt_amino_acid is the alternate AA of the variant evaluated
    
    # original_alt_amino is the alternate AA of the variant we know to
    # be oncogenic 
    
    # grantham score for the variant known to be oncogenic (because in CGI)
    grantham_score_orig = grantham_table.loc[ref_amino_acid, original_alt_amino]
    
    # grantham score for the variant evaluated 
    grantham_score_new = grantham_table.loc[ref_amino_acid, alt_amino_acid]
    return grantham_score_orig, grantham_score_new
        
## Input data 

input_data_dir = "/path/to/oncovi/testdata"

# -------------------------------------------------------------------------- ##

## Path for the output 
output_dir = os.path.join(input_data_dir, "OncoVI_eval")

if os.path.isdir(output_dir) == True:
    pass
else:
    os.mkdir(output_dir)

# -------------------------------------------------------------------------- ##

## Resources 
resources_dir = "/path/to/oncovi/resources"

# Tumor suppressor genes (TSGs) in Cancer Gene Census Tier 1 and OncoKB
tsg_bona_fide = pd.read_csv(os.path.join(resources_dir, "bona_fide_tsg.txt"), 
                            header = None, 
                            sep='\t',
                            names = ["Gene"])

tsg_bona_fide_list = tsg_bona_fide.Gene.tolist()

# TSGs in Cancer Gene Census and OncoKB
tsg = pd.read_csv(os.path.join(resources_dir, "tsg_list.txt"),
                  header = None, 
                  sep='\t',
                  names = ["Gene"])

tsg_list = tsg.Gene.tolist()

# Oncogenes list in Cancer Gene Census and OncoKB
oncogenes = pd.read_csv(os.path.join(resources_dir, "ogs_list.txt"), 
                        header = None, 
                        sep='\t',
                        names = ["Gene"])

oncogenes_list = oncogenes.Gene.tolist()

# Oncogenic variants [from Cancer Genome Interpreter (CGI)]
with open(os.path.join(resources_dir, "cgi_dictionary.txt"), 'r') as fp:
    onco_var_cgi = json.load(fp)

# Dictionary for MutSpliceDB 
with open(os.path.join(resources_dir, "mut_splice_dict.txt"), 'r') as fp:
    mut_splice_dict = json.load(fp)

# cancerhotspots.org database
# single residues from cancerhotspots.org database
with open(os.path.join(resources_dir, "single_residue_dict.txt"), 'r') as fp:
    single_residue_dict = json.load(fp)

# in-frame indel from cancerhotspots.org database
with open(os.path.join(resources_dir, "inframe_indel_dict.txt"), 'r') as fp:
    in_frame_indel_dict = json.load(fp)

# Grantham score
grantham_table = pd.read_csv(os.path.join(
    resources_dir, "grantham.tsv"), index_col=0, sep='\t')

# Dictionary for AA conversion from three letter convention to one 
# letter convention.
with open(os.path.join(resources_dir, "amino_dict.txt"), 'r') as fp:
    amino_dict = json.load(fp)
    
# Dictionary for UniProt domains annotation
with open(os.path.join(resources_dir, "domains_dictionary.txt"), 'r') as fp:
    domains_dict = json.load(fp)
    
# ClinVar ClinicalSignificance records manually selected for OS2 criterion
ClinicalSignificance = pd.read_csv(os.path.join(resources_dir, "os2_manually_selected.txt"), 
                                   header = None, 
                                   sep='\t',
                                   names = ["value"])

ClinicalSignificance_list = ClinicalSignificance.value.tolist()
    
# COSMIC dictionary for entries without HGVSp but with HGVSG identifier
with open(os.path.join(resources_dir, "cosmic_hgvsg_dictionary.txt"), 'r') as fp:
    cosmic_hgvsg_dict = json.load(fp)
    
# COSMIC dictionary for entries with genomic position GRCh38 identifier
with open(os.path.join(resources_dir, "cosmic_all_dictionary.txt"), 'r') as fp:
    cosmic_all_dict = json.load(fp)

## ------------------------------------------------------------------------ ##
    
print("\n")
print('#---------------------------------------------------------------------#')
print('#                ONCOGENICITY GUIDELINES INTERPRETATION               #')
print('#---------------------------------------------------------------------#\n')
print("\n") 

for filename in os.listdir(input_data_dir):
        
    file = os.path.join(input_data_dir, filename)
    
    if os.path.isfile(file) and filename.endswith("_prediction.xlsx"):
        
        name = os.path.basename(file).split("_prediction.xlsx")[0]
        
        print("\n")
        print('#---------------------------------------------------------------------#')
        print('#                SAMPLE PROCESSED: ' + name + '                       #')
        print('#---------------------------------------------------------------------#\n')
        print("\n")
        
        # -------------------------------------------------------------------------- ##
        
        # Read in data    
        data_init = pd.read_excel(file, na_filter=False)
        
        # criteria                                                  
        criteria_colnames = ["OVS1", "OS1", "OS2", "OS3", 
                             "OM1", "OM2", "OM3", "OM4", 
                             "OP1", "OP2", "OP3", "OP4",
                             "SBVS1", "SBS1", "SBS2", "SBP1", "SBP2"]
        
        # Data frame of triggered criteria
        # "yes" if the criterion is triggered by the variant
        # "no" otherwise
        criteria_df = pd.DataFrame(columns=criteria_colnames)
        
        criteria_colnames_points = ["OVS1_p", "OS1_p", "OS2_p", "OS3_p", 
                                    "OM1_p", "OM2_p", "OM3_p", "OM4_p",
                                    "OP1_p", "OP2_p", "OP3_p", "OP4_p",
                                    "SBVS1_p", "SBS1_p", "SBS2_p", "SBP1_p", "SBP2_p"]
        
        # Data frame of points associated to the triggered criteria
        points_df = pd.DataFrame(columns=criteria_colnames_points)
        
        # Dictionary of quantitative and qualitative interpretation
        classification_dict = {}
        
        # Data frame with classification (quantitative and qualitative)
        # and triggered criteria.
        triggered_criteria_df = pd.DataFrame(columns = ["Criteria"])
        
        # -------------------------------------------------------------------------- ##

        # Iterate over the variants of the input data  
        for index, row in data_init.iterrows():
            
            ## Variant features employed for criteria triggering
            chrom = data_init.loc[index, "CHROM"]
            
            pos = data_init.loc[index, "POS"]
            
            ref_allele = data_init.loc[index, "REF"]          # e.g. "C"
            
            alt_allele = data_init.loc[index, "ALT"]          # e.g. "T"
            
            gene = data_init.loc[index,"SYMBOL"]              # e.g. "MTOR"
            
            consequence = data_init.loc[index, "Consequence"] # e.g. "missense_variant"
            
            hgvsc = data_init.loc[index, "HGVSc"]             # e.g. "NM_004958.4:c.7255G>A"
            
            hgvsp = data_init.loc[index, "HGVSp"]             # e.g. "NP_004949.1:p.Glu2419Lys"
            
            #amino_acids = data_init.loc[index, "Amino_acids"] # e.g. "E/K"
            
            transcript = hgvsc.split(':')[0]                  # e.g. "NM_004958.4"
            
            # -------------------------------------------------------------------------- ##
            
            if hgvsc != "":
                
                cdna = hgvsc.split(':')[1]                    
                # e.g. "c.7255G>A"
                
            else:
                
                cdna = ""
            
            transcript_without_version = transcript.split(".")[0]    
            # e.g. "NM_004958" 
            # isolate the transcript from the transcript version
        
            # -------------------------------------------------------------------------- ##
            
            # %3D is the way VEP uses to indicate variants 
            # causing an alternate AA identical to the reference AA
            # (synonymous variant at the protein level)
            if "%3D" in hgvsp:
                
                # given prot_change "Val413%3D"
                prot_change = hgvsp.split('p.')[1]
                
                ref = prot_change[:3]  # "Val"
                alt = ref              # "Val"
                
                # ref and alt AA in single letters
                ref_amino_acid = amino_dict[ref]  # "V"
                alt_amino_acid = amino_dict[alt]  # "V"
                
                # replace %3D% with ref
                prot_change = prot_change.replace("%3D", ref)
                
                # variant identifier
                var_identifier = gene + ":" + "p." + prot_change
                
                amino_acid_position = '-'.join(re.findall("\d+", prot_change))
                
                # iterate over amino_dict keys and replace with values
                for key in amino_dict.keys():
                    if key in prot_change:
                        prot_change = prot_change.replace(key, amino_dict[key])
                
                # position in the protein in which the change occurs
                ref_amino_acid_and_position = prot_change[:-1]
                
                # to access CGI database 
                var_id = gene + ":" + "p." + prot_change
                
            # -------------------------------------------------------------------------- ##    
            
            elif hgvsp != "":
                
                prot_change = hgvsp.split('p.')[1]           
                
                if "_" and "delins" in prot_change:
                    
                    # given prot_change 'Gln86_Glu87delinsArgLys'                            
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = ""
                    alt = ""
                    
                    # ref and alt AA in single letters
                    ref_amino_acid = ""
                    alt_amino_acid = ""
                      
                    amino_acid_position = '-'.join(re.findall("\d+", prot_change))
                    
                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                            
                    ref_amino_acid_and_position = ""
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change
                    
                # -------------------------------------------------------------------------- ##
                
                elif "_" and "ins" in prot_change:
                    
                    # given prot_change                             
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = ""
                    alt = ""
                    
                    # ref and alt AA in single letters
                    ref_amino_acid = ""
                    alt_amino_acid = ""
                      
                    amino_acid_position = '-'.join(re.findall("\d+", prot_change))
                    
                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                            
                    ref_amino_acid_and_position = ""
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change
                
                # -------------------------------------------------------------------------- ##
                 
                elif "_" and "del" in prot_change:
                    
                    # given prot_change                            
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = ""
                    alt = ""
                    
                    # ref and alt AA in single letters
                    ref_amino_acid = ""
                    alt_amino_acid = ""
                      
                    amino_acid_position = '-'.join(re.findall("\d+", prot_change))
                    
                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                            
                    ref_amino_acid_and_position = ""
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change
                
                
                # -------------------------------------------------------------------------- ##
                
                elif "_" and "dup" in prot_change:
                    
                    # given prot_change 'Leu14_Leu16dup'
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = ""
                    alt = ""
                    
                    # ref and alt AA in single letters
                    ref_amino_acid = ""
                    alt_amino_acid = ""
                      
                    amino_acid_position = '-'.join(re.findall("\d+", prot_change))
                    
                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                            
                    ref_amino_acid_and_position = ""    
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change
                
                # -------------------------------------------------------------------------- ##
                
                elif "_" in prot_change:
                    
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = ""
                    alt = ""
                    
                    # ref and alt AA in single letters
                    ref_amino_acid = ""
                    alt_amino_acid = ""
                    
                    # given prot_change 
                    amino_acid_position = '-'.join(re.findall("\d+", prot_change))
                    
                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                            
                    ref_amino_acid_and_position = ""    
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change
                    
                # -------------------------------------------------------------------------- ##
                
                elif "ext" in prot_change:
                    
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = ""
                    alt = ""
                    
                    # ref and alt AA in single letters
                    ref_amino_acid = ""
                    alt_amino_acid = ""
                    
                    # given prot_change 'Ter730LeuextTer205'                        
                    amino_acid_position = ''.join(re.findall("\d+", prot_change.split('ext')[0]))
                    
                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                    
                    ref_amino_acid_and_position = ""    
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change
                
                # -------------------------------------------------------------------------- ##
                
                elif "dup" in prot_change:
                    
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = prot_change[:3]  # "Leu"
                    alt = ""
                    
                    # ref and alt AA in single letters 
                    ref_amino_acid = amino_dict[ref]  # "L"
                    alt_amino_acid = ""
                    
                    # given prot_change 'Leu15dup'                        
                    amino_acid_position = ''.join(re.findall("\d+", prot_change))
                    
                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                            
                    # position and ref AA in the prot in which the change occurs
                    ref_amino_acid_and_position = ref_amino_acid + amino_acid_position
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change
                    
                # -------------------------------------------------------------------------- ##
                
                elif "del" in prot_change:
                    
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = prot_change[:3]  # "Thr"
                    alt = ""
                    
                    # ref and alt AA in single letters 
                    ref_amino_acid = amino_dict[ref]  # "T"
                    alt_amino_acid = ""
                    
                    # given prot_change 'Thr949del'                        
                    amino_acid_position = ''.join(re.findall("\d+", prot_change))

                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                            
                    # position and ref AA in the prot in which the change occurs
                    ref_amino_acid_and_position = ref_amino_acid + amino_acid_position
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change 
                
                # -------------------------------------------------------------------------- ##
                
                elif "fs" in prot_change:
                    
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = prot_change[:3]  # "Gly"
                    alt = ""
                    
                    # ref and alt AA in single letters 
                    ref_amino_acid = amino_dict[ref]  # "G"
                    alt_amino_acid = ""
                                                
                    # given prot_change 'Gly645ValfsTer58'                        
                    amino_acid_position = ''.join(re.findall("\d+", prot_change.split('fs')[0]))

                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                    
                    # position and ref AA in the prot in which the change occurs
                    ref_amino_acid_and_position = ref_amino_acid + amino_acid_position
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change
                
                # -------------------------------------------------------------------------- ##
                
                elif "Ter" in prot_change:
                    
                    # prot_change is 'Arg651Ter'
                    var_identifier = gene + ":" + "p." + prot_change
                    
                    ref = prot_change[:3]  # "Arg"
                    alt = prot_change[-3:] # "Ter"
                    
                    # ref and alt AA in single letters 
                    ref_amino_acid = amino_dict[ref]  # "R"
                    alt_amino_acid = alt  # "Ter"
                    
                    # given prot_change 'Arg651Ter'                        
                    amino_acid_position = ''.join(re.findall("\d+", prot_change.split('Ter')[0]))
                    
                    # iterate over amino_dict keys and replace with values
                    for key in amino_dict.keys():
                        if key in prot_change:
                            prot_change = prot_change.replace(key, amino_dict[key])
                    
                    # position and ref AA in the prot in which the change occurs
                    ref_amino_acid_and_position = ref_amino_acid + amino_acid_position
                
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change
                
                # -------------------------------------------------------------------------- ##
                
                elif "?" in prot_change:
                    
                    ref = prot_change[:3]  # "M"
                    alt = ""
                    
                    # ref and alt AA in single letters 
                    ref_amino_acid = amino_dict[ref]  # "G"
                    alt_amino_acid = ""
                    
                    # for some variants, like NP_065812.1:p.Glu2?, is not 
                    # possible to retrieve the alt AA from alt_amino_acid
                    # because an X is reported (e.g. M/X).
                    try:
                        for key in amino_dict.keys():
                            if key in prot_change:
                                prot_change = prot_change.replace(key, amino_dict[key])
                                
                        # to access CGI database 
                        var_id = gene + ":" + "p." + prot_change
                        
                        var_identifier = gene + ":" + "p." + prot_change

                    except:
                        var_identifier = gene + ":" + "p." + prot_change
                        
                        # to access CGI database 
                        var_id = var_identifier
                    
                    # given prot_change 'Glu2?'                        
                    amino_acid_position = ''.join(re.findall("\d+", prot_change))
                    
                    # position and ref AA in the prot in which the change occurs
                    ref_amino_acid_and_position = ref_amino_acid + amino_acid_position

                # -------------------------------------------------------------------------- ##

                else:
                    
                    # given prot_change "Ala117Thr"
                    var_identifier = gene + ":" + "p." + prot_change      
                    
                    ref = prot_change[:3]  # "Ala"
                    alt = prot_change[-3:] # "Thr"
                    
                    # ref and alt AA in single letters 
                    ref_amino_acid = amino_dict[ref]  # "A"
                    alt_amino_acid = amino_dict[alt]  # "T"
                    
                    amino_acid_position = ''.join(re.findall("\d+", prot_change))
                    
                    prot_change = prot_change.replace(ref, ref_amino_acid)
                    prot_change = prot_change.replace(alt, alt_amino_acid)
                    
                    # position and ref AA in the prot in which the change occurs
                    ref_amino_acid_and_position = prot_change[:-1]
                    
                    # to access CGI database 
                    var_id = gene + ":" + "p." + prot_change

            # -------------------------------------------------------------------------- ##
            
            elif (hgvsc != "") and (hgvsp == ""):
                
                var_identifier = gene + ":" + hgvsc.split(':')[1]     
                
                ref = alt = ""
                
                # ref and alt AA in single letters 
                ref_amino_acid = alt_amino_acid = ""
                
                amino_acid_position = ""
                
                prot_change = ""
                
                # position and ref AA in the prot in which the change occurs
                ref_amino_acid_and_position = ""
                
                var_id = var_identifier
                
            else:
                
                var_identifier = gene + ":" + str(pos) + ref_allele + ">" + alt_allele
                
                ref = ""
                alt = ""
                
                # ref and alt AA in single letters
                ref_amino_acid = ""
                alt_amino_acid = ""
                
                amino_acid_position = ""
                
                prot_change = ""
                        
                ref_amino_acid_and_position = ""
                
                # to access CGI database 
                var_id = gene + ":" + "g." + str(pos) + ref_allele + ">" + alt_allele
            
                
            # -------------------------------------------------------------------------- ##
            
            # AF from gnomAD exome populations. 
            # 5 general continental populations: 
            # African, East Asian, European (non-Finnish), Latino and South Asian. 
            gnomADe_freq = data_init.loc[index,"gnomADe_AF"]
            gnomADe_AFR_AF = data_init.loc[index, "gnomADe_AFR_AF"]
            gnomADe_EAS_AF = data_init.loc[index, "gnomADe_EAS_AF"]
            gnomADe_NFE_AF = data_init.loc[index, "gnomADe_NFE_AF"]
            gnomADe_AMR_AF = data_init.loc[index, "gnomADe_AMR_AF"]
            gnomADe_SAS_AF = data_init.loc[index, "gnomADe_SAS_AF"]
            
            # AF from gnomAD genome populations. 
            # 5 general continental populations: 
            # African, East Asian, European (non-Finnish), Latino and South Asian. 
            gnomADg_freq = data_init.loc[index,"gnomADg_AF"]
            gnomADg_AFR_AF = data_init.loc[index, "gnomADg_AFR_AF"]
            gnomADg_EAS_AF = data_init.loc[index, "gnomADg_EAS_AF"]
            gnomADg_NFE_AF = data_init.loc[index, "gnomADg_NFE_AF"]
            gnomADg_AMR_AF = data_init.loc[index, "gnomADg_AMR_AF"]
            gnomADg_SAS_AF = data_init.loc[index, "gnomADg_SAS_AF"]
            
            # List of AF in gnomAD exome populations in any 5 general continental
            # populations . 
            gnomADe_freq_pop = [gnomADe_AFR_AF, 
                                gnomADe_EAS_AF, 
                                gnomADe_NFE_AF,
                                gnomADe_AMR_AF, 
                                gnomADe_SAS_AF]
            
            # List of AF in gnomAD genome populations in any 5 general continental
            # populations. 
            gnomADg_freq_pop = [gnomADg_AFR_AF, 
                                gnomADg_EAS_AF, 
                                gnomADg_NFE_AF, 
                                gnomADg_AMR_AF, 
                                gnomADg_SAS_AF]
            
            # -------------------------------------------------------------------------- ##

            # In silico algorithm predictions 
            phylop_score = data_init.loc[index, "phyloP100way_vertebrate_rankscore"]
            if phylop_score == "":
                pass
            else:
                phylop_score = float(phylop_score)
                
            phastCons_score = data_init.loc[index, "phastCons100way_vertebrate_rankscore"]
            if phastCons_score == "":
                pass
            else:
                phastCons_score = float(phastCons_score)
                
            spliceAI_score = data_init.loc[index, "SpliceAI_cutoff"]
            
            # -------------------------------------------------------------------------- ##

            
            # -------------------------------------------------------------------------- ##
            ####### CRITERIA FOR EVIDENCE OF ONCOGENICITY OF SOMATIC VARIANTS #######
            
            ## Category: VERY STRONG (8 points)
            ## Evidence: OVS1
            
            ## Criteria: 
            ## Null variant (nonsense, frameshift, 
            ## canonical ±1 or 2 splice sites, initiation
            ## codon, single-exon or multixon deletion) 
            ## in a bona fide tumor supressor gene. 
            
            # List of null variants for OVS1 criteria 
            null_variants = ["stop_gained",
                             "stop_lost",
                             "start_lost",
                             "frameshift_variant",
                             "splice_donor_variant",
                             "splice_acceptor_variant"]
            
            # Check that the gene is not in TSGs list
            if gene not in tsg_bona_fide_list:
                
                criteria_df.loc[var_identifier, "OVS1"] = "no"
                points_df.loc[var_identifier, "OVS1_p"] = 0
            
            # The gene is in TSGs list    
            else:
                
                # check if the variant is a null variant 
                if any(cons in consequence for cons in null_variants):
                    
                  criteria_df.loc[var_identifier, "OVS1"] = 'yes'
                  points_df.loc[var_identifier, "OVS1_p"] = 8 
                                    
                elif 'splice' in consequence:
                    
                    if (("+") or ("-")) in hgvsc:
                        
                        if "+" in hgvsc:
                            
                            position_and_change = hgvsc.split('+')[1]
                            position_changed = re.match(r'^\d{1,2}', position_and_change).group(0)
                            
                            # with this condition we are considering 4 bases (whose 2 are 
                            # in the canonical sites and the additional 2 bases are those
                            # in the intron - either 5' donor or 3' acceptor)
                            if float(position_changed) < 5:
                                criteria_df.loc[var_identifier, "OVS1"] = 'yes'
                                points_df.loc[var_identifier, "OVS1_p"] = 8
                                
                            else:
                                criteria_df.loc[var_identifier, "OVS1"] = 'no'
                                points_df.loc[var_identifier, "OVS1_p"] = 0
                                 
                        elif "-" in hgvsc:
                            
                            position_and_change = hgvsc.split('-')[1]
                            position_changed = re.match(r'^\d{1,2}', position_and_change).group(0)
                            
                            if  float(position_changed) < 5:
                                criteria_df.loc[var_identifier, "OVS1"] = 'yes'
                                points_df.loc[var_identifier, "OVS1_p"] = 8
                                
                            else:
                                criteria_df.loc[var_identifier, "OVS1"] = 'no'
                                points_df.loc[var_identifier, "OVS1_p"] = 0
                                
                    else:

                        # if no + or - appears in the HGVSc
                        criteria_df.loc[var_identifier, "OVS1"] = 'no'
                        points_df.loc[var_identifier, "OVS1_p"] = 0
                        
                else:
                    
                    criteria_df.loc[var_identifier, "OVS1"] = "no"
                    points_df.loc[var_identifier, "OVS1_p"] = 0  
                    
                    
            # Check if the gene mutated is reported in MutSpliceDB 
            # in case the OVS1 criterion was not triggered               
            if (criteria_df.loc[var_identifier, 'OVS1'] == "no") and gene in mut_splice_dict.keys():
                
                # Check if the splicing variant is in the DB
                gene_vars = [value for key, value in mut_splice_dict[gene].items() if transcript_without_version == key]
                    
                if len(gene_vars) > 0 and (cdna in gene_vars[0]):
                    
                    criteria_df.loc[var_identifier, "OVS1"] = 'yes'
                    points_df.loc[var_identifier, "OVS1_p"] = 8
                    
                else:
                    criteria_df.loc[var_identifier, "OVS1"] = 'no'
                    points_df.loc[var_identifier, "OVS1_p"] = 0
                
            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            ## Category: STRONG (4 points)
            ## Evidence: OS1
            
            ## Criteria: 
            ## Same amino acid change as a previously established 
            ## oncogenic variant (using this standard) regardless of 
            ## nucleotide change. Example: Val->Leu caused by
            ## either G>C or G>T in the same codon. 
            
            # To evaluate this criterion, the Cancer Genome Interpreter
            # database is used as an external set of previously 
            # identified variants. 
            
            if hgvsp != "":
                    
                # Check if the variant is in CGI 
                if var_id not in onco_var_cgi.keys():
                        
                        criteria_df.loc[var_identifier, 'OS1'] = 'no'
                        points_df.loc[var_identifier, 'OS1_p'] = 0
                        
                else:
                    
                    criteria_df.loc[var_identifier, 'OS1'] = 'yes'
                    points_df.loc[var_identifier, 'OS1_p'] = 4
                    
            else:

                #####################################
                # CHECK FOR VARIANTS WITHOUT p.CHANGE
                
                if var_identifier not in onco_var_cgi.keys():
                    
                    criteria_df.loc[var_identifier, 'OS1'] = 'no'
                    points_df.loc[var_identifier, 'OS1_p'] = 0
                    
                else:
                    
                    criteria_df.loc[var_identifier, 'OS1'] = 'yes'
                    points_df.loc[var_identifier, 'OS1_p'] = 4
                                
            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: STRONG (4 points)
            ## Evidence: OS2
            
            ## Criteria: 
            ## Well-established in vitro or in vivo functional studies, supportive of
            ## an oncogenic effect of the variant.  
            
            # The ClinVar_review_status contains the level of review supporting the 
            # assertion of clinical significance for the variant reported by ClinVar.
            review_status = data_init.loc[index, "ClinVar_germline_ReviewStatus"]
            clinvar_classification = data_init.loc[index, "ClinVar_germline"]
            
            # List of accepted review status from ClinVar to meet the OS1 criterion 
            review_level_list = ['practice guideline',
                                 'reviewed by expert panel',
                                 'criteria provided, single submitter',
                                 'criteria provided, multiple submitters, no conflicts']
            
            if ("Pathogenic" in clinvar_classification) and (review_status in review_level_list):
                
                criteria_df.loc[var_identifier, 'OS2'] = 'yes'
                points_df.loc[var_identifier, 'OS2_p'] = 4
                
            elif ("Likely pathogenic" in clinvar_classification) and (review_status in review_level_list):
                
                criteria_df.loc[var_identifier, 'OS2'] = 'yes'
                points_df.loc[var_identifier, 'OS2_p'] = 4
                
            elif (clinvar_classification in ClinicalSignificance_list) and (review_status in review_level_list):
                
                criteria_df.loc[var_identifier, 'OS2'] = 'yes'
                points_df.loc[var_identifier, 'OS2_p'] = 4
                        
            else:
                
                criteria_df.loc[var_identifier, 'OS2'] = 'no'
                points_df.loc[var_identifier, 'OS2_p'] = 0
                    
            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: STRONG (4 points)
            ## Evidence: OS3
            
            ## Criteria: 
            ## Located in one of the hotspots in cancerhotspots with at least 50 samples 
            ## with a somatic variant at the same amino acid position, and the same amino
            ## acid change count in cancerhotspots in at least 10 samples.   
            
            # Define indel consequence
            indels = ["inframe_insertion", 
                      "inframe_deletion"]
            
            # If the OS1 criterion is applicable, then 
            if criteria_df.loc[var_identifier, 'OS1'] == "yes":
                
                # This avoid verifying the presence of the variant in COSMIC
                res_cancerhotspots = "True"
                
                criteria_df.loc[var_identifier, 'OS3'] = "no"
                points_df.loc[var_identifier, 'OS3_p'] = 0
                
            else:
                # For inframe deletion or insertion variants the in_frame_indel_dict is 
                # investigated 
                if any(s for s in indels if s in consequence):
                    
                    # Check if the gene is reported in cancerhotspots
                    if gene in in_frame_indel_dict.keys():
                        
                        gene_records_inframe_indel_dict = in_frame_indel_dict[gene]  
                        
                        for key in gene_records_inframe_indel_dict.keys():
                                                        
                            if not '-' in key:
                                
                                if amino_acid_position == key:
                                    
                                    res_variant_inframe_indel = [value for key, value in gene_records_inframe_indel_dict[key].items() if prot_change == key]
                                    
                                    if len(res_variant_inframe_indel) == 0:
                                        
                                        # the residue is not in cancerhotspots
                                        res_cancerhotspots = "False" 
                                            
                                        criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                        points_df.loc[var_identifier, 'OS3_p'] = 0
                                            
                                    else:
                                        
                                        # Check if the requirements from the guidelines are satisfied
                                        # by cancerhotspots records
                                        if (sum(list(map(int, gene_records_inframe_indel_dict[key].values()))) >= 50) and (int(res_variant_inframe_indel[0]) >= 10):
                                            
                                            res_cancerhotspots = "True"
                                            
                                            criteria_df.loc[var_identifier, 'OS3'] = 'yes'
                                            points_df.loc[var_identifier, 'OS3_p'] = 4
                                            break;
                                                
                                        else:
                                            
                                            # the residue is in cancerhotspots but not 
                                            # satisfying OS3 requirements. This avoid 
                                            # verifying the presence of the variant in 
                                            # COSMIC, given that we give priority to
                                            # cancerhotspots.org
                                            res_cancerhotspots = "True"
                                            
                                            criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                            points_df.loc[var_identifier, 'OS3_p'] = 0
                                
                                else:
                                    
                                    # the residue is not in cancerhotspots and 
                                    # we want to check the presence of the 
                                    # variant in COSMIC
                                    res_cancerhotspots = "False"
                                    
                                    criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                    points_df.loc[var_identifier, 'OS3_p'] = 0
                                    
                            else:
                                
                                ## create range value 
                                min_range, max_range = key.split('-')                                    
                                
                                amino_acid_pos_list = re.findall("\d+", amino_acid_position)

                                if len(amino_acid_pos_list) >1:
                                    
                                    if (int(amino_acid_pos_list[0]) >= int(min_range)) and (int(amino_acid_pos_list[1]) <= int(max_range)):
                                        
                                        res_variant_inframe_indel = [value for key, value in gene_records_inframe_indel_dict[key].items() if prot_change == key]
                                        
                                        if len(res_variant_inframe_indel) == 0:
                                            
                                            # the residue is not in cancerhotspots
                                            res_cancerhotspots = "False" 
                                                
                                            criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                            points_df.loc[var_identifier, 'OS3_p'] = 0
                                                
                                        else:
                                                
                                           # Check if the requirements from the guidelines are satisfied
                                           # by cancerhotspots records
                                           if (sum(list(map(int, gene_records_inframe_indel_dict[key].values()))) >= 50) and (int(res_variant_inframe_indel[0]) >= 10):
                                               
                                               # the residue is in cancerhotspots
                                               res_cancerhotspots = "True"
                                               
                                               criteria_df.loc[var_identifier, 'OS3'] = 'yes'
                                               points_df.loc[var_identifier, 'OS3_p'] = 4
                                               break;
                                                   
                                           else:
                                               
                                               # the residue is in cancerhotspots but not 
                                               # satisfying OS3 requirements. This avoid 
                                               # verifying the presence of the variant in 
                                               # COSMIC, given that we give priority to
                                               # cancerhotspots.org
                                               res_cancerhotspots = "True" 
                                               
                                               criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                               points_df.loc[var_identifier, 'OS3_p'] = 0
                                                    
                                    else:
                                        
                                        # the residue is not in cancerhotspots
                                        res_cancerhotspots = "False" 
                                        
                                        criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                        points_df.loc[var_identifier, 'OS3_p'] = 0
                                        
                                                    
                                else:
                                    
                                    if int(min_range) <= int(amino_acid_pos_list[0]) <= int(max_range):
                                        
                                        res_variant_inframe_indel = [value for key, value in gene_records_inframe_indel_dict[key].items() if prot_change == key]
                                        
                                        if len(res_variant_inframe_indel) == 0:
                                            
                                            # the residue is not in cancerhotspots
                                            res_cancerhotspots = "False" 
                                                
                                            criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                            points_df.loc[var_identifier, 'OS3_p'] = 0
                                                
                                        else:
                                            
                                            if (sum(list(map(int, gene_records_inframe_indel_dict[key].values()))) >= 50) and (int(res_variant_inframe_indel[0]) >= 10):
                                                
                                                res_cancerhotspots = "True" 
                                                
                                                criteria_df.loc[var_identifier, 'OS3'] = 'yes'
                                                points_df.loc[var_identifier, 'OS3_p'] = 4
                                                break;
                                            
                                            else:
                                                
                                                # the residue is in cancerhotspots but not 
                                                # satisfying OS3 requirements. This avoid 
                                                # verifying the presence of the variant in 
                                                # COSMIC, given that we give priority to
                                                # cancerhotspots.org 
                                                res_cancerhotspots = "True" 
                                                
                                                criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                                points_df.loc[var_identifier, 'OS3_p'] = 0
                                    
                                    else:
                                        
                                        # the residue is not in cancerhotspots
                                        res_cancerhotspots = "False" 
                                        
                                        criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                        points_df.loc[var_identifier, 'OS3_p'] = 0

                    else:
                        
                        # there are no residues in cancerhotspots 
                        # belonging to the evaluated gene
                        res_cancerhotspots = "False"
                        
                        criteria_df.loc[var_identifier, 'OS3'] = 'no'
                        points_df.loc[var_identifier, 'OS3_p'] = 0
                    
                else:
                    
                    # Otherwise, the single_residue_dict is considered.
                    # Check if the gene is reported in cancerhotspots
                    if gene in single_residue_dict.keys():
                        
                        res_var_single_residue = [value for key, value in single_residue_dict[gene].items() if ref_amino_acid_and_position == key]
                        
                        # Check if the variant evaluated is reported as residue in cancerhotspots
                        # If not
                        if len(res_var_single_residue) == 0: 
                            
                            # the residue is not in cancerhotspots
                            res_cancerhotspots = "False" 
                            
                            criteria_df.loc[var_identifier, 'OS3'] = 'no'
                            points_df.loc[var_identifier, 'OS3_p'] = 0
                        
                        # If yes
                        else:            
                            
                            res_var = [value for key, value in res_var_single_residue[0].items() if alt_amino_acid == key]
                            
                            # Satisfying OS3 requirements 
                            if (len(res_var) != 0) and (sum(list(map(int, res_var_single_residue[0].values() ))) >= 50) and (int(res_var[0]) >= 10):
                                
                                # the residue is in cancerhotspots
                                res_cancerhotspots = "True" 
                                
                                criteria_df.loc[var_identifier, 'OS3'] = 'yes'
                                points_df.loc[var_identifier, 'OS3_p'] = 4
                            
                            # Not satisfying OS3 requirements
                            else:
                                
                                # the residue is in cancerhotspots but not 
                                # satisfying OS3 requirements. This avoid 
                                # verifying the presence of the variant in 
                                # COSMIC, given that we give priority to
                                # cancerhotspots.org
                                res_cancerhotspots = "True" 
                                
                                criteria_df.loc[var_identifier, 'OS3'] = 'no'
                                points_df.loc[var_identifier, 'OS3_p'] = 0 
                                    
                    else:
                        
                        # there are no residues in cancerhotspots
                        # belonging to the evaluated gene
                        res_cancerhotspots = "False" 
                        
                        criteria_df.loc[var_identifier, 'OS3'] = 'no'
                        points_df.loc[var_identifier, 'OS3_p'] = 0
                        
            # If the residue is not present at all in cancerhotspots.org, the
            # COSMIC database is considered
            if (res_cancerhotspots == "False") and (criteria_df.loc[var_identifier, 'OS3'] == "no") and (criteria_df.loc[var_identifier, 'OS1'] == "no"):
                
                if "-" in amino_acid_position:
                    ref_AA, alt_AA = amino_acid_position.split('-')
                    
                else:
                    ref_AA = amino_acid_position
                
                if "p." in var_id:
                    
                    if (gene in cosmic_all_dict.keys()) and (ref_AA in cosmic_all_dict[gene].keys()):
                        
                        selection = cosmic_all_dict[gene][ref_AA]
                        
                        # How many samples are mutated in the AA position
                        cosmic_entries_AA_pos = sum([selection[key1][key2] for key1 in selection.keys() for key2 in selection[key1].keys()])
                        
                        # How many samples are mutated with the same AA change
                        cosmic_entries_same_AA_change = sum([selection[key3][key4] for key3 in selection.keys() for key4 in selection[key3].keys() if key3 == var_id.split(":")[1]])
                        
                    else:
                        cosmic_entries_AA_pos = float('nan')
                        cosmic_entries_same_AA_change = float('nan')
                
                
                elif ("c." in var_id) or ("g." in var_id):
                    
                    # strings to query COSMIC Cancer Gene Mutations database 
                    key_AA_pos = str(chrom) + ":g." + str(pos)
                    
                    key_AA_change = str(chrom) + ":" + data_init.loc[index, "HGVSg"].split(':')[1]
                    
                    # How many samples are mutated in the AA position
                    cosmic_entries_AA_pos = sum([len(cosmic_hgvsg_dict[key1]) for key1 in cosmic_hgvsg_dict.keys() if key1.startswith(key_AA_pos)])
                        
                    # How many samples are mutated with the same AA change
                    cosmic_entries_same_AA_change = sum([len(cosmic_hgvsg_dict[key2]) for key2 in cosmic_hgvsg_dict.keys() if key2 == key_AA_change])
                    
                else:
                    cosmic_entries_AA_pos = float('nan')
                    cosmic_entries_same_AA_change = float('nan')

                    
                if (cosmic_entries_AA_pos >= 50) and (cosmic_entries_same_AA_change >= 10):
                        
                    criteria_df.loc[var_identifier, 'OS3'] = 'yes'
                    points_df.loc[var_identifier, 'OS3_p'] = 4
                
                else:
                    criteria_df.loc[var_identifier, 'OS3'] = 'no'
                    points_df.loc[var_identifier, 'OS3_p'] = 0
                    
            else:
                cosmic_entries_AA_pos = float('nan')
                cosmic_entries_same_AA_change = float('nan')
                
            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: MODERATE (2 points)
            ## Evidence: OM1
            
            ## Criteria: 
            ## Located in a critical and well-established part of a functional domain 
            ## (e.g. active site of an enzyme)
            
            # If the OS1 or OS3 criterion is applicable, then 
            if (criteria_df.loc[var_identifier, 'OS1'] == "yes") or (criteria_df.loc[var_identifier, 'OS3'] == "yes"):
                
                criteria_df.loc[var_identifier, 'OM1'] = "no"
                points_df.loc[var_identifier, 'OM1_p'] = 0
                
            else:
                
                if (gene not in domains_dict.keys()) or (hgvsp == ""):
                    
                    criteria_df.loc[var_identifier, 'OM1'] = "no"
                    points_df.loc[var_identifier, 'OM1_p'] = 0
                    
                else:
                    
                    domains_gene = domains_dict[gene]
                    
                    # protein change caused by the variant
                    if "-" not in amino_acid_position: 
                    
                        for domain in domains_gene:
                            
                            start_domain, end_domain = domain.split('-')
                            start_domain = ''.join(re.findall("\d+", start_domain))
                            end_domain = ''.join(re.findall("\d+", end_domain))
                            
                            if int(start_domain) <= int(amino_acid_position) <= int(end_domain):
                                
                                criteria_df.loc[var_identifier, 'OM1'] = "yes"
                                points_df.loc[var_identifier, 'OM1_p'] = 2
                                break;
                            
                            else:
                                
                                criteria_df.loc[var_identifier, 'OM1'] = "no"
                                points_df.loc[var_identifier, 'OM1_p'] = 0
                                
                    else:
                       
                       for domain in domains_gene:
                           
                           start_domain, end_domain = domain.split('-')
                           start_domain = ''.join(re.findall("\d+", start_domain))
                           end_domain = ''.join(re.findall("\d+", end_domain))
                           
                           start_end_aminoacid_change = re.findall("\d+", amino_acid_position)
                           
                           # the domain of the protein is completely
                           # affected by the variant
                           if (int(start_end_aminoacid_change[0]) >= int(start_domain)) and (int(start_end_aminoacid_change[1]) >= int(end_domain)):
                               
                               criteria_df.loc[var_identifier, 'OM1'] = "yes"
                               points_df.loc[var_identifier, 'OM1_p'] = 2
                               break;
                           
                           # the start position of the variant affects
                           # the domain of the protein
                           elif int(start_end_aminoacid_change[0]) in range(int(start_domain), int(end_domain)):
                            
                               criteria_df.loc[var_identifier, 'OM1'] = "yes"
                               points_df.loc[var_identifier, 'OM1_p'] = 2
                               break;
                           
                           # the end position of the variant affects
                           # the domain of the protein
                           elif int(start_end_aminoacid_change[1]) in range(int(start_domain), int(end_domain)):
                            
                               criteria_df.loc[var_identifier, 'OM1'] = "yes"
                               points_df.loc[var_identifier, 'OM1_p'] = 2
                               break;
                           
                           else:
                               
                               criteria_df.loc[var_identifier, 'OM1'] = "no"
                               points_df.loc[var_identifier, 'OM1_p'] = 0

                
            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: MODERATE (2 points)
            ## Evidence: OM2
            
            ## Criteria: 
            ## Protein length changes as a result of in-frame deletions/insertions in a  
            ## known oncogene or tumor supressor gene or stop-loss variants in a known 
            ## tumor supressor gene.
            
            # If the OVS1 criterion is applicable, then 
            if criteria_df.loc[var_identifier, 'OVS1'] == "yes":
                
                criteria_df.loc[var_identifier, 'OM2'] = "no"
                points_df.loc[var_identifier, 'OM2_p'] = 0
                
            else:
                
                if any(cons in consequence for cons in indels) and ((gene in oncogenes_list) or (gene in tsg_list)):
                            
                    criteria_df.loc[var_identifier, 'OM2'] = "yes"
                    points_df.loc[var_identifier, 'OM2_p'] = 2
                            
                # or consequence is "stop_lost" and in a known TSGs  
                elif ("stop_lost" in consequence) and (gene in tsg_list):
                    
                    criteria_df.loc[var_identifier, 'OM2'] = "yes"
                    points_df.loc[var_identifier, 'OM2_p'] = 2
                    
                else:
                    criteria_df.loc[var_identifier, 'OM2'] = "no"
                    points_df.loc[var_identifier, 'OM2_p'] = 0                        
                    

            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: MODERATE (2 points)
            ## Evidence: OM3
            
            ## Criteria: 
            ## Missense variant at an amino acid residue where a different missense variant 
            ## determined to be oncogenic (using this standard) has been documented. Example:
            ## p.Arg156His is oncogenic; now you observe p.Arg156Cys. Amino acid difference
            ## from reference amino acid should be grater or at least approximately the same
            ## as for missense change determined to be oncogenic.
            
            # Here, the exception to OM3 application due to OS1 activation was 
            # removed. Because the OM3 criterion needs OS1 to be applicable 
            # in order to be evaluated.
            
            # If the OS3 or OM1 criterion is applicable, then 
            if ((criteria_df.loc[var_identifier, 'OS3'] == "yes") or (criteria_df.loc[var_identifier, 'OM1'] == "yes")):
                
                criteria_df.loc[var_identifier, 'OM3'] = "no"
                points_df.loc[var_identifier, 'OM3_p'] = 0
                
            else:
                
                # check if the variant is missense
                if "missense_variant" in consequence:
                    
                    # variant identifier without the alternate allele is 
                    # used to search in CGI
                    var_id_without_alt_AA = var_id.replace(alt_amino_acid, "")
                    
                    # Check if the variant was previously classified as 
                    # oncogenic variant (i.e. in CGI database) 
                    # Example: p.Arg156His is oncogenic; now you observe 
                    # p.Arg156Cys.
                    # Grantham distance is for missense variants only 
                    res = [key for key in onco_var_cgi.keys() if (var_id_without_alt_AA in key) and (key != var_id)]                       
                    
                    if res == []:
                        
                        criteria_df.loc[var_identifier, 'OM3'] = "no"
                        points_df.loc[var_identifier, 'OM3_p'] = 0
                        
                    else:
                        
                        for item in res:
                            
                            original_alt_amino = item[-1]
                            
                            try:
                                original_grantham_score, new_grantham_score = compute_grantham_score(ref_amino_acid, alt_amino_acid, original_alt_amino)
                                
                                if new_grantham_score >= original_grantham_score:
                                    
                                    criteria_df.loc[var_identifier, 'OM3'] = "yes"
                                    points_df.loc[var_identifier, 'OM3_p'] = 2
                                    break;
                                
                                else:
                                    
                                    criteria_df.loc[var_identifier, 'OM3'] = "no"
                                    points_df.loc[var_identifier, 'OM3_p'] = 0
                            except:
                                
                                criteria_df.loc[var_identifier, 'OM3'] = "no"
                                points_df.loc[var_identifier, 'OM3_p'] = 0
            
                else:
                    
                    criteria_df.loc[var_identifier, 'OM3'] = "no"
                    points_df.loc[var_identifier, 'OM3_p'] = 0
                
            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: MODERATE (2 points)
            ## Evidence: OM4
            
            ## Criteria: 
            ## Located in one of the hotspots in cancerhotspots with at < 50 samples 
            ## with a somatic variant at the same amino acid change count in cancerhotspots
            ## in at least 10.  
            
            # If the OM1 or OM3 criterion is applicable, then 
            if criteria_df.loc[var_identifier, 'OS3'] == "yes":
                
                criteria_df.loc[var_identifier, 'OM4'] = "no"
                points_df.loc[var_identifier, 'OM4_p'] = 0
                
            elif (criteria_df.loc[var_identifier, 'OM1'] == "yes") or (criteria_df.loc[var_identifier, 'OM3'] == "yes"):
                
                criteria_df.loc[var_identifier, 'OM4'] = "no"
                points_df.loc[var_identifier, 'OM4_p'] = 0
                    
            else:
            
                # For inframe deletion or insertion variants the in_frame_indel_dict is 
                # investigated 
                if any(s for s in indels if s in consequence):
                    
                    # Check if the gene is reported in cancerhotspots
                    if gene in in_frame_indel_dict.keys():
                        
                        gene_records_inframe_indel_dict = in_frame_indel_dict[gene]  
                        
                        for key in gene_records_inframe_indel_dict.keys():
                            
                            if not '-' in key:
                                
                                if amino_acid_position == key:
                                    
                                    res_variant_inframe_indel = [value for key, value in gene_records_inframe_indel_dict[key].items() if prot_change == key]
                                    
                                    if len(res_variant_inframe_indel) == 0:
                                            
                                        criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                        points_df.loc[var_identifier, 'OM4_p'] = 0
                                            
                                    else:
                                        
                                        # Check if the requirements from the guidelines are satisfied
                                        # by cancerhotspots records
                                        if (sum(list(map(int, gene_records_inframe_indel_dict[key].values()))) < 50) and (int(res_variant_inframe_indel[0]) >= 10):
                                    
                                            criteria_df.loc[var_identifier, 'OM4'] = 'yes'
                                            points_df.loc[var_identifier, 'OM4_p'] = 2
                                            break;
                                                
                                        else:
                                            
                                            criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                            points_df.loc[var_identifier, 'OM4_p'] = 0
                                
                                else:
                                    
                                    criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                    points_df.loc[var_identifier, 'OM4_p'] = 0
                                    
                            else:
                                
                                ## create range value 
                                min_range, max_range = key.split('-')                                    
                                
                                amino_acid_pos_list = re.findall("\d+", amino_acid_position)

                                if len(amino_acid_pos_list) >1:
                                    
                                    if (int(amino_acid_pos_list[0]) >= int(min_range)) and (int(amino_acid_pos_list[1]) <= int(max_range)):
                                        
                                        res_variant_inframe_indel = [value for key, value in gene_records_inframe_indel_dict[key].items() if prot_change == key]
                                        
                                        if len(res_variant_inframe_indel) == 0:
                                                
                                            criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                            points_df.loc[var_identifier, 'OM4_p'] = 0
                                                
                                        else:
                                                
                                           # Check if the requirements from the guidelines are satisfied
                                           # by cancerhotspots records
                                           if (sum(list(map(int, gene_records_inframe_indel_dict[key].values()))) < 50) and (int(res_variant_inframe_indel[0]) >= 10):
                                       
                                               criteria_df.loc[var_identifier, 'OM4'] = 'yes'
                                               points_df.loc[var_identifier, 'OM4_p'] = 2
                                               break;
                                                   
                                           else:
                                               
                                               criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                               points_df.loc[var_identifier, 'OM4_p'] = 0
                                                    
                                    else:
                                        criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                        points_df.loc[var_identifier, 'OM4_p'] = 0
                                        
                                                    
                                else:
                                    
                                    if int(min_range) <= int(amino_acid_pos_list[0]) <= int(max_range):
                                        
                                        res_variant_inframe_indel = [value for key, value in gene_records_inframe_indel_dict[key].items() if prot_change == key]
                                        
                                        if len(res_variant_inframe_indel) == 0:
                                                
                                            criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                            points_df.loc[var_identifier, 'OM4_p'] = 0
                                                
                                        else:
                                            
                                            if (sum(list(map(int, gene_records_inframe_indel_dict[key].values()))) < 50) and (int(res_variant_inframe_indel[0]) >= 10):
                                    
                                                criteria_df.loc[var_identifier, 'OM4'] = 'yes'
                                                points_df.loc[var_identifier, 'OM4_p'] = 2
                                                break;
                                            
                                            else:
                                                criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                                points_df.loc[var_identifier, 'OM4_p'] = 0
                                    
                                    else:
                                        criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                        points_df.loc[var_identifier, 'OM4_p'] = 0

                    else:
                        
                        criteria_df.loc[var_identifier, 'OM4'] = 'no'
                        points_df.loc[var_identifier, 'OM4_p'] = 0
                    
                else:
                    
                    # Otherwise, the single_residue_dict is considered.
                    # Check if the gene is reported in cancerhotspots
                    if gene in single_residue_dict.keys():
                        
                        res_var_single_residue = [value for key, value in single_residue_dict[gene].items() if ref_amino_acid_and_position == key]
                        
                        # Check if the variant evaluated is reported as residue in cancerhotspots
                        if len(res_var_single_residue) == 0: 
                            
                            criteria_df.loc[var_identifier, 'OM4'] = 'no'
                            points_df.loc[var_identifier, 'OM4_p'] = 0
                        
                        else:            
                            
                            res_var = [value for key, value in res_var_single_residue[0].items() if alt_amino_acid == key]
                            
                            if (len(res_var) != 0) and (sum(list(map(int, res_var_single_residue[0].values() ))) < 50) and (int(res_var[0]) >= 10):
                                    
                                criteria_df.loc[var_identifier, 'OM4'] = 'yes'
                                points_df.loc[var_identifier, 'OM4_p'] = 2
                                
                            else:
                                criteria_df.loc[var_identifier, 'OM4'] = 'no'
                                points_df.loc[var_identifier, 'OM4_p'] = 0 
                                    
                    else:
                        
                        criteria_df.loc[var_identifier, 'OM4'] = 'no'
                        points_df.loc[var_identifier, 'OM4_p'] = 0
                
            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: SUPPORTING (1 points)
            ## Evidence: OP1
            
            ## Criteria: 
            ## All used lines of computational evidence support an oncogenic effect of a 
            ## variant (conservation/evolutionary, splicing effect, etc.).
            
            if (spliceAI_score == "") and (phylop_score == "") and (phastCons_score == ""):
                
                criteria_df.loc[var_identifier, 'OP1'] = 'no'
                points_df.loc[var_identifier, 'OP1_p'] = 0
            
            elif (spliceAI_score != "" and spliceAI_score == "PASS") or (phylop_score != "" and phylop_score >= 0.5) or (phastCons_score != "" and phastCons_score >= 0.5):
                    
                    criteria_df.loc[var_identifier, 'OP1'] = 'yes'
                    points_df.loc[var_identifier, 'OP1_p'] = 1

            else:
                
                criteria_df.loc[var_identifier, 'OP1'] = 'no'
                points_df.loc[var_identifier, 'OP1_p'] = 0

            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: SUPPORTING (1 points)
            ## Evidence: OP3
            
            ## Criteria: 
            ## Located in one of the hotspots in cancerhotspots and in particular amino
            ## acid change count in cancerhotspots is below 10.  
            
            # For inframe deletion or insertion variants the in_frame_indel_dict is 
            # investigated 
            if criteria_df.loc[var_identifier, 'OS3'] == "yes": 
                
                criteria_df.loc[var_identifier, 'OP3'] = 'no'
                points_df.loc[var_identifier, 'OP3_p'] = 0
                
            elif criteria_df.loc[var_identifier, 'OM4'] == "yes":
                
                criteria_df.loc[var_identifier, 'OP3'] = 'no'
                points_df.loc[var_identifier, 'OP3_p'] = 0
                
            else:
                if any(s for s in indels if s in consequence):
                    
                    # Check if the gene is reported in cancerhotspots
                    if gene in in_frame_indel_dict.keys():
                        
                        gene_records_inframe_indel_dict = in_frame_indel_dict[gene]  
                        
                        for key in gene_records_inframe_indel_dict.keys():
                            
                            if not '-' in key:
                                
                                if amino_acid_position == key:
                                    
                                    res_variant_inframe_indel = [value for key, value in gene_records_inframe_indel_dict[key].items() if prot_change == key]
                                    
                                    if len(res_variant_inframe_indel) == 0:
                                            
                                        criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                        points_df.loc[var_identifier, 'OP3_p'] = 0
                                            
                                    else:
                                        
                                        # Check if the requirements from the guidelines are satisfied
                                        # by cancerhotspots records
                                        if int(res_variant_inframe_indel[0]) < 10:
              
                                            criteria_df.loc[var_identifier, 'OP3'] = 'yes'
                                            points_df.loc[var_identifier, 'OP3_p'] = 1
                                            break;
                                                
                                        else:
                                            
                                            criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                            points_df.loc[var_identifier, 'OP3_p'] = 0
                                
                                else:
                                    
                                    criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                    points_df.loc[var_identifier, 'OP3_p'] = 0
                                    
                            else:
                                
                                ## create range value 
                                min_range, max_range = key.split('-')                                    
                                
                                amino_acid_pos_list = re.findall("\d+", amino_acid_position)

                                if len(amino_acid_pos_list) >1:
                                    
                                    if (int(amino_acid_pos_list[0]) >= int(min_range)) and (int(amino_acid_pos_list[1]) <= int(max_range)):
                                        
                                        res_variant_inframe_indel = [value for key, value in gene_records_inframe_indel_dict[key].items() if prot_change == key]
                                        
                                        if len(res_variant_inframe_indel) == 0:
                                                
                                            criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                            points_df.loc[var_identifier, 'OP3_p'] = 0
                                                
                                        else:
                                                
                                           # Check if the requirements from the guidelines are satisfied
                                           # by cancerhotspots records
                                           if int(res_variant_inframe_indel[0]) < 10:
                                               
                                               criteria_df.loc[var_identifier, 'OP3'] = 'yes'
                                               points_df.loc[var_identifier, 'OP3_p'] = 1
                                               break;
                                                   
                                           else:
                                               
                                               criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                               points_df.loc[var_identifier, 'OP3_p'] = 0
                                                    
                                    else:
                                        criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                        points_df.loc[var_identifier, 'OP3_p'] = 0
                                        
                                                    
                                else:
                                    
                                    if int(min_range) <= int(amino_acid_pos_list[0]) <= int(max_range):
                                        
                                        res_variant_inframe_indel = [value for key, value in gene_records_inframe_indel_dict[key].items() if prot_change == key]
                                        
                                        if len(res_variant_inframe_indel) == 0:
                                                
                                            criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                            points_df.loc[var_identifier, 'OP3_p'] = 0
                                                
                                        else:
                                            
                                            if int(res_variant_inframe_indel[0]) < 10:
                                               
                                                criteria_df.loc[var_identifier, 'OP3'] = 'yes'
                                                points_df.loc[var_identifier, 'OP3_p'] = 1
                                                break;
                                            
                                            else:
                                                criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                                points_df.loc[var_identifier, 'OP3_p'] = 0
                                    
                                    else:
                                        criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                        points_df.loc[var_identifier, 'OP3_p'] = 0

                    else:
                        
                        criteria_df.loc[var_identifier, 'OP3'] = 'no'
                        points_df.loc[var_identifier, 'OP3_p'] = 0
                    
                else:
                    
                    # Otherwise, the single_residue_dict is considered.
                    # Check if the gene is reported in cancerhotspots
                    if gene in single_residue_dict.keys():
                        
                        res_var_single_residue = [value for key, value in single_residue_dict[gene].items() if ref_amino_acid_and_position == key]
                        
                        # Check if the variant evaluated is reported as residue in cancerhotspots
                        if len(res_var_single_residue) == 0: 
                            
                            criteria_df.loc[var_identifier, 'OP3'] = 'no'
                            points_df.loc[var_identifier, 'OP3_p'] = 0
                        
                        else:            
                            
                            res_var = [value for key, value in res_var_single_residue[0].items() if alt_amino_acid == key]
                            
                            if (len(res_var) != 0) and (int(res_var[0]) < 10):
                           
                                criteria_df.loc[var_identifier, 'OP3'] = 'yes'
                                points_df.loc[var_identifier, 'OP3_p'] = 1
                                
                            else:
                                criteria_df.loc[var_identifier, 'OP3'] = 'no'
                                points_df.loc[var_identifier, 'OP3_p'] = 0 
                                    
                    else:
                        
                        criteria_df.loc[var_identifier, 'OP3'] = 'no'
                        points_df.loc[var_identifier, 'OP3_p'] = 0

            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: SUPPORTING (1 points)
            ## Evidence: OP4
            
            ## Criteria: 
            ## Absent from controls (or at an extremely low frequency) in gnomAD. 
            
            if (gnomADe_freq == "") and (gnomADg_freq == ""): 

                criteria_df.loc[var_identifier, 'OP4'] = 'yes'
                points_df.loc[var_identifier, 'OP4_p'] = 1
            
            elif (gnomADg_freq != "") and (float(gnomADg_freq) <= 0.01): 

                criteria_df.loc[var_identifier, 'OP4'] = 'yes'
                points_df.loc[var_identifier, 'OP4_p'] = 1
               
            elif (gnomADe_freq != "") and (float(gnomADe_freq) <= 0.01):

                criteria_df.loc[var_identifier, 'OP4'] = 'yes'
                points_df.loc[var_identifier, 'OP4_p'] = 1
                
            else:

                criteria_df.loc[var_identifier, 'OP4'] = 'no'
                points_df.loc[var_identifier, 'OP4_p'] = 0
            
         
            ###############################################################################
            # -------------------------------------------------------------------------- ##
            
            ## CRITERIA FOR EVIDENCE OF BENIGN EFFECT OF SOMATIC VARIANTS
            
            ## Category: VERY STRONG (-8 points)
            ## Evidence: SBVS1
            
            ## Criteria: 
            ## Minor allele frequency is >5% in gnomAD in any 5 general continental 
            ## populations: African, East Asian, European (non-Finnish), Latino, and South
            ## Asian.
            
            # exclude empty strings from gnomADg_freq_pop and gnomADe_freq_pop
            real_gnomADg_freq_pop = list(filter(None, gnomADg_freq_pop))
            real_gnomADe_freq_pop = list(filter(None, gnomADe_freq_pop))
            
            # exclude zeros from gnomADg_freq_pop and gnomADe_freq_pop
            real_gnomADg_freq_pop = [i for i in real_gnomADg_freq_pop if i != '0']
            real_gnomADe_freq_pop = [i for i in real_gnomADe_freq_pop if i != '0']
            
            # If we don´t have MAF in 5 general continental populations
            # from exome and genome studies the criteria is not triggered. 
            if (real_gnomADg_freq_pop == []) and (real_gnomADe_freq_pop == []):
                
                criteria_df.loc[var_identifier, 'SBVS1'] = 'no'
                points_df.loc[var_identifier, 'SBVS1_p'] = 0
                
            elif any(float(item) > 0.05 for item in real_gnomADg_freq_pop):   
                 
                criteria_df.loc[var_identifier, 'SBVS1'] = 'yes'
                points_df.loc[var_identifier, 'SBVS1_p'] = -8
                   
            elif any(float(item) > 0.05 for item in real_gnomADe_freq_pop):

                criteria_df.loc[var_identifier, 'SBVS1'] = 'yes'
                points_df.loc[var_identifier, 'SBVS1_p'] = -8
                    
            else: 

                criteria_df.loc[var_identifier, 'SBVS1'] = 'no'
                points_df.loc[var_identifier, 'SBVS1_p'] = 0


            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: STRONG (-4 points)
            ## Evidence: SBS1
            
            ## Criteria: 
            ## Minor allele frequency is >1% in gnomAD in any 5 general continental 
            ## populations: African, East Asian, European (non-Finnish), Latino, and South
            ## Asian.
            
            # If minor allele frequency is absent in any 5 general continental pop
            # the criteria is not triggered. 
            if (real_gnomADg_freq_pop == []) and (real_gnomADe_freq_pop == []):
                
                criteria_df.loc[var_identifier, 'SBS1'] = 'no'
                points_df.loc[var_identifier, 'SBS1_p'] = 0
                
            elif any(float(item) > 0.01 for item in real_gnomADg_freq_pop):   
                 
                criteria_df.loc[var_identifier, 'SBS1'] = 'yes'
                points_df.loc[var_identifier, 'SBS1_p'] = -4
                   
            elif any(float(item) > 0.01 for item in real_gnomADe_freq_pop):

                criteria_df.loc[var_identifier, 'SBS1'] = 'yes'
                points_df.loc[var_identifier, 'SBS1_p'] = -4
                    
            else: 

                criteria_df.loc[var_identifier, 'SBS1'] = 'no'
                points_df.loc[var_identifier, 'SBS1_p'] = 0
                      

            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: STRONG (-4 points)
            ## Evidence: SBS2
            
            ## Criteria: 
            ## Well-established in vitro or in vivo functional studies show no oncogenic 
            ## effects. 
            
            if ("Benign" in clinvar_classification) and (review_status in review_level_list):
                
                criteria_df.loc[var_identifier, 'SBS2'] = 'yes'
                points_df.loc[var_identifier, 'SBS2_p'] = -4
                
            elif ("Likely benign" in clinvar_classification) and (review_status in review_level_list):
                
                criteria_df.loc[var_identifier, 'SBS2'] = 'yes'
                points_df.loc[var_identifier, 'SBS2_p'] = -4
                
            else:
                
                criteria_df.loc[var_identifier, 'SBS2'] = 'no'
                points_df.loc[var_identifier, 'SBS2_p'] = 0
                
         
            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: SUPPORTING (-1 points)
            ## Evidence: SBP1
            
            ## Criteria: 
            ## All used lines of computational evidence suggest no effect of a variant
            ## (conservation/evolutionary, splicing effect, etc.).
            
            if (spliceAI_score == "") and (phylop_score == "") and (phastCons_score == ""):
                
                criteria_df.loc[var_identifier, 'SBP1'] = 'no'
                points_df.loc[var_identifier, 'SBP1_p'] = 0
            
            elif (spliceAI_score != "" and spliceAI_score == "FAIL") and (phylop_score != "" and phylop_score < 0.5) and (phastCons_score != "" and phastCons_score < 0.5):
                    
                criteria_df.loc[var_identifier, 'SBP1'] = 'yes'
                points_df.loc[var_identifier, 'SBP1_p'] = -1
                    
            else:
                criteria_df.loc[var_identifier, 'SBP1'] = 'no'
                points_df.loc[var_identifier, 'SBP1_p'] = 0     
         
            #-----------------------------------------------------------------------------#
            # -------------------------------------------------------------------------- ##
            
            ## Category: SUPPORTING (-1 points)
            ## Evidence: SBP2
            
            ## Criteria: 
            ## A synonymous (silent) variant for which splicing prediction algorithms 
            ## predict no effect on the splice consensus sequence nor the creation of a new
            ## splice site AND the nucleotide is not highly conserved. 
             
            if "synonymous_variant" in consequence:
                
                if (spliceAI_score == "") and (phylop_score == "") and (phastCons_score == ""):
                    
                    criteria_df.loc[var_identifier, 'SBP2'] = 'no'
                    points_df.loc[var_identifier, 'SBP2_p'] = 0
                
                elif (spliceAI_score != "" and spliceAI_score == "FAIL") and (phylop_score != "" and phylop_score < 0.5) and (phastCons_score != "" and phastCons_score < 0.5):
                    
                   criteria_df.loc[var_identifier, 'SBP2'] = 'yes'
                   points_df.loc[var_identifier, 'SBP2_p'] = -1
                    
                else:
                   criteria_df.loc[var_identifier, 'SBP2'] = 'no'
                   points_df.loc[var_identifier, 'SBP2_p'] = 0
                   
            else:
                criteria_df.loc[var_identifier, 'SBP2'] = 'no'
                points_df.loc[var_identifier, 'SBP2_p'] = 0
                
            
            ###############################################################################
            # -------------------------------------------------------------------------- ##
            
            # Sum of the points
            sum_points = points_df.loc[var_identifier].sum()
            
            # Qualitative classification
            if (sum_points <= -7):
                
                classification = "Benign"
                
            elif (sum_points >= -6) and (sum_points <= -1):
                
                classification = "Likely Benign"
                
            elif (sum_points >= 0) and (sum_points <= 5):
                
                classification = "VUS"
                
            elif (sum_points >= 6) and (sum_points <= 9):
                
                classification = "Likely Oncogenic"
                
            elif (sum_points >= 10):
                
                classification = "Oncogenic"
                
            if var_identifier not in classification_dict.keys():
                
                classification_dict[var_identifier] = {}
                classification_dict[var_identifier] = {'Points' : sum_points,
                                                       'Classification' : classification}
                
            else:
                
                classification_dict[var_identifier] = {'Points' : sum_points,
                                                       'Classification' : classification}
            
            # -------------------------------------------------------------------------- ##
            
            # triggered criteria
        
            triggered_criteria_list = criteria_df.loc[var_identifier][criteria_df.loc[var_identifier] == "yes"].index.format()
            
            # triggered criteria stored in a data frame
            triggered_criteria_df.loc[var_identifier, 'Criteria'] = ','.join(triggered_criteria_list)
            
        # -------------------------------------------------------------------------- ##
        
        # Data frame wit the classification (quantitative and qualitative) of the assessed variants 
        classification_df = pd.DataFrame.from_dict(classification_dict,  orient='index', columns = ["Points", "Classification"])
            
        # Data frame wit the classification, triggered criteria, points
        classification_df_final = pd.concat([classification_df, triggered_criteria_df, criteria_df, points_df], axis=1)
        
        # Save to file
        output_data_dir = os.path.join(output_dir, str(name) + "_OncoVI_eval.csv")
        classification_df_final.to_csv(output_data_dir, sep = ',', index = True)
