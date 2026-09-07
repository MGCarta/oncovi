# Oncogenicity Variant Interpreter (OncoVI)
OncoVI is a fully-automated Python implementation of the [oncogenicity guidelines](https://pubmed.ncbi.nlm.nih.gov/35101336/) by Horak et al. (Genetics in Medicine, 2022). 

Starting from the genomic location of the variants, OncoVI:
1. performs functional annotation based on the [Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) from Ensembl;
2. collects biological evidences from the implemented publicly available resources;
3. classifies the oncogenicity of somatic variants, based on the point-based system for combining pieces of evidence defined by Horak et al.

More detailed information on the resources used by OncoVI, the implementation of the oncogenicity guidelines, and its application to real-world data can be found in [our paper](https://doi.org/10.1016/j.jmoldx.2026.03.004).

Workflow of OncoVI: 
![alt text][logo]

[logo]: https://github.com/MGCarta/oncovi/blob/main/figures/oncovi_v2.JPG "Logo Title Text 2"

###### The figure shows the implemented criteria in OncoVI (11 and five criteria for evidence of oncogenic and benign effect respectively), the public resources utilised to assess each criterion, the points associated with each criterion and the classification of oncogenicity into one of five classes on the basis of the variant-specific score, obtained as the sum of the points associated to the criteria triggered by OncoVI for the variant: score≥10:Oncogenic (O), 6≤score≤9:Likely Oncogenic (LO), 0≤score≤5:Variant of uncertain significance (VUS), -6≤score≤-1:Likely Benign (LB), score≤-7:Benign (B). Blue: resources suggested by the Standard Operating Procedure by Horak et al.; black: resources identified by the authors of OncoVI.

## Software requirements
OncoVI was implemented and tested on a dedicated conda enviroment running on a remote server based on Ubuntu 20.04.4 long-term support (LTS) operating system. To use this repository make sure to have:

* conda >= 24.11.1
* python >=3.8.8

OncoVI has been tested on:
* Ensembl VEP conda package >=v.111
* VEP Plugins (dbNSFP and spliceAI)

OncoVI supports only variants in human genome assembly GRCh38. No other genome assemblies are currently supported.

## Get started
Clone the GitHub repository:
```rb
git clone https://github.com/MGCarta/oncovi.git
```

### Create the conda environment and install the required packages
1. Create the conda environment envpy310
```rb
conda create --name envpy310 python=3.10
```
2. Activate the conda environment
```rb
conda activate envpy310
```
3. Install pandas
```rb
pip install pandas
```
4. Install colorama
```rb
pip install colorama
```
5. Install openpyxl
```rb
pip install openpyxl
```

## Resources
OncoVI requires a set of pre-processed resources for the oncogenicity classification.

Most resources were generated from the corresponding databases using the preparation script:
```resources/01_prepare_resources.py``` 

and are distributed as ready-to-use resources in the ```resources/``` directory of this repository. The original input files are not included in this repository.

In case the user wants to collect the most updated resources, the ```01_prepare_resources.py``` can be used.
1. Create a folder with the original resources to process:
    ```text
    01_original_resources/
    ├── Census_....csv
    ├── cancerGeneList.tsv
    ├── MutSpliceDB_....csv
    ├── amino_conversion.txt
    ├── ...
    ```
2. Run the script:
```rb
python 01_prepare_resources.py \
    --input-dir /path/to/original_resources \
    --output-dir /path/to/resources
```

The resources currently provided include:

| Resource  | Source |
| ------------- | ------------- |
| mut_splice_dict.txt  | [MutSpliceDB](https://brb.nci.nih.gov/splicing/index.html)  |
| amino_dict.txt  | Amino-acid conversion table  |
| single_residue_dict.txt  | [Cancer Hotspots](https://www.cancerhotspots.org/#/home)  |
| inframe_indel_dict.txt  | [Cancer Hotspots](https://www.cancerhotspots.org/#/home)  |
| cgi_dictionary.txt  | [Cancer Genome Interpreter](https://www.cancergenomeinterpreter.org/)  |
| domains_dictionary.txt  | [UniProt](https://www.uniprot.org/)  |
| cosmic_hgvsg_dictionary.txt.gz  | [COSMIC](https://cancer.sanger.ac.uk/cosmic)  |
| cosmic_all_dictionary.txt.gz  | [COSMIC](https://cancer.sanger.ac.uk/cosmic)  |
| oncogenes_cgc.csv  | [COSMIC Cancer Gene Census](https://cancer.sanger.ac.uk/cmc/home) |
| tsg_cgc.csv  | [COSMIC Cancer Gene Census](https://cancer.sanger.ac.uk/cmc/home) |
| tsg_tier1.csv  | [COSMIC Cancer Gene Census](https://cancer.sanger.ac.uk/cmc/home) |
| og_oncokb.csv  | [OncoKB](https://www.oncokb.org/cancer-genes)  |
| tsg_oncokb.csv  | [OncoKB](https://www.oncokb.org/cancer-genes)  |
| bona_fide_tsg.txt  | [COSMIC](https://cancer.sanger.ac.uk/cosmic) + [OncoKB](https://www.oncokb.org/cancer-genes)  |
| ogs_list.txt  | [COSMIC](https://cancer.sanger.ac.uk/cosmic) + [OncoKB](https://www.oncokb.org/cancer-genes)  |
| tsg_list.txt  | [COSMIC](https://cancer.sanger.ac.uk/cosmic) + [OncoKB](https://www.oncokb.org/cancer-genes)  |

The corresponding source/version and access information for each resource are documented in ```resources/README.txt```.

### ClinVar resources
ClinVar is handled separately because its database can be downloaded and prepared locally.

The download and preparation of the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) resource utilised by the functional annotation STEP are automated by the script ```01_clinvar_resource_manager.sh```.

The script downloads the ClinVar ```variant_summary.txt.gz``` file, extracts the GRCh38 variants, reduces the dataset to the fields required by OncoVI, and generates the ClinVar dictionary using:
```src/create_clinvar_dict.py```

To prepare the ClinVar resource:

```rb
# Move to the folder where the script is located
cd oncovi/src/
```
```rb
# Run the bash script
bash 01_clinvar_resource_manager.sh
```

The ClinVar resource is therefore not included as a ready-to-use file in the repository.

### Instructions to set up VEP manually
> [!NOTE]
Make sure to have enough space in your current directory. Cache files require almost 30 GB (~26 GB for ```homo_sapiens_refseq_vep_114_GRCh38.tar.gz``` and ~1 GB for ```homo_sapiens_vep_114_GRCh38.tar.gz```).

1. Install the conda package Ensembl Vep 
```rb
conda install bioconda::ensembl-vep
```
2. Run the installer (the latest version of the package will be automatically installed)
```rb
vep_install --NO_HTSLIB -a fp -s homo_sapiens -y GRCh38 --PLUGINS all -c '../.vep' -r '../.vep/Plugins/'
```
3. Install the cache files
```rb
# Move to the folder generated to save cache files
cd ../.vep
```
```rb
# Download the cache files
wget ftp://ftp.ensembl.org/pub/release-114/variation/indexed_vep_cache/homo_sapiens_refseq_vep_114_GRCh38.tar.gz
```
```rb
# Extract the archive
tar -xzf homo_sapiens_refseq_vep_114_GRCh38.tar.gz
```
```rb
# Remove the archive
rm homo_sapiens_refseq_vep_114_GRCh38.tar.gz
```

### dbNSFP plugin
The dbNSFP plugin is used by the the functional annotation STEP. VEP reports detailed information on how to set up the dbNSFP plugin [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp).

About dbNSFP data files:

* Download dbNSFP files from https://sites.google.com/site/jpopgen/dbNSFP.
* There are two distinct branches of the files provided for academic and  commercial usage. Please use the appropriate files for your use case.
* The file must be processed depending on dbNSFP release version and assembly  (see commands below).
* The resulting file must be indexed with tabix before use by this plugin  (see commands below).

For release 4.9c:
```rb
version=4.9c
wget https://dbnsfp.s3.amazonaws.com/dbNSFP${version}.zip
unzip dbNSFP${version}.zip
zcat dbNSFP${version}_variant.chr1.gz | head -n1 > h
```
GRCh38/hg38 data
```rb
zgrep -h -v ^#chr dbNSFP${version}_variant.chr* | sort -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP${version}_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP${version}_grch38.gz
```
Replace in the script ```vep.sh``` of this repository the full path to the dbNSFP resource:
```rb
--plugin dbNSFP,/path/to/dbNSFP4.9c_grch38.gz
```

### spliceAI plugin
The spliceAI plugin is used during the the functional annotation STEP.

To download the data:

1. Go to https://basespace.illumina.com/s/otSPW8hnhaZR
2. Login with an Illumina account (a free account can be created if needed).
3. Accept the invitation to get access to the project Predicting splicing from primary sequence
4. Make sure that the Context to accept invite is set to Personal and click Accept. The project data will now be accessible in your account.
5. On the navigation bar at the top, click Projects and then click on Predicting splicing from primary sequence.
6. Below the project title, click Other datasets.
7. Click genome_scores_v1.3 in the list of files.
8. On the new page, click genome_scores_v1.3 again, which will expand to show the folder contents.
9. Click on the files (from those with raw in the filename) to download from the list: ```spliceai_scores.raw.indel.hg38.vcf.gz```, ```spliceai_scores.raw.indel.hg38.vcf.gz.tbi```, ```spliceai_scores.raw.snv.hg38.vcf.gz```, ```spliceai_scores.raw.snv.hg38.vcf.gz.tbi```.

Replace in the script ```vep.sh``` of this repository the full path to the spliceAI resources:
```rb
--plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/path/to/spliceaiceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5
```

## Test OncoVI on an exemplary set of variants

#### Input data
Variants in both text format and in variant call format (VCF) are accepted by OncoVI. No other formats are currently supported.

An exemplary input data in text format is available under:
```rb
# ~/oncovi/testdata/SOP_table_union.txt
```
An exemplary input data in VCF format is available under:
```rb
# ~/oncovi/testdata/onco_som_GRCh38_ClinVar_vep_input.vcf
```

#### Run OncoVI

```rb
# Move to the directory where the python script 02_VEP_based_pipeline.py is located
cd ../src
```
```rb
# Run OncoVI
python 02_VEP_based_pipeline.py -i ../testdata/SOP_table_union.txt
```

#### OncoVI output
The current output format produced by OncoVI is a tab-delimited excel file. It includes all the variants provided in input for which the functional annotation STEP and the OncoVI prediction STEP were successful. 
An exemplary output file is available here for inspection:
```rb
# ~/oncovi/exampleoutput/SOP_table_union_OncoVI_prediction.xlsx
```

The output file includes:
* Features characterising the variants that were collected during either the VEP annotation step (i.e., ```SYMBOL```,	```HGVSc```,	```HGVSp```,	```MANE_SELECT```,	```EXON```,	```CHROM```,	```POS```,	```REF```,	```ALT```,	```Amino_acids```,	```Existing_variation```,	```Consequence```,	```Transcript #1```,	```gnomADe_AF```,	```gnomADg_AF```) or the ClinVar interrogation step (i.e., ```ClinVar_germline```,	```ClinVar_germline_ReviewStatus```,	```ClinVar_Oncogenicity```,	```ClinVar_Oncogenicity_ReviewStatus```,	```ClinVar_Somatic_Clinical_Impact```,	```ClinVar_Somatic_ReviewStatus```,	```Clinvar_url```)
  
* Features calculated by OncoVI: the variant-specific score, i.e., the sum of the points associated with the criteria triggered by OncoVI (```Points```), the oncogenicity classification associated with the variant-specific score according to the point-based system provided by Horak et al. (```Classification```), the criteria triggered by OncoVI based on the evidences collected	(```Criteria```)

## OncoVI issues
Please, help us to improve OncoVI by describing your bug/issue in detail


## License
The MIT license file applies to only the scripts within this repository.

OncoVI makes use of resources that require a dedicated license agreement, such as [OncoKB](https://www.oncokb.org/terms) and [COSMIC](https://www.cosmickb.org/licensing/). 
For these reasons, OncoVI is intended for research purposes only and its use outside of this context is under the responsibility of the user. 
Please, visit the relative websites and verify that you are part of an academic institution to freely use OncoVI. 
It is the user's responsibility to carefully check and comply with the licenses of the resources that need to be additionally installed to use OncoVI.

## Using OncoVI without OncoKB
OncoVI can be used without the OncoKB resource. In this configuration, the role of tumour suppressor genes and oncogenes is defined using COSMIC only and the other resources available to OncoVI.

We evaluated this configuration using the 7,802 somatic variants from the Molecular Tumour Board (MTB) data set, comprising variants detected in more than 500 tumours analysed with the Illumina TruSight Oncology 500 (TSO500) gene panel.
We compared the results obtained with OncoVI version 1 (which uses OncoKB). and version 2, which does not use OncoKB. 

Overall, the two configurations showed high agreement, with an accuracy of 0.984. The oncogenicity classification changed for 126 of 7,802 variants, while the variant-specific score changed for 162 variants. For the variants with a changed score, the score was consistently lower in the configuration without OncoKB and was associated with the loss of either the OVS1 or OM2 criterion. 

The complete analysis and results are available in the following [report](https://github.com/MGCarta/oncovi/blob/main/docs/OncoVI_without_oncoKB_MTB_public.pdf).

The resources and scripts required to run OncoVI without OncoKB are provided in this repository. This would allow users to run OncoVI without requiring access to the licensed OncoKB resource. 

## References
Please cite our paper [Oncogenicity Variant Interpreter (OncoVI): oncogenicity guidelines implementation to support somatic variants interpretation in precision oncology](https://doi.org/10.1016/j.jmoldx.2026.03.004) if you decide to use OncoVI.
