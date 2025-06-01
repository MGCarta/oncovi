# Oncogenicity Variant Interpreter (OncoVI)
OncoVI is a fully-automated Python implementation of the [oncogenicity guidelines](https://pubmed.ncbi.nlm.nih.gov/35101336/) by Horak et al. (Genetics in Medicine, 2022). 

Starting from the genomic location of the variants, OncoVI:
1. performs functional annotation based on the [Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) from Ensembl;
2. collects biological evidences from the implemented publicly available resources;
3. classifies the oncogenicity of somatic variants, based on the point-based system for combining pieces of evidence defined by Horak et al.

More detailed information on the resources used by OncoVI, the implementation of the oncogenicity guidelines, and its application to real-world data can be found in [our pre-print](https://www.medrxiv.org/content/10.1101/2024.10.10.24315072v1).

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

### ClinVar resources
The download and preparation of the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) resource utilised by the functional annotation STEP is handled by the script ```01_clinvar_resource_manager.sh```.
```rb
# Move to the folder where the script is located
cd oncovi/src/
```
```rb
# Run the bash script
bash 01_clinvar_resource_manager.sh
```

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
The spliceAI plugin is used during the the functional annotation STEP. Detailed information on how to set up the spliceAI plugin for VEP can be found [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#spliceAI). The spliceAI Plugin must be enabled in the script ```vep.sh``` according to the Plugin instructions.  

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

## References
Please cite our preprint [Oncogenicity Variant Interpreter (OncoVI): oncogenicity guidelines implementation to support somatic variants interpretation in precision oncology](https://www.medrxiv.org/content/10.1101/2024.10.10.24315072v1) if you decide to use OncoVI.
