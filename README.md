# Oncogenicity Variant Interpreter (OncoVI)
OncoVI is a fully-automated Python implementation of the [oncogenicity guidelines](https://pubmed.ncbi.nlm.nih.gov/35101336/) by Horak et al. (Genetics in Medicine, 2022). 

Starting from the genomic location of the variants in human genome assembly GRCh38, OncoVI:
1. performs functional annotation based on the [Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) from Ensembl;
2. collects biological evidences from the implemented publicly available resources;
3. classifies the oncogenicity of somatic variants, based on the point-based system for combining pieces of evidence defined by Horak et al.

More detailed information on the resources used by OncoVI, the implementation of the oncogenicity guidelines, and its application to real-world data can be found in [our pre-print](https://www.medrxiv.org/content/10.1101/2024.10.10.24315072v1).

Workflow of OncoVI: 
![alt text][logo]

[logo]: https://github.com/MGCarta/oncovi/blob/main/figures/Figure1_oncovi.PNG "Logo Title Text 2"

###### The figure shows the implemented criteria in OncoVI (11 and five criteria for evidence of oncogenic and benign effect respectively), the public resources utilised to assess each criterion, the points associated with each criterion and the classification of oncogenicity into one of five classes on the basis of the variant-specific score, obtained as the sum of the points associated to the criteria triggered by OncoVI for the variant: score≥10:Oncogenic (O), 6≤score≤9:Likely Oncogenic (LO), 0≤score≤5:Variant of uncertain significance (VUS), -6≤score≤-1:Likely Benign (LB), score≤-7:Benign (B). Blue: resources suggested by the Standard Operating Procedure by Horak et al.; black: resources identified by the authors of OncoVI.

## Software requirements
OncoVI was implemented and tested on a dedicated conda enviroment running on a remote server based on Ubuntu 20.04.4 long-term support (LTS) operating system. To use this repository make sure to have:

* conda >= 24.11.1
* python >=3.8.8
* python site-packages (numpy, pandas, subprocess)

OncoVI has been tested on:
* Ensembl VEP (v. 111)
* VEP Plugins (dbNSFP and spliceAI)

## Create the conda environment and install the required package
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

## Get started
Clone the GitHub repository:
```rb
git clone https://github.com/MGCarta/oncovi.git
```
### ClinVar resources
The download and preparation of the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) resource utilised by the functional annotation STEP is handled by the script ```01_clinvar_resource_manager.sh```.
```rb
# Move to the folder where the script is located
cd oncovi/src/
# Run the bash script
bash 01_clinvar_resource_manager.sh
# Deactivate the envpy310 environment
conda deactivate
```
Create the conda environment:
```rb
# Create the conda environment for oncovi
conda env create -n oncovi -f /path/to/OncoVIenvFile.yml
```
```rb
# Activate the conda environment
conda activate oncovi
```
## Set up VEP
```rb
# Run the installer (v. 111) available in the conda environment
vep_install --NO_HTSLIB -c '/path/to/.vep' -r '/path/to/.vep/Plugins/'
```
Then:
  * select homo_sapiens_refseq_vep_111_GRCh38.tar.gz as cache
  * homo_sapiens_vep_111_GRCh38.tar.gz as reference genome
  * install all Plugins

### dbNSFP plugin
The dbNSFP plugin is used by the the functional annotation STEP. Detailed information on how to set up the dbNSFP plugin for VEP can be found [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp). The dbNSFP Plugin must be enabled in the script ```vep.sh``` according to the Plugin instructions.

### spliceAI plugin
The spliceAI plugin is used during the the functional annotation STEP. Detailed information on how to set up the spliceAI plugin for VEP can be found [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#spliceAI). The spliceAI Plugin must be enabled in the script ```vep.sh``` according to the Plugin instructions.  

## Prepare your variants
Both variants in text format and in variant call format (VCF) are accepted by [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html). Please refer to VEP [official documentation](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#input) for a detailed description of input formats.
A test data is available under:

```rb
# /oncovi/testdata/SOP_table_union.txt
```

#### Perform functional annotation via VEP

```rb
# Navigate to the directory in which the python script 02_VEP_based_pipeline.py is located

# Run the functional annotation
python 02_VEP_based_pipeline.py -i /path/to/oncovi/testdata/SOP_table_union.txt
```

#### Run OncoVI

```rb
# Navigate to the directory in which the python script 03_OncoVI_SOP.py is located

# Run OncoVI
python 03_OncoVI_SOP.py
```

## OncoVI issues
Please, help us to improve OncoVI by describing your bug/issue in detail


## License
The MIT license file applies to only the scripts within this repository.

OncoVI makes use of resources that require a dedicated license agreement, such as [OncoKB](https://www.oncokb.org/terms) and [COSMIC](https://www.cosmickb.org/licensing/). 
For these reasons, OncoVI is intended for research purposes only and its use outside of this context is under the responsibility of the user. 
Please, visit the relative websites and verify that you are part of an academic institution to freely use OncoVI. 
It is the user´s responsibility to carefully check and comply with the licenses of the resources that need to be additionally installed to use OncoVI.

## References
Please cite our preprint [Oncogenicity Variant Interpreter (OncoVI): oncogenicity guidelines implementation to support somatic variants interpretation in precision oncology](https://www.medrxiv.org/content/10.1101/2024.10.10.24315072v1) if you decide to use OncoVI.
