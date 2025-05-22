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
OncoVI was implemented and tested on a dedicated conda enviroment running on a remote server based on Ubuntu 20.04.4 long-term support (LTS) operating system. To run OncoVI the following packages are required:

* python and site-packages (numpy, pandas, subprocess)
* Ensembl VEP and VEP plugins dbNSFP and spliceAI
* external resources, such as COSMIC


### COSMIC resources
Due to size constraints, the [COSMIC](https://cancer.sanger.ac.uk/cosmic/download/cosmic) resources utilised by OncoVI could not be uploaded on the GitHub repo. How to download and handle COSMIC data to make them usable by OncoVI is described here below.

#### Cancer Mutation Census 

1. First, [All Data CMC](https://cancer.sanger.ac.uk/cosmic/download/cancer-mutation-census/v100/alldata-cmc) for genome GRCh38 was downloaded
2. Then, the data set was reduced to the columns: ```GENE_NAME```, ```Mutation CDS```, ```Mutation AA```, ```AA_MUT_START```, ```Mutation genome position GRCh38```, and ```COSMIC_SAMPLE_MUTATED```
3. The reduced data set was converted into a dictionary with the python script ``` /src/prepare_cosmic_resources.py ```
4. The resulted dictionary was saved under the name ```cosmic_all_dictionary.txt```
5. The path to the dictionary must be provided to the python script ```03_OncoVI_SOP.py```

#### Census Genes Mutations 

1. First, [Census Genes Mutations](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v100/mutantcensus) for genome GRCh38 was downloaded
2. Then, the data set was reduced to the columns: ```GENE_SYMBOL```, ```MUTATION_CDS```, ```MUTATION_AA```, and ```HGVSG```
3. The reduced data set was converted into a dictionary with the python script ``` /src/prepare_cosmic_resources.py ```
4. The resulted dictionary was saved under the name ```cosmic_hgvsg_dictionary.txt```
5. The path to the dictionary must be provided to the python script ```03_OncoVI_SOP.py```

### ClinVar resources
The download and preparation of [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) resources utilised by the functional annotation STEP is handled by the script ```01_clinvar_resource_manager.sh```.

## To download and prepare the ClinVar resource
```rb
# Run the bash script
bash 01_clinvar_resource_manager.sh
```

## Get started
Clone the GitHub repository:
```rb
git clone https://github.com/MGCarta/oncovi.git
```
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
