# Oncogenicity Variant Interpreter (OncoVI)
OncoVI is a fully-automated Python implementation of the [oncogenicity guidelines](https://pubmed.ncbi.nlm.nih.gov/35101336/) published in 2022 by Horak et al. 

Starting from the genomic location of the variants OncoVI:
1. Performs functional annotation based on the [Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) from Ensembl 
2. Collects biological evidences from the implemented publicly available resources
3. Provides a classification of oncogenicity based on the oncogenicity guidelines.

Workflow of OncoVI: 
![alt text][logo]

[logo]: https://github.com/MGCarta/oncovi/blob/main/figures/Figure1.PNG "Logo Title Text 2"

###### The figure shows the implemented criteria in OncoVI (11 and five criteria for evidence of oncogenic and benign effect respectively), the public resources utilised to assess each criterion, the points associated with each criterion and the classification of oncogenicity into one of five classes on the basis of the variant-specific score, obtained as the sum of the points associated to the criteria triggered by OncoVI for the variant: score≥10:Oncogenic (O), 6≤score≤9:Likely Oncogenic (LO), 0≤score≤5:Variant of uncertain significance (VUS), -6≤score≤-1:Likely Benign (LB), score≤-7:Benign (B). Blue: resources suggested by the Standard Operating Procedure, black: resources identified by the authors of this study.

## Software requirements
OncoVI was implemented and tested on a dedicated conda enviroment running on a remote server based on Ubuntu's 20.04.4 long-term support (LTS) operating system. To run OncoVI the following packages are required:

* python
* numpy
* pandas
* subprocess

### COSMIC resources
Due to size constraints the [COSMIC](https://cancer.sanger.ac.uk/cosmic/download/cosmic) resources utilised by OncoVI cannot be uploaded on the GitHub repo. Here, a user-friendly guide on how to download and handle cosmic data to make it usable by OncoVI.

#### Cancer Mutation Census 

1. First, [All Data CMC](https://cancer.sanger.ac.uk/cosmic/download/cancer-mutation-census/v100/alldata-cmc) for genome GRCh38 was downloaded;
2. Then, the data set was reduced to columns: ```GENE_NAME```, ```MUTATION CDS```, ```MUTATION AA```, ```COSMIC_SAMPLE_MUTATED```, ```Mutation genome position GRCh38```;
3. The reduced data set was converted into a dictionary with the python script ``` /src/prepare_resources.py ```;
4. The resulted dictionary is saved under the name ```cosmic_hgvsg_dictionary.txt``` and must be provided to OncoVI for the evaluation of the oncogenicity guidelines.

#### Census Genes Mutations 

1. First, [Census Genes Mutations](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v100/mutantcensus) for genome GRCh38 was downloaded;
2. Then, the data set was reduced to columns: ```GENE_SYMBOL```, ```MUTATION_CDS```, ```MUTATION_AA```, and ```HGVSG```;
3. The reduced data set was converted into a dictionary with the python script ``` /src/prepare_resources.py ```;
4. The resulted dictionary is saved under the name ```cosmic_all_dictionary.txt``` and must be provided to OncoVI for the evaluation of the oncogenicity guidelimes. 

### ClinVar resources
Due to size constraints the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) resources utilised by the functional annotation STEP cannot be uploaded on the GitHub repo. Here, a user-friendly guide on how to download and handle ClinVar data to make it usable by the functional annotation STEP.

1. First, ```variant_summary.txt.gz``` was downloaded from the [ftp site](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/);
2. Then, the data set was reduced to columns: ```GeneSymbol```, ```ClinicalSignificance```, ```Chromosome```, ```Start```, ```VariationID```, ```ReferenceAlleleVCF```, ```AlternateAlleleVCF```, ```ReviewStatus```, ```NumberSubmitters```;
3. The reduced data set was converted into a dictionary with the python script ``` /src/create_clinvar_dict.py ```;
4. The resulted dictionary is saved under the name ```clinvar_all_dictionary.txt``` and must be provided to the functional annotation STEP. 

## Get started
Clone the GitHub repository:
```rb
git clone https://github.com/MGCarta/oncovi.git

# Create the conda environment oncovi
conda env create -f /path/to/OncoVIenvFile.yml

# Activate the conda environment
conda activate oncovi
```

## Prepare your variants
Variants in text and variant call format (VCF) are both accepted by [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html). Please refere to VEP´s [official documentation](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#input) for a detailed description of input formats.
A test data is available under:

```rb
# /testdata/SOP_table_union.txt
```

#### Perform functional annotation with VEP

```rb
# Navigate to the directory in which the python script 02_VEP_based_pipeline.py is located
cd /path/to/02_VEP_based_pipeline.py

# Run the functional annotation
python -i /testdata/SOP_table_union.txt
```

#### Run OncoVI

```rb
# Navigate to the directory in which the python script 03_OncoVI_SOP.py is located
cd /path/to/03_OncoVI_SOP.py

# Run OncoVI
python 03_OncoVI_SOP.py -i /testdata/VEP/SOP_table_union_prediction.xlsx -r /resources
```


## License
OncoVI is intended for research purposes only and its use outside of this context is under the responsibility of the user, who should also comply with licences of the resources utilised.

## References
Please cite our preprint [Oncogenicity Variant Interpreter (OncoVI): oncogenicity guidelines implementation to support somatic variants interpretation in precision oncology](https://www.medrxiv.org/content/10.1101/2024.10.10.24315072v1) if you decide to use OncoVI.
