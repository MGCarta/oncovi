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
OncoVI was implemented and tested on a dedicated conda enviroment running on a remote server based on Ubuntu's 20.04.4 long-term support (LTS) operating system. To run OncoVI be sure that the following packages are installed:

* python
* numpy
* pandas
* subprocess

### COSMIC resources
Due to size constraints the [COSMIC](https://cancer.sanger.ac.uk/cosmic/download/cosmic) resources utilised by OncoVI cannot be uploaded on the GitHub repo. However, here is a user-friendly guide on how to download and handle cosmic data to make it usable by OncoVI.

#### Cancer Mutation Census 

1. [All Data CMC](https://cancer.sanger.ac.uk/cosmic/download/cancer-mutation-census/v100/alldata-cmc) was first downloaded;
2. The data set was reduced to columns: GENE_NAME, MUTATION CDS, MUTATION AA, COSMIC_SAMPLE_MUTATED, Mutation genome position GRCh38;
3. The reduced version of the data set was read and converted to a dictionary through the python script ```rb # /src/prepare_resources.py ```


## Get started
Clone the GitHub repository:
```rb
git clone XXX

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
python -i /testdata/SOP_table_union.txt -r /resources/VEP
```

#### Run OncoVI

```rb
# Navigate to the directory in which the python script 03_OncoVI_SOP.py is located
cd /path/to/03_OncoVI_SOP.py

# Run OncoVI
python 03_OncoVI_SOP.py -i /testdata/VEP/SOP_table_union_prediction.xlsx -r /resources/OncoVI
```


## License
OncoVI is intended for research purposes only and its use outside of this context is under the responsibility of the user, who should also comply with licences of the resources utilised.

## References
Please cite our preprint [Oncogenicity Variant Interpreter (OncoVI): oncogenicity guidelines implementation to support somatic variants interpretation in precision oncology](https://www.medrxiv.org/content/10.1101/2024.10.10.24315072v1) if you decide to use OncoVI.


# H1
## H2
### H3
#### H4
##### H5
###### H6

Alternatively, for H1 and H2, an underline-ish style:

Alt-H1
======

Alt-H2
------

Emphasis, aka italics, with *asterisks* or _underscores_.

Strong emphasis, aka bold, with **asterisks** or __underscores__.

Combined emphasis with **asterisks and _underscores_**.

Strikethrough uses two tildes. ~~Scratch this.~~

1. First ordered list item
2. Another item
⋅⋅* Unordered sub-list. 
1. Actual numbers don't matter, just that it's a number
⋅⋅1. Ordered sub-list
4. And another item.

⋅⋅⋅You can have properly indented paragraphs within list items. Notice the blank line above, and the leading spaces (at least one, but we'll use three here to also align the raw Markdown).

⋅⋅⋅To have a line break without a paragraph, you will need to use two trailing spaces.⋅⋅
⋅⋅⋅Note that this line is separate, but within the same paragraph.⋅⋅
⋅⋅⋅(This is contrary to the typical GFM line break behaviour, where trailing spaces are not required.)

* Unordered list can use asterisks
- Or minuses
+ Or pluses

[I'm an inline-style link](https://www.google.com)

[I'm an inline-style link with title](https://www.google.com "Google's Homepage")

[I'm a reference-style link][Arbitrary case-insensitive reference text]

[I'm a relative reference to a repository file](../blob/master/LICENSE)

[You can use numbers for reference-style link definitions][1]

Or leave it empty and use the [link text itself].

URLs and URLs in angle brackets will automatically get turned into links. 
http://www.example.com or <http://www.example.com> and sometimes 
example.com (but not on Github, for example).

Some text to show that the reference links can follow later.

Here's our logo (hover to see the title text):

Inline-style: 
![alt text](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Logo Title Text 1")

Reference-style: 
![alt text][logo]

[logo]: https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Logo Title Text 2"

[arbitrary case-insensitive reference text]: https://www.mozilla.org
[1]: http://slashdot.org
[link text itself]: http://www.reddit.com
