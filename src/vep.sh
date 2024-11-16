#!/bin/bash


# shell script which takes as input the txt/VCF file provided by the 
# 02_VEP_based_pipeline.py script and runs VEP

# Move to the directory in which the shell script is located 

cd "$(dirname "$0")"
echo "VEP annotation of file $1"
echo " "

# Verify that the file exists and has the correct extension
# Modify here ".vcf" according to file extension

if [ -f $1 ] && [ ${1: -4} == ".txt" ]

then

	vep --offline \
	--cache \
	--fork 4 \
	--assembly GRCh38 \
	--refseq \
	--force_overwrite \
	--symbol \
	--canonical \
	--mane \
	--hgvs \
	--hgvsg \
	--hgvsg_use_accession \
	--numbers \
	--no_intergenic \
	--clin_sig_allele 1 \
	--af \
	--af_1kg \
	--af_gnomade \
	--af_gnomadg \
	--pubmed \
	--domains \
	--vcf \
	--dir $4 \
	--fasta /path/to/.vep/homo_sapiens/111_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
	--input_file $1 \
	--plugin dbNSFP,/path/to/.vep/Plugins/dbNSFP/dbNSFP4.5a_grch38.gz,pep_match=0,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate_rankscore \
	--plugin SpliceAI,snv=/path/to/.vep/Plugins/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/path/to/.vep/Plugins/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
	--output_file $2/vep_${3}

else 

	echo "THE FILE HAS NOT THE EXPECTED EXTENSION!"
	
fi 

echo " "
echo "process done!"