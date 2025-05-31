#!/bin/bash

# The script takes in input the txt/VCF file and runs VEP

FULL_PATH_INPUT_DATA=$1
PARENT_DIR_INPUT_DATA=$2
OUTPUT_NAME=$3
VEP_DIR=$4

# Move to the directory in which the shell script is located
cd "$(dirname "$0")"
echo "VEP annotation of file: $FULL_PATH_INPUT_DATA"
echo " "

# Verify that the file exists and has the correct extension

if [ -f "$FULL_PATH_INPUT_DATA" ] && ([[ "$FULL_PATH_INPUT_DATA" == *.txt ]] || [[ "$FULL_PATH_INPUT_DATA" == *.vcf ]]);
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
	--fasta $VEP_DIR/homo_sapiens/114_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
	--input_file $FULL_PATH_INPUT_DATA \
	--plugin dbNSFP,$VEP_DIR/Plugins/dbNSFP/dbNSFP4.5a_grch38.gz,pep_match=0,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate_rankscore \
	--plugin SpliceAI,snv=$VEP_DIR/Plugins/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=$VEP_DIR/Plugins/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
	--output_file $PARENT_DIR_INPUT_DATA/vep_${OUTPUT_NAME}

else
	echo "ERROR: the input file has not either .txt or .vcf extension!"

fi

echo ""
echo "process done!"
