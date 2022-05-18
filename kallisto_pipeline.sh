#!/bin/bash
Trinity_fasta=$1 # fasta file with transcripts
fastq1=$2 # Right pair of paired end reads
fastq2=$3 #L eft pair of paired end reads

TRINITY_PATH=/home/escobar/bin/trinityrnaseq-v2.14.0



# Build kallisto index files from Trinity contigs
echo 'Building kallisto index files from Trinity contigs...'
kallisto index --make-unique -i ${Trinity_fasta}_kallisto.index ${Trinity_fasta}

# Quantifying expression using paired reads on the kallisto index files
echo 'Quantifying expression using paired reads on the kallisto index files...'
kallisto quant -i ${Trinity_fasta}_kallisto.index -o ${Trinity_fasta}_Trinity_kallisto -b 100 ${fastq1} ${fastq2} -t 40


####### Trinity Transcript Quantification #######
# Trinity script for abundance Quantification

## Just prepare the reference for alignment and abundance estimation
echo 'Preparing reference for alignment...'

${TRINITY_PATH}/util/align_and_estimate_abundance.pl --transcripts ${Trinity_fasta} --est_method kallisto --trinity_mode --prep_reference --output_dir ${Trinity_fasta}_Trinity_kallisto

## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
echo 'Aligning and estimating abundance...'
${TRINITY_PATH}/util/align_and_estimate_abundance.pl --transcripts ${Trinity_fasta} --seqType fq --left  ${fastq1} --right ${fastq2} --est_method kallisto --trinity_mode --output_dir ${Trinity_fasta}_Trinity_kallisto

## Use the --gene_trans_map or --trinity_mode parameters in order to get a gene counts matrix in addition to the isoform counts matrix
echo 'Getting counts matrix...'
echo ${Trinity_fasta}_Trinity_kallisto/abundance.tsv
${TRINITY_PATH}/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir ${Trinity_fasta}_Trinity_kallisto/abundance.tsv --gene_trans_map ${Trinity_fasta}.gene_trans_map


#/home/escobar/bin/trinityrnaseq-v2.14.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir kallisto_output/abundance.tsv  --gene_trans_map Choloepus_didactylus_trinity_output.Trinity.fasta.longest.fasta.gene_trans_map
#http://sepsis-omics.github.io/tutorials/modules/kallisto/
