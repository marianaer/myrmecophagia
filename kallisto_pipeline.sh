#!/bin/bash
Trinity_fasta=$1 # fasta file with transcripts
fastq1=$2 # Right pair of paired end reads
fastq2=$3 #L eft pair of paired end reads
output_dir=$4 # Path for outputs
TRINITY_PATH=/home/escobar/bin/trinityrnaseq-v2.14.0



# Build kallisto index files from Trinity contigs
kallisto index --make-unique -i ${Trinity_fasta}_kallisto.index ${Trinity_fasta}; done

# Quantifying expression using paired reads on the kallisto index files
kallisto quant -i ${Trinity_fasta}_kallisto.index -o ${Trinity_fasta}_Trinity_kallisto -b 100 ${fastq1} ${fastq2} -t 40


####### Trinity Transcript Quantification #######
# Trinity script for abundance Quantification

## Just prepare the reference for alignment and abundance estimation
${TRINITY_PATH}/util/align_and_estimate_abundance.pl --transcripts ${Trinity_fasta} --est_method kallisto --trinity_mode --prep_reference --output_dir ${output_dir}

## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
${TRINITY_PATH}/util/align_and_estimate_abundance.pl --transcripts ${Trinity_fasta} --seqType fq --left  ${fastq1} --right ${fastq2} --est_method kallisto --trinity_mode --output_dir ${output_dir}

## Use the --gene_trans_map or --trinity_mode parameters in order to get a gene counts matrix in addition to the isoform counts matrix
/media/bigvol/fdelsuc/bin/trinityrnaseq-v2.9.0/util/abundance_estimates_to_matrix.pl --est_method kallisto --gene_trans_map Tamandua_tetradactyla_G2525_Spleen_TAMtetFC04_Trinity.fasta.gene_trans_map --out_prefix Tamandua_tetradactyla_G2525_Spleen_TAMtetFC04 --name_sample_by_basedir Tamandua_tetradactyla_G2525_Spleen_TAMtetFC04/abundance.tsv


${TRINITY_PATH}/util/align_and_estimate_abundance.pl --est_method kallisto --gene_trans_map ${Trinity_fasta}.gene.trans.map --out_prefix ${Trinity_fasta}_kallisto --name_sample_by_base_dir  ${Trinity_fasta}_kallisto/abundance.tsv
