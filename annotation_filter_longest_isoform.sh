#!/bin/bash

Trinity_fasta=$1 #fasta file with transcripts
#Trinity_gene_map=$2 #Gene_trans_map output by trinity
PATH_TO_DB=/home/escobar/Desktop/databases/trinotate_db
TRINOTATE_PATH=/home/escobar/bin/Trinotate
RNAMMER_PATH=/home/escobar/bin/rnammer-1.2
TRINITY_PATH=/home/escobar/bin/trinityrnaseq-v2.14.0
PATH_TO_TRANSDECODER=/home/escobar/bin/TransDecoder
PATH_TO_TMHMM=/home/escobar/bin/tmhmm-2.0c/bin

# Filter trinity assembly by the longest idoform for each gene (has the most exons)
${TRINITY_PATH}/util/misc/get_longest_isoform_seq_per_trinity_gene.pl ${Trinity_fasta} > ${Trinity_fasta}.longest.fasta

# Filter the gene_trans map as well (2nd column)


${TRINITY_HOME}/util/support_scripts/get_Trinity_gene_to_trans_map.pl ${Trinity_fasta}.longest.fasta >  ${Trinity_fasta}.longest.fasta.gene_trans_map


# 1. Extract long ORFs (at least 100 amino acids long)
echo 'TransDecoder LongORFs...'
echo 'command:' "${PATH_TO_TRANSDECODER}"/TransDecoder.LongOrfs -t ${Trinity_fasta}.longest.fasta
${PATH_TO_TRANSDECODER}/TransDecoder.LongOrfs -t ${Trinity_fasta}.longest.fasta

echo 'TransDecoder Predict...'
${PATH_TO_TRANSDECODER}/TransDecoder.Predict -t ${Trinity_fasta}.longest.fasta

# 2.
### Capturing Blast homologies with Trinotate ###

# search Trinity transcripts
blastx -query ${Trinity_fasta}.longest.fasta -db ${PATH_TO_DB}/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

# search Transdecoder-predicted proteins
blastp -query ${Trinity_fasta}.longest.fasta.transdecoder_dir/longest_orfs.pep  -db ${PATH_TO_DB}/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

# Running HMMER to identify protein domains
hmmsearch --cpu 40 --domtblout TrinotatePFAM.out  ${PATH_TO_DB}/Pfam-A.hmm ${Trinity_fasta}.longest.fasta.transdecoder_dir/longest_orfs.pep > pfam.log

## SignalP to predict signal peptides
signalp6 --output_dir $(pwd) --mode fast --format txt --organism eukarya --fastafile ${Trinity_fasta}.longest.fasta.transdecoder.pep

#DeepTMHMM (tmHMM) to predict transmembrane regions
# conda activate tmhmm
#biolib run DTU/DeepTMHMM --fasta ${Trinity_fasta}.longest.fasta.transdecoder.pep
${PATHO_TO_TMHMM}/tmhmm --short < ${Trinity_fasta}.longest.fasta.transdecoder.pep > tmhmm.out
#conda deactivate

#RNAMMER to identify rRNA transcripts

# First, concatenate all transcripts together into a super-scaffold. Then run rnammer to identify rRNA homologies and transform rRNA feature cordinates in the big scaffold back to the transcriptome reference coordinates.
${TRINOTATE_PATH}/util/rnammer_support/RnammerTranscriptome.pl --transcriptome ${Trinity_fasta}.longest.fasta --path_to_rnammer ${RNAMMER_PATH}/rnammer


# Load resulting gff file into trinotate in SQLite databases

# Getting boilerplate Trinotate sqlite db and populating it with our data
#d 1. Load transcripts and coding regions into sqlite db
conda activate trinotate
${TRINOTATE_PATH}/Trinotate ${PATH_TO_DB}/Trinotate.sqlite init --gene_trans_map ${Trinity_fasta}.longest.fasta.gene_trans_map --transcript_fasta ${Trinity_fasta}.longest.fasta --transdecoder_pep ${Trinity_fasta}.longest.fasta.transdecoder.pep

# 2. Load Blast homologies
# load protein hits
${TRINOTATE_PATH}/Trinotate ${PATH_TO_DB}/Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6

# load transcript hits
${TRINOTATE_PATH}/Trinotate ${PATH_TO_DB}/Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6

# 3. Load Pfam domain entries
${TRINOTATE_PATH}/Trinotate ${PATH_TO_DB}/Trinotate.sqlite LOAD_pfam TrinotatePFAM.out

# 4. Load transmembrane domains
${TRINOTATE_PATH}/Trinotate ${PATH_TO_DB}/Trinotate.sqlite  LOAD_tmhmm tmhmm.out

# 5. Load signal peptide predictions
${TRINOTATE_PATH}/Trinotate ${PATH_TO_DB}/Trinotate.sqlite LOAD_signalp prediction_results.txt

# Trinotate output and annotation Report
${TRINOTATE_PATH}/Trinotate ${PATH_TO_DB}/Trinotate.sqlite report --incl_pep --incl_trans > trinotate_annotation_report.xls

conda deactivate
