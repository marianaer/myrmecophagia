#!/bin/bash

# Transcriptome annotation with Trinnotate

# Follow installation and database acquisition instructions from  https://github.com/Trinotate/Trinotate.github.io/blob/master/index.asciidoc

### Identify candidate coding regions from Trinity outputs with TransDecoder (ie. looking for ORFs) ###

Trinity_fasta=$1 #fasta file with transcripts
Trinity_gene_map=$2 #Gene_trans_map output by trinity
PATH_TO_DB=/home/escobar/Desktop/databases/trinotate_db
TRINOTATE_PATH=/home/escobar/bin/Trinotate
RNAMMER_PATH=/home/escobar/bin/rnammer-1.2

# 1. Extract long ORFs (at least 100 amino acids long)
#${PATH_TO_TRANSDECODER}/TransDecoder.LongOrfs -t $Trinity_fasta

#TransDecoder.Predict -t ${Trinity_fasta}

# 2.
### Capturing Blast homologies with Trinotate ###

# search Trinity transcripts
blastx -query ${Trinity_fasta} -db ${PATH_TO_DB}/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

# search Transdecoder-predicted proteins
blastp -query ${Trinity_fasta}.transdecoder_dir/longest_orfs.pep  -db ${PATH_TO_DB}/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

# Running HMMER to identify protein domains
hmmsearch --cpu 40 --domtblout TrinotatePFAM.out  ${PATH_TO_DB}/Pfam-A.hmm ${Trinity_fasta}.transdecoder_dir/longest_orfs.pep > pfam.log

## SignalP to predict signal peptides
signalp6 -f short -n signalp.out ${Trinity_fasta}.transdecoder.pep

#DeepTMHMM (tmHMM) to predict transmembrane regions
# conda activate tmhmm
biolib run DTU/DeepTMHMM --fasta ${Trinity_fasta}.transdecoder.pep
#tmhmm --short < transdecoder.pep > tmhmm.out
#conda deactivate

#RNAMMER to identify rRNA transcripts

# First, concatenate all transcripts together into a super-scaffold. Then run rnammer to identify rRNA homologies and transform rRNA feature cordinates in the big scaffold back to the transcriptome reference coordinates.
${TRINOTATE_PATH}/util/rnammer_support/RnammerTranscriptome.pl --transcriptome ${Trinity_fasta} --path_to_rnammer ${RNAMMER_PATH}/rnammer


# Load resulting gff file into trinotate in SQLite databases

# Import transcriptome and protein data
init --gene_trans_map ${Trinity_gene_map} --transcript_fasta ${Trinity_fasta} --transdecoder_pep ${Trinity_fasta}.transdecoder.pep

# TransDecoder protein search results
#LOAD_swissprot_blastp blastp.outfmt6
#LOAD_pfam TrinotatePFAM.out
#LOAD_tmhmm <file>
#LOAD_signalp <file>

#${TRINOTATE_PATH}/Trinotate
