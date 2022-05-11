#!/bin/bash

# Transcriptome annotation with Trinnotate

# Follow installation and database acquisition instructions from  https://github.com/Trinotate/Trinotate.github.io/blob/master/index.asciidoc

### Identify candidate coding regions from Trinity outputs with TransDecoder (ie. looking for ORFs) ###

PATH_TO_DB=/home/escobar/Desktop/databases/trinotate_db
PATH_TO_TRANSDECODER=/home/escobar/bin/TransDecoder
Trinity_fasta=$1 #fasta file with transcripts

# 1. Extract long ORFs (at least 100 amino acids long)
#${PATH_TO_TRANSDECODER}/TransDecoder.LongOrfs -t $Trinity_fasta


# 2. Identidy ORFs with homology to known proteins via blast or pfam searches
# Searches candidate peptides from transDecoder for homology using blastp anf pfam

#  BlastP search

#blastp -query ${Trinity_fasta}.transdecoder_dir/longest_orfs.pep -db ${PATH_TO_DB}/uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 24 > blastp.outfmt6

# Pfam search for protein domains
#hmmsearch --cpu 8 --domtblout pfam.domtblout ${PATH_TO_DB}/Pfam-A.hmm transdecoder_dir/longest_orfs.pep


# Integrating Blast and Pfam searches into coding region selection
# Make use of blastp and hmmer searches  and report them as likely coding regions.
#TransDecoder.Predict -t $Trinity_fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

#TransDecoder.Predict -t ${Trinity_fasta}

### Capturing Blast homologies with Trinotate ###

# search Trinity transcripts
blastx -query ${Trinity_fasta} -db ${PATH_TO_DB}/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6

# search Transdecoder-predicted proteins
blastp -query ${Trinity_fasta}.transdecoder_dir/longest_orfs.pep  -db ${PATH_TO_DB}/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6

# Running HMMER to identify protein domains
hmmsearch --cpu 40 --domtblout TrinotatePFAM.out  ${PATH_TO_DB}/Pfam-A.hmm ${Trinity_fasta}.transdecoder_dir/longest_orfs.pep > pfam.log



#Every evening from 9:30 pm on, the  Cinema de la Plage provides open-air screenings for all to enjoy
