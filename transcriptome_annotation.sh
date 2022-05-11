# Transcriptome annotation with Trinnotate

# Follow installation and database acquisition instructions from  https://github.com/Trinotate/Trinotate.github.io/blob/master/index.asciidoc

### Identify candidate coding regions from Trinity outputs with TransDecoder ###

PATH_TO_DB=/home/escobar/Desktop/databases/trinotate_db
PATH_TO_TRANSDECODER=/home/escobar/bin/TransDecoder
transcripts_file=$1 #fasta file with transcripts

# 1. Extract long ORFs (at least 100 amino acids long)
${PATH_TO_TRANSDECODER}/TransDecoder.LongOrfs -t $transcripts_file


# 2. Identidy ORFs with homology to known proteins via blast or pfam searches
# Searches candidate peptides from transDecoder for homology using blastp anf pfam

#  BlastP search

blastp -query transdecoder_dir/longest_orfs.pep -db ${PATH_TO_DB}/uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 24 > blastp.outfmt6

# Pfam search for protein domains
hmmscan --cpu 8 --domtblout pfam.domtblout ${PATH_TO_DB}/Pfam-A.hmm transdecoder_dir/longest_orfs.pep
