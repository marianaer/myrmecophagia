#!/bin/bash

# This script looks for the complete and/or longest isoform based on ORF finding from BORF, and keeps the aminoacid sequence in a file.

#b=$(echo $(basename "$PWD")| sed 's/_kallisto//g')
b=$(echo $(basename "$PWD"))
echo ${b}
TRINITY_PATH=/home/escobar/bin/trinityrnaseq-v2.14.0
blast_results=${b}_muc7_tblastx.outfmt6

cut -f2 ${blast_results} |sort -u > "${b}"_muc7_isoforms_blast.list

# Keep just one isoform (complete and longest, or at least longest or complete)

for i in $(cat *muc7_isoforms_blast.list); do echo $i|sed "s/_i[0-9]//g" ;done



###################33
# look for the genes (only gene id) found with blast and keep all isoform peptides in a document
for i in $(cut -f1-4 ${b}_muc7_isoforms_blast.list -d'_' |sort -u |sed "s/_i[0-9]//g" ); do grep -E -A 1 $i ${b}_trinity_output.Trinity.pep >> ${b}_MUC7.pep;done

# Find if there are incomplete isoforms for such gene identifiers, if so, look for a complete one and delete the incomplete one (then use trinity script to keep the longest )

# List of gene identifiers withat least one incomplete isoform
grep incomplete ${b}_MUC7.pep | cut -f1 -d' '|cut -d'_' -f1-4|sort -u > incomplete_isoforms_muc7.txt

# look if there is also a complete isoform for these genes
# If at least one match containing the word complete was found


# If at least one match containing the word complete was found, remove the incomplete isoforms
for i in $(cat incomplete_isoforms_muc7.txt); do if [ "$(grep -c "$i.* complete" ${b}_trinity_output.Trinity.pep)" -ge 1 ]; then echo $i >> keep_complete.tmp ;fi;done

# remove the incomplete isoforms

for i in $(cat keep_complete.tmp); do grep -A1 "$i.* incomplete" ${b}_MUC7.pep >> remove_seqs.tmp;done

grep -v -f remove_seqs.tmp ${b}_MUC7.pep

rm remove_seqs.tmp; rm keep_complete.tmp

# keep only the longest isoforms
${TRINITY_PATH}/util/misc/get_longest_isoform_seq_per_trinity_gene.pl ${b}_MUC7.pep > ${b}_MUC7_longest_complete.pep
