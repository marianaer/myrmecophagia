#!/bin/bash
TRINITY_PATH=/home/escobar/bin/trinityrnaseq-v2.14.0

#b=$(echo $(basename "$PWD")| sed 's/_kallisto//g')
b=$(echo $(basename "$PWD"))
echo ${b}
blast_results=${b}_MUC5B_tblastx.outfmt6
echo ${blast_results}
#sed -i "s/>/>${b}_/g" ${b}.fasta
# remove line breaks from Fasta transcriptome files ins equences
#awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${b}*.Trinity.fasta > ${b}*.Trinity.fasta.tmp

#mv ${b}*.Trinity.fasta.tmp ${b}*.Trinity.fasta
# extract Annotated mucin sequences
# Get isoforms found via blast
cut -f2 ${blast_results} |sort -u > "${b}"_MUC5B_isoforms_blast.list

echo ${b}

rm ${b}_MUC5B.fasta
#for i in $(cat *MUC5B_isoforms*list); do grep -E -A 1 "$i len" ${b}_trinity_output.Trinity.fasta >> ${b}_MUC5B.fasta;done

for i in $(cat *MUC5B_isoforms*list); do grep -E -A 1 $i. ${b}_trinity_output.Trinity.pep >> ${b}_MUC5B.pep;done

#add name of organism to sequence name
#sed -i "s/TRINITY/${b}_TRINITY/g" ${b}_MUC5B.fasta

# Get longest isoforms per gene
~/bin/trinityrnaseq-v2.14.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl ${b}_MUC5B.fasta  > ${b}_MUC5B_longest.fasta


# Keep just one isoform (complete and longest, or at least longest or complete)

for i in $(cat *muc7_isoforms*list); do echo $i|sed "s/_i[0-9]//g" ;done

# Finding the exact gene identifier (no isoform)
cut -f1-4 Carollia_sowelli_muc7_isoforms_blast.list -d'_' |sort -u



###################33
# look for the genes found with blast and keep all isoform peptides in a document
for i in $(cut -f1-4 Carollia_sowelli_muc7_isoforms_blast.list -d'_' |sort -u |sed "s/_i[0-9]//g" ); do grep -E -A 1 $i ${b}_trinity_output.Trinity.pep >> ${b}_MUC7.pep;done

# Find if there are incomplete isoforms for such gene identifiers, if so, look for a complete one and delete the incomplete one (then use trinity script to keep the longest )

# List of gene identifiers withat least one incomplete isoform
grep incomplete Carollia_sowelli_MUC7.pep | cut -f1 -d' '|cut -d'_' -f1-4|sort -u >incomplete_isoforms_muc7.txt

# look if there is also a complete isoform for these genes
# If at least one match containing the word complete was found
for i in $(cat incomplete_isoforms_muc7.txt); do if [ "$(grep -c "$i.* complete" Carollia_sowelli_trinity_output.Trinity.pep)" -ge 1 ]; then echo $i;fi;done

# If at least one match containing the word complete was found, remove the incomplete isoforms
for i in $(cat incomplete_isoforms_muc7.txt); do if [ "$(grep -c "$i.* complete" Carollia_sowelli_trinity_output.Trinity.pep)" -ge 1 ]; then echo $i >>keep_complete.tmp ;fi;done

# remove the incomplete isoforms

for i in $(cat keep_complete.tmp); do grep -A1 "$i.* incomplete" Carollia_sowelli_MUC7.pep >> remove_seqs.tmp;done

grep -v -f remove_seqs.tmp Carollia_sowelli_MUC7.pep

rm remove_seqs.tmp; rm keep_complete.tmp

# keep only the longest isoforms
${TRINITY_PATH}/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Carollia_sowelli_MUC7.pep > Carollia_sowelli_MUC7_longest_complete.pep
