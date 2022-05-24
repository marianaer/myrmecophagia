#!/bin/bash

b=$(basename "$PWD")
# remove line breaks from Fasta transcriptome files ins equences
#awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${b}*.Trinity.fasta > ${b}*.Trinity.fasta.tmp

#mv ${b}*.Trinity.fasta.tmp ${b}*.Trinity.fasta
# extract Annotated mucin sequences

rm ${b}_MUC5B.fasta
for i in $(cat *muc5b_isoforms*list); do grep -E -A 1 "$i len" ${b}_trinity_output.Trinity.fasta >> ${b}_MUC5B.fasta;done

#add name of organism to sequence name
sed -i "s/TRINITY/${b}_TRINITY/g" ${b}_MUC5B.fasta
