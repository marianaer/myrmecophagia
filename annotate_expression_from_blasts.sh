#!/bin/bash
#rm kallisto_annotated_muc5b.results
# Extract names of MUC* isoforms found with blast
cut -f1 muc5b_blasts.tsv | sort -u|egrep -o TRINITY_.* > muc5b_isoforms_blast.list

#get kallisto lines for our found mucs
for i in $(cat muc5b_isoforms_blast.list);do grep -E $i"(\^|\s)" kallisto.isoform.TPM.not_cross_norm_annotated >> kallisto_annotated_muc5b.results; done

# remove these lines from kallito file


grep -Fvxf kallisto_annotated_muc5b.results kallisto.isoform.TPM.not_cross_norm_annotated > kallisto.isoform.TPM.not_cross_norm_annotated.tmp

rm kallisto.isoform.TPM.not_cross_norm_annotated; mv kallisto.isoform.TPM.not_cross_norm_annotated.tmp kallisto.isoform.TPM.not_cross_norm_annotated

#ADD MUC7 suffix to transcript names

awk '{ print $1"_MUC5B_BLAST "$2; }' kallisto_annotated_muc5b.results >> kallisto.isoform.TPM.not_cross_norm_annotated
#TRINITY_DN1337_c0_g1_i10


# Extract sequences and store them in a fasta files
grep -w -A 2 -Ff muc5b_isoforms.list $(basename "$PWD")_trinity_output.Trinity.fasta --no-group-separator > $(basename "$PWD")_MUC5B.fasta



cd $i; mv kallisto.isoform.TPM.not_cross_norm_annotated "$i"_kallisto.isoform.TPM.not_cross_norm_annotated


for i in $(ls | grep -v txt|grep -v sh |grep -v txt |grep -v tsv|grep -v fa|grep -v dbs); do cd $i; mv kallisto.isoform.TPM.not_cross_norm_annotated "$i"_kallisto.isoform.TPM.not_cross_norm_annotated; cd ..;done
