#!/bin/bash

for i in $(ls|grep -v Carollia|grep -v sh)
  do
    cd $i

    echo $i
    trinity_fa=$(ls|grep Trinity.fasta|grep -v gene)
    trinity_map=$(ls|grep Trinity.fasta.gene)

    #sh /home/escobar/Desktop/RNA_seq_data/trinity_assemblies/test.sh "$trinity_fa" "$trinity_map"
    sh /home/escobar/github/myrmecophagia/trinotate_annotation.sh "$trinity_fa" "$trinity_map"
    cd ..
  done
