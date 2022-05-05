#!/bin/bash
# Use this script in a loop. Works for individual files

# eg. loop for i in *1_fastp.fq.gz; do sh script.sh "$(echo $i|sed 's/_[0-9]_fastp.fq.gz//g')"_1_fastp.fq.gz  "$(echo $i|sed 's/_[0-9]_fastp.fq.gz//g')"_2_fastp.fq.gz; done

#$1 = xxx_1_fastq.gz file n1
#$2 = xxx_2_fastq.gz file n2

CORENAME=echo $1| sed 's/_[0-9]_fastp.fq.gz//g'

# eg. sh thranscriptome_Assembly. sh xxxx.fastq.gz

# Quality control. FastQC

# Create dorectpry for fastqc outputs if it doesnt already exist
mkdir -p fastq_outputs
mkdir -p fastq_outputs/$CORENAME_fastqc

fastqc $1 --outdir fastqc_outputs/$CORENAME_fastqc --threads 5
fastqc $2 --outdir fastqc_outputs/$CORENAME_fastqc --threads 5

# MultiQC of both ends
python3 -m multiqc fastqc_outputs/*/*


# remove erroneous kmers with r Corrector


#perl /home/escobar/bin/rcorrector/run_rcorrector.pl -1 $1 -2 $2 -t 12


# Discard pairs for which one of the reads is deemed unfixable (also strips 'cor' from the header)

#cor1 = $(echo $1 | sed 's/fastp.fq.gz/fastp.cor.fq.gz/g')
#cor2 =$(echo $2 | sed 's/fastp.fq.gz/fastp.cor.fq.gz/g')

#python ~/bin/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 $cor1 -2 $cor2 -s $CORENAME


# Trinity Assembly

dir_list=$(ls|grep -v fast|grep -v .sh)

for dir in $dir_list
  do
    cd $dir
    echo $dir

    gdown $(cat fq1.gid)
    gdown $(cat fq2.gid)

    fq1=$(ls|grep *1_fastp.fq.gz)
    fq2=$(ls|grep *2_fastp.fq.gz)

    echo 'Starting assembly...'
    Trinity --max_memory 30G --seqType fq  --left $fq1 --right $fq2 --CPU 24 --full_cleanup --output "$dir"_trinity_output  > "$dir"_trinity.log

    rm $fq1
    rm $fq2
    rm -rf "$dir"_trinity_output

    # Get trinity assembly metrics

    echo 'Calculating assembly metrics...'

    ~/bin/trinityrnaseq-v2.14.0/util/TrinityStats.pl "$dir"_trinity_output.Trinity.fasta > "$dir"_trinity_assemlby.metrics

    cd ..
  done



  ## Running busco

  ./busco_5.1.3.sif busco
