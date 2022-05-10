#!/bin/bash
# Use this script in a loop. Works for individual files

# eg. loop for i in *1_fastp.fq.gz; do sh script.sh "$(echo $i|sed 's/_[0-9]_fastp.fq.gz//g')"_1_fastp.fq.gz  "$(echo $i|sed 's/_[0-9]_fastp.fq.gz//g')"_2_fastp.fq.gz; done

#$1 = xxx_1_fastq.gz file n1
#$2 = xxx_2_fastq.gz file n2

CORENAME=echo $1| sed 's/_[0-9]_fastp.fq.gz//g'

# eg. sh thranscriptome_Assembly. sh xxxx.fastq.gz

# fastp
#fastp -i Carollia_sowelli_SalivaryGland_SRR7064957_1.fastq.gz -o Carollia_sowelli_SalivaryGland_SRR7064957_1_fastp.fq.gz -I Carollia_sowelli_SalivaryGland_SRR7064957_2.fastq.gz -O Carollia_sowelli_SalivaryGland_SRR7064957_2_fastp.fq.gz -w 40 -h Carollia_sowelli_SalivaryGland_SRR7064957_fastp.html -V

#Quality control. FastQC


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
singularity pull busco_5.1.3.sif docker://quay.io/biocontainers/busco:5.1.3--pyhdfd78af_0

# use link to where we downloaded busco originally (previous command)
./busco_5.1.3.sif busco -i "$dir"_trinity_output.Trinity.fasta -l mammalia_odb10 -o "$dir" -m tran

  # Generate Summary plot. IMPORTANT: all summary for all transcriptomes whould be in the same folder. This should be run there

./busco_5.1.3.sif generate_plot.py -wd ${PATH_TO_SUMMARIES}



  #### Get ENSEMBL metazoa peptide database for transcriptome annotation with assembly2ORF
 # We only need the pep  (peptide) files. We cant do ftp with regex insid so we need to do it like this
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-53/fasta/
grep -o '".*"' index.html > metazoa.list
sed 's/"//g' metazoa.list
sed 's/"//g' metazoa.list | sed 's/\///g' > metazoa_ensembl.list

for i in $(cat metazoa_ensembl.list); do wget -r "ftp://ftp.ensemblgenomes.org/pub/metazoa/release-53/fasta/"$i"/pep"


#cat them all together
zcat *pep.all.fa.gz > all_metazoa.pep.all.fa

# Some issues with ascii characters, must be removed
cat all_metazoa.pep.all.fa | perl -ne 's/[^\x00-\x7F]+/ /g; print;' > all_metazoa_db.pep.all.fa


# Create blastdb to blast all vs all_metazoa
/home/escobar/bin/ncbi-blast-2.13.0+/bin/makeblastdb  -in all_metazoa_db.pep.all.fa  -input_type fasta  -dbtype prot

# All vs. all blast
/home/escobar/bin/ncbi-blast-2.13.0+/bin/blastp -db all_metazoa_db.pep.all.fa -query all_metazoa.pep.all.fa -out all_sv_all.tsv -outfmt 6  -num_threads 20





# O tal vez si se podia solo habia que usar comillas dobles en la direccion ftp :(
