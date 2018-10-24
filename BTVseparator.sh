#!/bin/bash

##run in folder with sequencing files, give root of sequencing files as first arg, all BTV reference segment fasta as second arg

logfile="$1"_assembly.log

exec > $logfile 2>&1

if [ $# -gt 0 ]; then
    echo "Your command line contains $# arguments. Two expected."
else
    echo "Your command line contains no arguments. Two expected."
fi

##make subdirectory for working, copy required files and convert from bz2 to gz
mkdir "$1"_sep

bunzip2 -ck "$1"_R1_001.fastq.bz2 > "$1"_1.fq
bunzip2 -ck "$1"_R2_001.fastq.bz2 > "$1"_2.fq

##remove adaptors, short and low quality sequence, then run fastqc
trim_galore -o "$1"_sep -q 30 --length 50 -r1 100 -r2 100 --fastqc --illumina --paired --retain_unpaired --phred33 "$1"_1.fq "$1"_2.fq
rm -f "$1"_1.fq "$1"_2.fq

##use FLASH to join overlapping paired end reads
cd "$1"_sep
mv "$1"_1_val_1.fq "$1"_1.fq
mv "$1"_2_val_2.fq "$1"_2.fq

##mapping cleaned reads to all BTV subtypes using bowtie2 - eXpress requires sort by name
/software/bowtie2-2.3.4.2/bowtie2-build -f ../"$2" "$2".ref
/software/bowtie2-2.3.4.2/bowtie2 -a --phred33 --no-unal --very-sensitive --rg-id "$1" -x "$2".ref -1 "$1"_1.fq -2 "$1"_2.fq -S "$1"_allBTV_bowtie_paired.sam --al-conc "$1"_BTV 
rm -f *fastqc* *bt2 "$1"_1.fq_trimming_report.txt "$1"_2.fq_trimming_report.txt "$1"_1_unpaired_1.fq  "$1"_2_unpaired_2.fq "$1"_allBTV_bowtie_paired.sam "$1"_1.fq "$1"_2.fq
mv "$1"_BTV.1 "$1"_BTV.1.fq
mv "$1"_BTV.2 "$1"_BTV.2.fq
mkdir Assembly
/software/SPAdes-3.9.0-Linux/bin/metaspades.py -1 "$1"_BTV.1.fq -2 "$1"_BTV.2.fq -o ./Assembly/
rm -f "$1"_BTV.1.fq "$1"_BTV.2.fq
