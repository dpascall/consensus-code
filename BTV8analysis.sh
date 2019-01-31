#!/bin/bash

##run in folder with sequencing files, give root of sequencing files as first arg, all BTV reference segment fasta as second arg, BTV-8 reference segment fasta as third arg

#########
###NOTE: assumes no mixed infection and all classical BTV-8 segments - check all BTV mapping results carefully to confirm this before using output
#########

logfile="$1".log

exec > $logfile 2>&1

if [ $# -gt 0 ]; then
    echo "Your command line contains $# arguments. Three expected."
else
    echo "Your command line contains no arguments. Three expected."
fi

##make subdirectory for working, copy required files and convert from bz2 to gz
mkdir "$1"_working
cp $3 "$1"_working/"$3"_working.fasta
cp $3 "$1"_working/"$3"

bunzip2 -ck "$1"_R1_001.fastq.bz2 > "$1"_1.fq
bunzip2 -ck "$1"_R2_001.fastq.bz2 > "$1"_2.fq

##remove adaptors, short and low quality sequence, then run fastqc
trim_galore -o "$1"_working -q 30 --length 50 -r1 100 -r2 100 --fastqc --illumina --paired --retain_unpaired --phred33 "$1"_1.fq "$1"_2.fq
rm -f "$1"_1.fq "$1"_2.fq

##use FLASH to join overlapping paired end reads
cd "$1"_working
mv "$1"_1_val_1.fq "$1"_1.fq
mv "$1"_2_val_2.fq "$1"_2.fq
/software/FLASH-1.2.11/flash -M 130 "$1"_1.fq "$1"_2.fq

##merge singletons
cat out.extendedFrags.fastq "$1"_1_unpaired_1.fq "$1"_2_unpaired_2.fq > "$1".mergedsingletons.fq
cat "$1"_1_unpaired_1.fq "$1"_2_unpaired_2.fq > "$1".mergedoriginalsingletons.fq
rm -f out.extendedFrags.fastq "$1"_1_unpaired_1.fq "$1"_2_unpaired_2.fq

##mapping cleaned reads to all BTV subtypes using bowtie2 - eXpress requires sort by name
/software/bowtie2-2.3.4.2/bowtie2-build -f ../"$2" "$2".ref
/software/bowtie2-2.3.4.2/bowtie2 -a --phred33 --no-unal --very-sensitive --rg-id "$1" -x "$2".ref -1 out.notCombined_1.fastq -2 out.notCombined_2.fastq -S "$1"_allBTV_bowtie_paired.sam
/software/bowtie2-2.3.4.2/bowtie2 -a --phred33 --no-unal --very-sensitive --rg-id "$1" -x "$2".ref -U "$1".mergedsingletons.fq -S "$1"_allBTV_bowtie_unpaired.sam
samtools view -bS "$1"_allBTV_bowtie_paired.sam > "$1"_allBTV_bowtie_paired.bam
samtools view -bS "$1"_allBTV_bowtie_unpaired.sam > "$1"_allBTV_bowtie_unpaired.bam
rm -f "$1"_allBTV_bowtie_paired.sam "$1"_allBTV_bowtie_unpaired.sam
samtools sort "$1"_allBTV_bowtie_unpaired.bam -o "$1"_allBTV_bowtie_unpaired.bam
samtools index "$1"_allBTV_bowtie_unpaired.bam
samtools sort "$1"_allBTV_bowtie_paired.bam -o "$1"_allBTV_bowtie_paired.bam
samtools index "$1"_allBTV_bowtie_paired.bam
samtools merge "$1"_allBTV_bowtie.bam "$1"_allBTV_bowtie_paired.bam "$1"_allBTV_bowtie_unpaired.bam
samtools sort -n "$1"_allBTV_bowtie.bam -o "$1"_allBTV_bowtie.bam
samtools index "$1"_allBTV_bowtie.bam
rm -f "$1"_allBTV_bowtie_paired.bam "$1"_allBTV_bowtie_unpaired.bam "$1"_allBTV_bowtie_paired.bam.bai "$1"_allBTV_bowtie_unpaired.bam.bai

##statistics and coverage plot of mapped reads all BTV + per segment mapping stats for original and de-duped
mkdir allBTV

##requires eXpress - calculate per segment mapping information
express --no-update-check -L 500 -o "./allBTV/" ../$2 "$1"_allBTV_bowtie.bam

##re-sort by location
samtools sort "$1"_allBTV_bowtie.bam -o "$1"_allBTV_bowtie.bam
samtools index "$1"_allBTV_bowtie.bam

##generate weeSAM stats
cd allBTV
Rscript ../../figure.R
weeSAM --bam ../"$1"_allBTV_bowtie.bam --html ./allBTV/"$1"_allBTV.html --out output.txt
cd ..

##remove all BTV mapping files
rm -f "$1"_allBTV_bowtie.bam

for i in {1..5}
do
	if [ "$i" -eq 1 ]
	then
		echo "1"
		##mapping cleaned reads to BTV-8 reference genome using tanoti
		tanoti -u 0 -r $3 -i out.notCombined_1.fastq out.notCombined_2.fastq -o "$1"_BTV8_tanoti_paired.sam -p 1
		tanoti -u 0 -r $3 -i "$1".mergedsingletons.fq -o "$1"_BTV8_tanoti_unpaired.sam
		samtools view -bS "$1"_BTV8_tanoti_paired.sam > "$1"_BTV8_tanoti_paired.bam
		samtools view -bS "$1"_BTV8_tanoti_unpaired.sam > "$1"_BTV8_tanoti_unpaired.bam
		rm -f "$1"_BTV8_tanoti_paired.sam "$1"_BTV8_tanoti_unpaired.sam
		samtools sort "$1"_BTV8_tanoti_paired.bam -o "$1"_BTV8_tanoti_paired.bam
		samtools index "$1"_BTV8_tanoti_paired.bam
		samtools sort "$1"_BTV8_tanoti_unpaired.bam -o "$1"_BTV8_tanoti_unpaired.bam
		samtools index "$1"_BTV8_tanoti_unpaired.bam
		samtools merge "$1"_BTV8_tanoti.bam "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam
		samtools sort "$1"_BTV8_tanoti.bam -o "$1"_BTV8_tanoti.bam
		samtools index "$1"_BTV8_tanoti.bam
		rm -f "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam "$1"_BTV8_tanoti_paired.bam.bai "$1"_BTV8_tanoti_unpaired.bam.bai
	
		##call SNPs
		lofreq call --no-default-filter -f ./"$3" -o BTV8a.vcf "$1"_BTV8_tanoti.bam
		lofreq filter --cov-min 1 -i BTV8a.vcf -o BTV8.vcf
		rm -f BTV8a.vcf
		
		##generate copy for bitwise checks
		cp "$3"_working.fasta temp.fasta

		##generate consensus and rebuild reference
		Rscript ../consensusinitial.R BTV8.vcf ./"$3"_working.fasta
		cat "$3"_working.fasta
		cat BTV8.vcf

		cp "$3"_working.fasta $3
	else
		if cmp -s "$3"_working.fasta temp.fasta
		then
			echo "Same."
		else
			echo "Different."
			
			##remove copy and current bam file
			rm -f temp.fasta BTV8.vcf "$1"_BTV8_tanoti.bam "$1"_BTV8_tanoti.bam.bai "$3".fai
			
			##mapping cleaned reads to BTV-8 reference genome using tanoti
			tanoti -u 0 -r $3 -i out.notCombined_1.fastq out.notCombined_2.fastq -o "$1"_BTV8_tanoti_paired.sam -p 1
			tanoti -u 0 -r $3 -i "$1".mergedsingletons.fq -o "$1"_BTV8_tanoti_unpaired.sam
			samtools view -bS "$1"_BTV8_tanoti_paired.sam > "$1"_BTV8_tanoti_paired.bam
			samtools view -bS "$1"_BTV8_tanoti_unpaired.sam > "$1"_BTV8_tanoti_unpaired.bam
			rm -f "$1"_BTV8_tanoti_paired.sam "$1"_BTV8_tanoti_unpaired.sam
			samtools sort "$1"_BTV8_tanoti_paired.bam -o "$1"_BTV8_tanoti_paired.bam
			samtools index "$1"_BTV8_tanoti_paired.bam
			samtools sort "$1"_BTV8_tanoti_unpaired.bam -o "$1"_BTV8_tanoti_unpaired.bam
			samtools index "$1"_BTV8_tanoti_unpaired.bam
			samtools merge "$1"_BTV8_tanoti.bam "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam
			samtools sort "$1"_BTV8_tanoti.bam -o "$1"_BTV8_tanoti.bam
			samtools index "$1"_BTV8_tanoti.bam
			rm -f "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam "$1"_BTV8_tanoti_paired.bam.bai "$1"_BTV8_tanoti_unpaired.bam.bai
			
			##call SNPs
			lofreq call --no-default-filter -f ./"$3" -o BTV8a.vcf "$1"_BTV8_tanoti.bam
			lofreq filter --cov-min 1 -i BTV8a.vcf -o BTV8.vcf
			rm -f BTV8a.vcf
		
			##generate copy for bitwise checks
			cp "$3"_working.fasta temp.fasta

			##generate consensus and rebuild bowtie reference
			Rscript ../consensusinitial.R BTV8.vcf ./"$3"_working.fasta
			cat "$3"_working.fasta
                	cat BTV8.vcf			

			cp "$3"_working.fasta $3
		fi
	fi
done

##remove copy and current bam file
rm -f temp.fasta BTV8.vcf "$1"_BTV8_tanoti.bam "$1"_BTV8_tanoti.bam.bai "$3".fai

##mapping cleaned reads to BTV-8 reference genome using tanoti
tanoti -u 0 -r $3 -i out.notCombined_1.fastq out.notCombined_2.fastq -o "$1"_BTV8_tanoti_paired.sam -p 1
tanoti -u 0 -r $3 -i "$1".mergedsingletons.fq -o "$1"_BTV8_tanoti_unpaired.sam
samtools view -bS "$1"_BTV8_tanoti_paired.sam > "$1"_BTV8_tanoti_paired.bam
samtools view -bS "$1"_BTV8_tanoti_unpaired.sam > "$1"_BTV8_tanoti_unpaired.bam
rm -f "$1"_BTV8_tanoti_paired.sam "$1"_BTV8_tanoti_unpaired.sam
samtools sort "$1"_BTV8_tanoti_paired.bam -o "$1"_BTV8_tanoti_paired.bam
samtools index "$1"_BTV8_tanoti_paired.bam
samtools sort "$1"_BTV8_tanoti_unpaired.bam -o "$1"_BTV8_tanoti_unpaired.bam
samtools index "$1"_BTV8_tanoti_unpaired.bam
samtools merge "$1"_BTV8_tanoti.bam "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam
samtools sort "$1"_BTV8_tanoti.bam -o "$1"_BTV8_tanoti.bam
samtools index "$1"_BTV8_tanoti.bam
rm -f "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam "$1"_BTV8_tanoti_paired.bam.bai "$1"_BTV8_tanoti_unpaired.bam.bai

##call SNPs
lofreq call --no-default-filter -f ./"$3" -o BTV8a.vcf "$1"_BTV8_tanoti.bam
lofreq filter --cov-min 1 -i BTV8a.vcf -o BTV8.vcf
rm -f BTV8a.vcf

##generate consensus via R script and via mpileup for bitwise comparison
samtools mpileup -uf "$3"_working.fasta "$1"_BTV8_tanoti.bam | bcftools call -c | vcfutils.pl vcf2fq > "$3"_mpileup.fq
seqtk seq -A "$3"_mpileup.fq > "$3"_mpileup.fasta
Rscript ../consensusfinal.R BTV8.vcf ./"$3"_working.fasta
cat "$3"_working.fasta
cat BTV8.vcf

cp "$3"_working.fasta $3

##keep final vcf
mkdir BTV8final

mv ./BTV8.vcf ./BTV8final/BTV8.vcf

##remove copy
rm -f temp.fasta "$1"_BTV8_tanoti.bam "$1"_BTV8_tanoti.bam.bai "$3"_working.fasta.fai "$3".fai "$3"_mpileup.fq

##map to original unFLASHed reads for final consensus generation
tanoti -u 0 -r $3 -i "$1"_1.fq "$1"_2.fq -o "$1"_BTV8_tanoti_paired.sam -p 1
tanoti -u 0 -r $3 -i "$1".mergedoriginalsingletons.fq -o "$1"_BTV8_tanoti_unpaired.sam
samtools view -bS "$1"_BTV8_tanoti_paired.sam > "$1"_BTV8_tanoti_paired.bam
samtools view -bS "$1"_BTV8_tanoti_unpaired.sam > "$1"_BTV8_tanoti_unpaired.bam
rm -f "$1"_BTV8_tanoti_paired.sam "$1"_BTV8_tanoti_unpaired.sam
samtools sort "$1"_BTV8_tanoti_paired.bam -o "$1"_BTV8_tanoti_paired.bam
samtools index "$1"_BTV8_tanoti_paired.bam
samtools sort "$1"_BTV8_tanoti_unpaired.bam -o "$1"_BTV8_tanoti_unpaired.bam
samtools index "$1"_BTV8_tanoti_unpaired.bam
samtools merge "$1"_BTV8_tanoti.bam "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam
samtools sort "$1"_BTV8_tanoti.bam -o "$1"_BTV8_tanoti.bam
samtools index "$1"_BTV8_tanoti.bam
rm -f "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam "$1"_BTV8_tanoti_paired.bam.bai "$1"_BTV8_tanoti_unpaired.bam.bai

##convert any areas with less than 5 reads to Ns for both R script...
samtools depth -aa "$1"_BTV8_tanoti.bam > depth.txt
Rscript ../finaldepthcorrection.R depth.txt ./"$3"_working.fasta
cp "$3"_working.fasta $3
##...and mpileup
Rscript ../finaldepthcorrection.R depth.txt ./"$3"_mpileup.fasta

rm -f "$1"_BTV8_tanoti.bam "$1"_BTV8_tanoti.bam.bai "$1"_1.fq "$1"_2.fq

##check if masked mpileup and R script consensuses are different, if different save locations that do not match
Rscript ../sequencecomparison.R ./"$3" ./"$3"_mpileup.fasta ../"$3"

##map for final statistics
tanoti -u 0 -r $3 -i out.notCombined_1.fastq out.notCombined_2.fastq -o "$1"_BTV8_tanoti_paired.sam -p 1
tanoti -u 0 -r $3 -i "$1".mergedsingletons.fq -o "$1"_BTV8_tanoti_unpaired.sam
samtools view -bS "$1"_BTV8_tanoti_paired.sam > "$1"_BTV8_tanoti_paired.bam
samtools view -bS "$1"_BTV8_tanoti_unpaired.sam > "$1"_BTV8_tanoti_unpaired.bam
rm -f "$1"_BTV8_tanoti_paired.sam "$1"_BTV8_tanoti_unpaired.sam
samtools sort "$1"_BTV8_tanoti_paired.bam -o "$1"_BTV8_tanoti_paired.bam
samtools index "$1"_BTV8_tanoti_paired.bam
samtools sort "$1"_BTV8_tanoti_unpaired.bam -o "$1"_BTV8_tanoti_unpaired.bam
samtools index "$1"_BTV8_tanoti_unpaired.bam
samtools merge "$1"_BTV8_tanoti.bam "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam
samtools sort "$1"_BTV8_tanoti.bam -o "$1"_BTV8_tanoti.bam
samtools index "$1"_BTV8_tanoti.bam
rm -f "$1"_BTV8_tanoti_paired.bam "$1"_BTV8_tanoti_unpaired.bam "$1"_BTV8_tanoti_paired.bam.bai "$1"_BTV8_tanoti_unpaired.bam.bai

##statistics and coverage plot of mapped reads BTV8
cd BTV8final
weeSAM --bam ../"$1"_BTV8_tanoti.bam --html ./BTV8final/"$1"_BTV8.html --out output.txt
cd ..
mv "$1"_BTV8_tanoti.bam ./BTV8final/"$1"_BTV8_tanoti.bam
mv "$3"_working.fasta ./BTV8final/"$1"_consensus.fasta
mv "$3"_mpileup.fasta ./BTV8final/"$1"_mpileupconsensus.fasta
mv differences.csv ./BTV8final/differences.csv
mv depth.txt ./BTV8final/depth.txt
Rscript ../namecorrection.R ./BTV8final/"$1"_consensus.fasta "$1"
##clean-up
rm -f $3 "$2".ref* "$1"_1_val_1_fastqc.html "$1"_2_val_2_fastqc.html "$1"_1.fq_trimming_report.txt "$1"_2.fq_trimming_report.txt *.bai *.fq *.fastq
