#!/bin/bash

main="/home"
###
fasta="$main/TBwiz/TB-wiz-main/Genome/GCF_000195955.2_ASM19595v2_genomic.fna"
Sample="$main/Sample"
mkdir $main/Result
mkdir $main/Result/Analysis
threads="3"


for file in "$main/Sample/"*1.fastq.gz
do
file_basename=$(basename "$file")
echo $file
file_prefix="${file_basename%*1.fastq.gz}" 
file2="${file_prefix}2.fastq.gz"
mkdir $main/Result/Analysis/$file_prefix
done

mkdir $main/Result/Analysis/Trimming
mkdir $main/Result/Analysis/Trimming/Unpaired
for file in "$main/Sample/"*1.fastq.gz
do
file_basename=$(basename "$file")
file_prefix="${file_basename%*1.fastq.gz}" 
file2="${file_prefix}2.fastq.gz"
java -jar $main/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads $threads -trimlog $main/Result/Analysis/Trimming/$file_prefix.txt $file $main/Sample/$file2 $main/Result/Analysis/Trimming/${file_prefix}1_paired.fastq $main/Result/Analysis/Trimming/${file_prefix}1_unpaired.fastq $main/Result/Analysis/Trimming/${file_prefix}2_paired.fastq $main/Result/Analysis/Trimming/${file_prefix}2_unpaired.fastq SLIDINGWINDOW:20:30 MINLEN:20
done

for file in  "$main/Result/Analysis/Trimming/"*unpaired.fastq; do
 mv "$file"  $main/Result/Analysis/Trimming/Unpaired/
done

mkdir $main/Result/Analysis/Trimming/FastQC
mkdir $main/Result/Analysis/Trimming/MultiQC
fastqc $main/Result/Analysis/Trimming/*.fastq -o $main/Result/Analysis/Trimming/FastQC
multiqc $main/Result/Analysis/Trimming/FastQC/* -o $main/Result/Analysis/Trimming/MultiQC

mkdir $main/Result/Analysis/SAM
for file in "$main/Result/Analysis/Trimming/"*1_paired.fastq
do
file_basename=$(basename "$file")
file_prefix="${file_basename%*1_paired.fastq}" 
file2="${file_prefix}2_paired.fastq"
bwa mem -t $threads $main/TBwiz/TB-wiz-main/Genome/TBC $file $main/Result/Analysis/Trimming/$file2 > $main/Result/Analysis/SAM/$file_prefix.sam
done

 
mkdir $main/Result/Analysis/BAM
mkdir $main/Result/Analysis/BAM/BAM
mkdir $main/Result/Analysis/BAM/noduplicates
for file in "$main/Result/Analysis/SAM/"*sam
do
file_basename=$(basename "$file")
file_prefix="${file_basename%.sam}" 
samtools view -@ $threads -S -b -T $fasta $main/Result/Analysis/SAM/$file_prefix.sam > $main/Result/Analysis/BAM/BAM/$file_prefix.bam
done


mkdir $main/Result/Analysis/BAM/BAM_Sorted
for file in "$main/Result/Analysis/BAM/BAM/"*bam
do
echo $file
file_basename=$(basename "$file")
file_prefix="${file_basename%.bam}" 
samtools sort -@ $threads -o $main/Result/Analysis/BAM/BAM_Sorted/$file_prefix.sorted.bam $main/Result/Analysis/BAM/BAM/$file_prefix.bam
samtools index $main/Result/Analysis/BAM/BAM_Sorted/$file_prefix.sorted.bam
java -jar $main/Tools/picard.jar MarkDuplicates -I $main/Result/Analysis/BAM/BAM_Sorted/$file_prefix.sorted.bam -O $main/Result/Analysis/BAM/noduplicates/$file_prefix.nodup.bam -METRICS_FILE $main/Result/Analysis/BAM/noduplicates/${file_prefix}_metrics.txt
samtools index $main/Result/Analysis/BAM/noduplicates/$file_prefix.nodup.bam
mkdir $main/Result/$file_prefix
$main/Tools/qualimap_v2.3/qualimap bamqc -bam $main/Result/Analysis/BAM/noduplicates/$file_prefix.nodup.bam -outdir $main/Result/$file_prefix/Qualimap
done

string=""
string2=""
mkdir $main/Result/Analysis/VCF
mkdir $main/Result/Analysis/Mpileup
for file in $(ls "$main/Result/Analysis/BAM/noduplicates/"*nodup.bam)
do 
string+="${file} "
done  
for file in "$main/Result/Analysis/BAM/noduplicates/"*nodup.bam; do
    file_basename=$(basename "$file")
    file_prefix="${file_basename%.nodup.bam}"
    printf "%s\n" "$file_prefix" >> "$main/Result/Analysis/Mpileup/name.txt"
done

samtools mpileup -B -f $fasta $string > $main/Result/Analysis/Mpileup/sample.mpileup 

java -jar $main/Tools/VarScan.v2.3.9.jar mpileup2cns  $main/Result/Analysis/Mpileup/sample.mpileup  --variants --output-vcf > $main/Result/Analysis/VCF/sample.vcf

bcftools reheader -s $main/Result/Analysis/Mpileup/name.txt -o $main/Result/Analysis/VCF/sample1.vcf $main/Result/Analysis/VCF/sample.vcf 
bcftools reheader -f $main/TBwiz/TB-wiz-main/Genome/GCF_000195955.2_ASM19595v2_genomic.fna.fai -o $main/Result/Analysis/VCF/sample2.vcf $main/Result/Analysis/VCF/sample1.vcf
MULTISAMPLEVCF="$main/Result/Analysis/VCF/sample2.vcf"
for sample in $(bcftools query -l "$MULTISAMPLEVCF"); do
    echo "$sample"
    bcftools view -c 1 -s $sample $MULTISAMPLEVCF -Oz -o "$main/Result/Analysis/$sample/$sample.vcf"
    fast-lineage-caller $main/Result/Analysis/$sample/$sample.vcf > $main/Result/$sample/$sample.lineage.txt
    bcftools annotate --rename-chrs $main/TBwiz/TB-wiz-main/tb/annotate $main/Result/Analysis/$sample/$sample.vcf > $main/Result/Analysis/$sample/$sample.corrected.vcf
    java -jar  $main/Tools/snpEff/snpEff.jar -s $main/Result/Analysis/$sample/$sample.txt -v Mycobacterium_tuberculosis_h37rv $main/Result/Analysis/$sample/$sample.corrected.vcf > $main/Result/Analysis/$sample/$sample.annotated.vcf
    /usr/bin/ruby $main/TBwiz/TB-wiz-main/tb/resistancetb.sh $main/Result/Analysis/$sample/$sample.annotated.vcf 
    cp $main/Result/Analysis/$sample/$sample.annotated.vcf.resistance.txt $main/Result/$sample/$sample.resistance.txt
done



