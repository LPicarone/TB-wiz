#!/bin/bash

set -e

# Function to check if the given path exists and has at least two FASTQ files
check_sample_path() {
    read -p "Enter the path to the Sample directory: " sample_path
    if [ ! -d "$sample_path" ]; then
        echo "Error: The provided path does not exist or is not a directory."
        exit 1
    fi

    # Count the number of FASTQ files in the directory
    fastq_count=$(find "$sample_path" -maxdepth 1 -type f -name "*.fastq.gz" | wc -l)

    # Check if there are at least two FASTQ files
    if [ "$fastq_count" -lt 2 ]; then
        echo "Error: The provided directory does not contain at least two FASTQ files."
        exit 1
    fi

    Sample="$sample_path"
}

check_main_path() {
    read -p "Enter the path to the main directory: " main_path
    if [ ! -d "$main_path" ]; then
        echo "Error: The provided path does not exist or is not a directory."
        exit 1
    fi

    main="$main_path"
}


# Function to create necessary directories
create_directories() {
  xargs -I {} mkdir -p "$main/{}" < "$directories_file"


}

# Function for Trimming
perform_trimming() {
    for file in "$Sample/"*1.fastq.gz; do
        file_basename=$(basename "$file")
        echo "$file"
        file_prefix="${file_basename%*1.fastq.gz}" 
        file2="${file_prefix}2.fastq.gz"
        mkdir -p "$main/Result/Analysis/$file_prefix"
        # Add error handling for the trimming process
        if ! java -jar /home/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads "$threads" -trimlog "$main/Result/Analysis/Trimming/$file_prefix.txt" "$file" "$Sample/$file2" "$main/Result/Analysis/Trimming/${file_prefix}1_paired.fastq" "$main/Result/Analysis/Trimming/${file_prefix}1_unpaired.fastq" "$main/Result/Analysis/Trimming/${file_prefix}2_paired.fastq" "$main/Result/Analysis/Trimming/${file_prefix}2_unpaired.fastq" SLIDINGWINDOW:20:30 MINLEN:20; then
            echo "Error occurred during trimming for file: $file"
            exit 1
        fi
    done
}

# Function to move unpaired files
move_unpaired_files() {
    mkdir -p "$main/Result/Analysis/Trimming/Unpaired"
    for file in  "$main/Result/Analysis/Trimming/"*unpaired.fastq; do
 mv "$file"  $main/Result/Analysis/Trimming/Unpaired/
done 
}

# Function to run FastQC and MultiQC
run_quality_control() {
    fastqc $main/Result/Analysis/Trimming/*.fastq -o $main/Result/Analysis/Trimming/FastQC || echo "Error occurred during FastQC"
    multiqc $main/Result/Analysis/Trimming/FastQC/* -o $main/Result/Analysis/Trimming/MultiQC || echo "Error occurred during MultiQC"
}

# Function for SAM alignment
perform_sam_alignment() {
    for file in "$main/Result/Analysis/Trimming/"*1_paired.fastq; do
        file_basename=$(basename "$file")
        file_prefix="${file_basename%*1_paired.fastq}" 
        file2="${file_prefix}2_paired.fastq"
        bwa mem -t "$threads" /home/TBwiz/TB-wiz-main/Genome/TBC "$file" "$main/Result/Analysis/Trimming/$file2" > "$main/Result/Analysis/SAM/$file_prefix.sam" || echo "Error occurred during SAM alignment for file: $file"
    done
}

# Function for BAM processing
process_bam() {
    for file in "$main/Result/Analysis/SAM/"*sam; do
        file_basename=$(basename "$file")
        file_prefix="${file_basename%.sam}" 
        samtools view -@ "$threads" -S -b -T "$fasta" "$main/Result/Analysis/SAM/$file_prefix.sam" > "$main/Result/Analysis/BAM/BAM/$file_prefix.bam" || echo "Error occurred during BAM processing for file: $file"
    done &&

    mkdir -p "$main/Result/Analysis/BAM/BAM_Sorted"
    for file in "$main/Result/Analysis/BAM/BAM/"*bam; do
        file_basename=$(basename "$file")
        file_prefix="${file_basename%.bam}" 
        samtools sort -@ "$threads" -o "$main/Result/Analysis/BAM/BAM_Sorted/$file_prefix.sorted.bam" "$main/Result/Analysis/BAM/BAM/$file_prefix.bam" || echo "Error occurred during BAM sorting for file: $file"
        samtools index "$main/Result/Analysis/BAM/BAM_Sorted/$file_prefix.sorted.bam" || echo "Error occurred while indexing BAM: $file"
        java -jar /root/Scaricati/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar MarkDuplicates -I "$main/Result/Analysis/BAM/BAM_Sorted/$file_prefix.sorted.bam" -O "$main/Result/Analysis/BAM/noduplicates/$file_prefix.nodup.bam" -METRICS_FILE "$main/Result/Analysis/BAM/noduplicates/${file_prefix}_metrics.txt" || echo "Error occurred during BAM deduplication for file: $file"
        samtools index "$main/Result/Analysis/BAM/noduplicates/$file_prefix.nodup.bam" || echo "Error occurred while indexing deduplicated BAM: $file"
        mkdir -p "$main/Result/$file_prefix"
        /home/Tools/qualimap_v2.3/qualimap bamqc -bam "$main/Result/Analysis/BAM/noduplicates/$file_prefix.nodup.bam" -outdir "$main/Result/$file_prefix/Qualimap" || echo "Error occurred during Qualimap analysis for file: $file"
    done
}

# Function for VCF generation
generate_vcf() {
    string=""
    mkdir -p "$main/Result/Analysis/VCF"
    mkdir -p "$main/Result/Analysis/Mpileup"

     name_file="$main/Result/Analysis/Mpileup/name.txt"

    # Clear the content of the existing 'name.txt' file or create a new one if it doesn't exist
    > "$name_file"

    for file in "$main/Result/Analysis/BAM/noduplicates/"*nodup.bam; do
        file_basename=$(basename "$file")
        file_prefix="${file_basename%.nodup.bam}"
        echo "$file_prefix" >> "$name_file"
        string+="${file} "
    done

    samtools mpileup -B -f "$fasta" $string > "$main/Result/Analysis/Mpileup/sample.mpileup" || echo "Error occurred during Mpileup"

    java -jar /home/Tools/VarScan.v2.3.9.jar mpileup2cns "$main/Result/Analysis/Mpileup/sample.mpileup" --variants --output-vcf --min-coverage 20 --min-reads2 10 --min-avg-qual 30 > "$main/Result/Analysis/VCF/sample.vcf" || echo "Error occurred during VCF generation"

    bcftools reheader -s "$main/Result/Analysis/Mpileup/name.txt" -o "$main/Result/Analysis/VCF/sample1.vcf" "$main/Result/Analysis/VCF/sample.vcf" || echo "Error occurred during VCF reheader"
    bcftools reheader -f "/home/TBwiz/TB-wiz-main/Genome/GCF_000195955.2_ASM19595v2_genomic.fna.fai" -o "$main/Result/Analysis/VCF/sample2.vcf" "$main/Result/Analysis/VCF/sample1.vcf" || echo "Error occurred during VCF reheader"

    MULTISAMPLEVCF="$main/Result/Analysis/VCF/sample2.vcf"
    for sample in $(bcftools query -l "$MULTISAMPLEVCF"); do
        echo "$sample"
        bcftools view -c 1 -s "$sample" "$MULTISAMPLEVCF" -Oz -o "$main/Result/Analysis/$sample/$sample.vcf" || echo "Error occurred during VCF subsetting for sample: $sample"
        fast-lineage-caller "$main/Result/Analysis/$sample/$sample.vcf" > "$main/Result/$sample/$sample.lineage.txt" || echo "Error occurred during lineage calling for sample: $sample"
        bcftools annotate --rename-chrs /home/TBwiz/TB-wiz-main/tb/annotate "$main/Result/Analysis/$sample/$sample.vcf" > "$main/Result/Analysis/$sample/$sample.corrected.vcf" || echo "Error occurred during VCF annotation for sample: $sample"
        java -jar /home/Tools/snpEff/snpEff.jar -s "$main/Result/Analysis/$sample/$sample.txt" -v Mycobacterium_tuberculosis_h37rv "$main/Result/Analysis/$sample/$sample.corrected.vcf" > "$main/Result/Analysis/$sample/$sample.annotated.vcf" || echo "Error occurred during VCF annotation with snpEff for sample: $sample"
        /usr/bin/ruby /home/TBwiz/TB-wiz-main/tb/resistancetb.sh "$main/Result/Analysis/$sample/$sample.annotated.vcf" || echo "Error occurred during resistance prediction for sample: $sample"
        cp "$main/Result/Analysis/$sample/$sample.annotated.vcf.resistance.txt" "$main/Result/$sample/$sample.resistance.txt" || echo "Error occurred while copying resistance file for sample: $sample"
    done
}

# Ask the user to input the value of threads

threads=$(nproc)
echo "Number of CPU cores/threads: $threads"
max_threads="$threads"
echo "Advisory: It is recommended to use 1-$max_threads threads based on your system's capabilities."
read -p "Enter the number of threads to use (1-$max_threads, default=$max_threads): " threads

# If the user doesn't provide any input, use the default value
if [[ -z "$threads" ]]; then
    threads=$max_threads
fi

# Check if the provided value is within the allowed range
if (( threads < 1 || threads > max_threads )); then
    echo "Error: Please provide a valid number of threads (1-$max_threads)."
    exit 1
fi





# Check the main path and ensure it exists
check_main_path

# Check the Sample path and ensure it exists with at least two FASTQ files
check_sample_path

fasta="/home/TBwiz/TB-wiz-main/Genome/GCF_000195955.2_ASM19595v2_genomic.fna"
directories_file="/home/TBwiz/TB-wiz-main/Genome/directories.txt"
# Call the functions
create_directories
perform_trimming
move_unpaired_files
run_quality_control
perform_sam_alignment
process_bam
generate_vcf

echo "Pipeline completed successfully."
