#!/bin/bash

# run_pipeline.sh
# Automated pipeline for RNA-Seq analysis using HISAT2

# Activate Conda Environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

# Directory setup
mkdir -p fastp_qc
mkdir -p alignments
mkdir -p counts

# Number of threads to use
THREADS=12

echo "Starting Pipeline..."

# Loop over each of the 4 samples
for i in 1 2 3 4; do
    SAMPLE="NPFM${i}"
    echo "======================================"
    echo "Processing Sample: $SAMPLE"
    echo "======================================"
    
    # STEP 3: Quality Control & Trimming (fastp)
    echo "Running fastp for $SAMPLE..."
    fastp -i ${SAMPLE}_R1.fastq.gz -I ${SAMPLE}_R2.fastq.gz \
          -o fastp_qc/${SAMPLE}_trimmed_R1.fastq.gz -O fastp_qc/${SAMPLE}_trimmed_R2.fastq.gz \
          --thread $THREADS \
          --html fastp_qc/${SAMPLE}_fastp_report.html \
          --json fastp_qc/${SAMPLE}_fastp_report.json
          
    # STEP 4 & 5: Alignment & BAM Conversion (HISAT2 -> SAMtools)
    echo "Running HISAT2 alignment for $SAMPLE..."
    hisat2 -p $THREADS -x reference_index \
           -1 fastp_qc/${SAMPLE}_trimmed_R1.fastq.gz \
           -2 fastp_qc/${SAMPLE}_trimmed_R2.fastq.gz | \
    samtools sort -@ $THREADS -o alignments/${SAMPLE}_sorted.bam
    
    # Index the BAM file
    samtools index alignments/${SAMPLE}_sorted.bam
    
done

# STEP 6: MultiQC to aggregate all fastp QC and HISAT2 logs
echo "Running MultiQC..."
multiqc fastp_qc alignments -o multiqc_report

# STEP 7: Read Counting (featureCounts)
echo "Running featureCounts on all BAM files..."
featureCounts -T $THREADS -p -a annotation.gtf -o counts/gene_counts.txt alignments/*_sorted.bam

echo "======================================"
echo "Pipeline Completed Successfully!"
echo "======================================"
