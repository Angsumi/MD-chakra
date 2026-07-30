#!/bin/bash

# resume_pipeline.sh
# Resumes the pipeline for RNA-Seq analysis from NPFM4 HISAT2 alignment

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

mkdir -p fastp_qc
mkdir -p alignments
mkdir -p counts

THREADS=4

echo "Resuming Pipeline for NPFM4..."

SAMPLE="NPFM4"
echo "======================================"
echo "Processing Sample: $SAMPLE"
echo "======================================"

# Skip fastp as it's already done
echo "Skipping fastp for $SAMPLE (already completed)..."

# STEP 4 & 5: Alignment & BAM Conversion (HISAT2 -> SAMtools)
echo "Running HISAT2 alignment for $SAMPLE..."
hisat2 -p $THREADS -x reference_index \
       -1 fastp_qc/${SAMPLE}_trimmed_R1.fastq.gz \
       -2 fastp_qc/${SAMPLE}_trimmed_R2.fastq.gz | \
samtools sort -@ $THREADS -o alignments/${SAMPLE}_sorted.bam

# Index the BAM file
samtools index alignments/${SAMPLE}_sorted.bam

# STEP 6: MultiQC to aggregate all fastp QC and HISAT2 logs
echo "Running MultiQC..."
multiqc fastp_qc alignments -o multiqc_report

# STEP 7: Read Counting (featureCounts)
echo "Running featureCounts on all BAM files..."
featureCounts -T $THREADS -p -a annotation.gtf -o counts/gene_counts.txt alignments/*_sorted.bam

echo "======================================"
echo "Pipeline Completed Successfully!"
echo "======================================"
