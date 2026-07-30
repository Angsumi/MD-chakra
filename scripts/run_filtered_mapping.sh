#!/bin/bash
# run_filtered_mapping.sh
# Script to map rRNA-filtered reads back to the genome to get Requirement 4 stats

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

mkdir -p filtered_alignments
THREADS=4

echo "======================================"
echo "Starting HISAT2 mapping on rRNA filtered reads..."
echo "======================================"

for i in 1 2 3 4; do
    SAMPLE="NPFM${i}"
    echo "Running HISAT2 alignment for filtered $SAMPLE..."
    
    # Map the clean mRNA reads to the genome
    # We use THREADS=4 to avoid the out-of-memory error from before
    hisat2 -p $THREADS -x reference_index \
           -1 rrna_filtered/${SAMPLE}_mRNA_R1.fastq.gz \
           -2 rrna_filtered/${SAMPLE}_mRNA_R2.fastq.gz 2> filtered_alignments/${SAMPLE}_mapping_stats.txt | \
    samtools sort -@ $THREADS -o filtered_alignments/${SAMPLE}_filtered_sorted.bam
    
    samtools index filtered_alignments/${SAMPLE}_filtered_sorted.bam
done

# Run MultiQC to aggregate the new mapping stats
echo "Generating MultiQC Report for Requirement 4..."
multiqc filtered_alignments -o filtered_multiqc_report

echo "======================================"
echo "Filtered Mapping Completed Successfully!"
echo "======================================"
