#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

echo "Re-running featureCounts with -t gene to capture protein-coding genes..."
featureCounts -T 8 -p -t gene -g gene_id -a annotation.gtf -o counts/gene_counts.txt alignments/*_sorted.bam

echo "Re-running DESeq2..."
Rscript run_deseq2.R
