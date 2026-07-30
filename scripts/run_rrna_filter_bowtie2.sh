#!/bin/bash
# run_rrna_filter_bowtie2.sh
# Script to install Bowtie2, index the rRNA database, and filter rRNA

export CONDA_PKGS_DIRS=/home/angsuman/extra_spac/conda_pkgs
source ~/miniconda3/etc/profile.d/conda.sh

echo "======================================"
echo "1. Installing Bowtie2..."
echo "======================================"
conda create -y --prefix /home/angsuman/extra_spac/conda_envs/bowtie2 -c bioconda bowtie2
conda activate /home/angsuman/extra_spac/conda_envs/bowtie2

echo "======================================"
echo "2. Downloading and indexing rRNA database..."
echo "======================================"
if [ ! -d "database" ]; then
    wget -nc https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
    tar -xzf database.tar.gz
fi

mkdir -p rrna_index
bowtie2-build --threads 4 smr_v4.3_default_db.fasta rrna_index/rRNA_db

mkdir -p rrna_filtered
mkdir -p bowtie2_logs

THREADS=4

echo "======================================"
echo "3. Starting rRNA Filtration..."
echo "======================================"

for i in 1 2 3 4; do
    SAMPLE="NPFM${i}"
    echo "Filtering rRNA for $SAMPLE..."
    
    # Run Bowtie2 mapping against the rRNA index
    # --un-conc-gz saves the read pairs that DO NOT align (the pure mRNA)
    # -S /dev/null discards the aligned SAM file to save space
    # 2> redirects the summary statistics to a log file
    bowtie2 -p $THREADS -x rrna_index/rRNA_db \
            -1 fastp_qc/${SAMPLE}_trimmed_R1.fastq.gz \
            -2 fastp_qc/${SAMPLE}_trimmed_R2.fastq.gz \
            --un-conc-gz rrna_filtered/${SAMPLE}_mRNA_R%.fastq.gz \
            -S /dev/null 2> bowtie2_logs/${SAMPLE}_rRNA_mapping.log
            
    # Rename outputs to match expected R1/R2 naming gracefully
    mv rrna_filtered/${SAMPLE}_mRNA_R1.fastq.gz rrna_filtered/${SAMPLE}_mRNA_R1.fastq.gz
    mv rrna_filtered/${SAMPLE}_mRNA_R2.fastq.gz rrna_filtered/${SAMPLE}_mRNA_R2.fastq.gz
    
    echo "$SAMPLE filtration complete!"
done

echo "======================================"
echo "rRNA Filtration Completed Successfully!"
echo "======================================"
