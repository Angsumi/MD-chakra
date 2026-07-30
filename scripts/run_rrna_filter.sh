#!/bin/bash
# run_rrna_filter.sh
# Script to install SortMeRNA, download databases, and filter rRNA

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/angsuman/extra_spac/conda_envs/sortmerna

echo "======================================"
echo "2. Downloading rRNA Reference Database..."
echo "======================================"
if [ ! -d "sortmerna_repo" ]; then
    git clone https://github.com/biocore/sortmerna.git sortmerna_repo
fi
DB_PATH="sortmerna_repo/data/rRNA_databases/smr_v4.3_default_db.fasta"

mkdir -p rrna_filtered
mkdir -p sortmerna_logs

THREADS=4

echo "======================================"
echo "3. Starting rRNA Filtration..."
echo "======================================"

for i in 1 2 3 4; do
    SAMPLE="NPFM${i}"
    echo "Filtering rRNA for $SAMPLE..."
    
    # Run SortMeRNA
    # --fastx outputs fasta/fastq
    # --other outputs non-aligned reads (mRNA)
    # --aligned outputs aligned reads (rRNA)
    # --paired_in and --out2 handle paired-end correctly
    sortmerna --ref $DB_PATH \
              --reads fastp_qc/${SAMPLE}_trimmed_R1.fastq.gz \
              --reads fastp_qc/${SAMPLE}_trimmed_R2.fastq.gz \
              --workdir sortmerna_workdir_${SAMPLE} \
              --aligned rrna_filtered/${SAMPLE}_rRNA \
              --other rrna_filtered/${SAMPLE}_mRNA \
              --fastx \
              --paired_in \
              --out2 \
              --threads $THREADS
              
    # Move log for QC report and clean up huge temp directory
    mv sortmerna_workdir_${SAMPLE}/out/aligned.log sortmerna_logs/${SAMPLE}_sortmerna.log
    rm -rf sortmerna_workdir_${SAMPLE}
    
    echo "$SAMPLE filtration complete!"
done

echo "======================================"
echo "rRNA Filtration Completed Successfully!"
echo "======================================"
