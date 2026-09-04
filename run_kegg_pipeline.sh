#!/bin/bash
set -e

echo "=========================================================="
echo "Step 1: Decompressing eggNOG Diamond database & Taxa data"
echo "=========================================================="
cd kegg_analysis

if [ ! -f "eggnog_proteins.dmnd" ]; then
    echo "Decompressing eggnog_proteins.dmnd.gz..."
    gunzip -k eggnog_proteins.dmnd.gz
fi

if [ ! -d "eggnog.taxa" ]; then
    echo "Extracting eggnog.taxa.tar.gz..."
    tar -zxf eggnog.taxa.tar.gz
fi

if [ ! -f "eggnog.db" ]; then
    echo "Decompressing eggnog.db.gz..."
    gunzip -k eggnog.db.gz
fi

echo "=========================================================="
echo "Step 2: Running Diamond Blast against eggNOG (KEGG) DB"
echo "=========================================================="
# Run Diamond blastp with high sensitivity and 12 threads
./diamond blastp \
  --db eggnog_proteins.dmnd \
  --query proteins.fasta \
  --out diamond_hits.seed_orthologs \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  --max-target-seqs 1 \
  --evalue 0.001 \
  --threads 12

echo "Diamond search complete! Parsing KEGG pathways..."
cd ..

python3 scripts/parse_kegg_annotations.py

echo "=========================================================="
echo "KEGG Pipeline Finished Successfully!"
echo "=========================================================="
