#!/bin/bash
set -e

echo "=========================================================="
echo "1. Downloading Nephila pilipes protein, CDS, and GTF sequences"
echo "=========================================================="
FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/974/015/GCA_019974015.1_Npil_1.0"

aria2c -c -x 8 -s 8 -d "reference_spider_data" \
  "${FTP_BASE}/GCA_019974015.1_Npil_1.0_protein.faa.gz" \
  "${FTP_BASE}/GCA_019974015.1_Npil_1.0_cds_from_genomic.fna.gz" \
  "${FTP_BASE}/GCA_019974015.1_Npil_1.0_rna_from_genomic.fna.gz" \
  "${FTP_BASE}/GCA_019974015.1_Npil_1.0_genomic.gtf.gz"

echo "=========================================================="
echo "2. Downloading eggNOG Taxa Database"
echo "=========================================================="
aria2c -c -x 8 -s 8 -d "Kegg databasse" \
  "http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz"

echo "All required sequences and taxa files downloaded successfully!"
