#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate eggnog

echo "Resuming downloads..."
cd "eggnog_db"
wget -c -O eggnog.db.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
gunzip -f eggnog.db.gz

wget -c -O eggnog.taxa.tar.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
tar -zxf eggnog.taxa.tar.gz
rm eggnog.taxa.tar.gz

wget -c -O eggnog_proteins.dmnd.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
gunzip -f eggnog_proteins.dmnd.gz
cd ..

echo "Starting emapper.py annotation..."
emapper.py -i cds.fasta -o eggnog_out --data_dir eggnog_db -m diamond --translate

echo "eggNOG-mapper pipeline completed successfully!"
