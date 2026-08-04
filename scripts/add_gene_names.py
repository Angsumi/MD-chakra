import re
import os

def map_gene_names():
    gtf_path = 'reference_genome/annotation.gtf'
    counts_path = 'counts/gene_counts.txt'
    output_path = 'counts/gene_counts_with_names.csv'
    
    if not os.path.exists(gtf_path):
        print(f"Error: {gtf_path} not found.")
        return
    if not os.path.exists(counts_path):
        print(f"Error: {counts_path} not found.")
        return
        
    print("Parsing GTF file for gene name mappings...")
    gene_map = {} # gene_id -> gene_name
    
    # Regex to find gene_id and gene name
    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene "([^"]+)"')
    
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            attributes = parts[8]
            gid_match = gene_id_re.search(attributes)
            name_match = gene_name_re.search(attributes)
            
            if gid_match:
                gid = gid_match.group(1)
                # If there's a gene name, use it; otherwise default to locus tag/Gene ID
                gname = name_match.group(1) if name_match else gid
                # Keep the first or any valid mapping found
                if gid not in gene_map or gene_map[gid] == gid:
                    gene_map[gid] = gname

    print("Parsing gene counts and merging gene names...")
    with open(counts_path, 'r') as infile, open(output_path, 'w') as outfile:
        # Write CSV header
        outfile.write("GeneID,GeneName,Chromosome,Start,End,Strand,Length,NPFM1,NPFM2,NPFM3,NPFM4\n")
        
        for line in infile:
            if line.startswith('#') or line.startswith('Geneid'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            
            gene_id = parts[0]
            chr_num = parts[1]
            start = parts[2]
            end = parts[3]
            strand = parts[4]
            length = parts[5]
            
            npfm1 = parts[6]
            npfm2 = parts[7]
            npfm3 = parts[8]
            npfm4 = parts[9]
            
            gene_name = gene_map.get(gene_id, gene_id)
            
            # Escape fields just in case of commas in gene names
            if ',' in gene_name:
                gene_name = f'"{gene_name}"'
                
            outfile.write(f"{gene_id},{gene_name},{chr_num},{start},{end},{strand},{length},{npfm1},{npfm2},{npfm3},{npfm4}\n")

    print(f"Successfully created: {output_path}")

if __name__ == '__main__':
    map_gene_names()
