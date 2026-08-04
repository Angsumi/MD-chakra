import os

def calculate_totals():
    gene_file = 'counts/gene_counts.txt'
    transcript_file = 'counts/temp_counts_transcript.txt'
    
    gene_totals = [0, 0, 0, 0]
    transcript_totals = [0, 0, 0, 0]
    
    # Read gene counts
    with open(gene_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Geneid'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                for i in range(4):
                    gene_totals[i] += int(parts[6 + i])
                    
    # Read transcript counts
    with open(transcript_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Geneid'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                for i in range(4):
                    transcript_totals[i] += int(parts[6 + i])
                    
    # Write to CSV
    output_file = 'counts/total_counts_summary.csv'
    with open(output_file, 'w') as f:
        f.write("Feature,NPFM1,NPFM2,NPFM3,NPFM4\n")
        f.write(f"Total Transcript,{transcript_totals[0]},{transcript_totals[1]},{transcript_totals[2]},{transcript_totals[3]}\n")
        f.write(f"Total Gene,{gene_totals[0]},{gene_totals[1]},{gene_totals[2]},{gene_totals[3]}\n")
        
    print(f"Created {output_file}")
    print("Feature,NPFM1,NPFM2,NPFM3,NPFM4")
    print(f"Total Transcript,{transcript_totals[0]},{transcript_totals[1]},{transcript_totals[2]},{transcript_totals[3]}")
    print(f"Total Gene,{gene_totals[0]},{gene_totals[1]},{gene_totals[2]},{gene_totals[3]}")

if __name__ == '__main__':
    calculate_totals()
