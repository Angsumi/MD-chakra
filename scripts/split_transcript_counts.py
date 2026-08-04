import os

def split_transcript_counts():
    input_file = 'counts/temp_counts_transcript.txt'
    output_dir = 'counts'
    
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found!")
        return

    # Map output file handles
    samples = ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']
    out_files = {}
    
    for sample in samples:
        out_path = os.path.join(output_dir, f"{sample}_transcript_counts.txt")
        out_files[sample] = open(out_path, 'w')
        # Write header
        out_files[sample].write("TranscriptID\tCount\n")

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Geneid'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            
            transcript_id = parts[0]
            # Columns 6, 7, 8, 9 are NPFM1, NPFM2, NPFM3, NPFM4
            counts = {
                'NPFM1': parts[6],
                'NPFM2': parts[7],
                'NPFM3': parts[8],
                'NPFM4': parts[9]
            }
            
            for sample in samples:
                out_files[sample].write(f"{transcript_id}\t{counts[sample]}\n")
                
    for sample in samples:
        out_files[sample].close()
        print(f"Created: {os.path.join(output_dir, f'{sample}_transcript_counts.txt')}")

if __name__ == '__main__':
    split_transcript_counts()
