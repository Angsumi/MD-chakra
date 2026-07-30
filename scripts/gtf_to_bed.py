import sys
import re

with open('annotation.gtf') as f, open('cds.bed', 'w') as out:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) > 8 and parts[2] == 'CDS':
            chrom = parts[0]
            start = int(parts[3]) - 1
            end = parts[4]
            strand = parts[6]
            
            m = re.search(r'transcript_id "([^"]+)"', parts[8])
            transcript_id = m.group(1) if m else "unknown"
            
            # Format: chrom, start, end, name, score, strand
            out.write(f"{chrom}\t{start}\t{end}\t{transcript_id}\t0\t{strand}\n")
