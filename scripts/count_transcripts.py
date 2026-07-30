import sys
import re

transcripts = set()
with open('annotation.gtf') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) > 8:
            # We only care about CDS lines or lines that have transcript_id with an actual ID
            m = re.search(r'transcript_id "([^"]+)"', parts[8])
            if m:
                tid = m.group(1)
                if tid:  # exclude empty ones like transcript_id ""
                    transcripts.add(tid)

print(f"Total unique transcripts in GTF: {len(transcripts)}")
