import sys

def calculate_n50(lengths):
    if not lengths:
        return 0
    lengths.sort(reverse=True)
    total_length = sum(lengths)
    cum_sum = 0
    for length in lengths:
        cum_sum += length
        if cum_sum >= total_length / 2:
            return length
    return 0

with open('counts/temp_counts_gene.txt') as f:
    lines = f.readlines()

# find header
header = None
for line in lines:
    if line.startswith('Geneid'):
        header = line.strip().split('\t')
        break

sample_indices = {}
for i, col in enumerate(header):
    if 'NPFM' in col:
        # Extract sample name, e.g., 'alignments/NPFM1_sorted.bam' -> 'NPFM1'
        sample_name = col.split('/')[-1].split('_')[0]
        sample_indices[sample_name] = i

sample_lengths = {sample: [] for sample in sample_indices}

for line in lines:
    if line.startswith('#') or line.startswith('Geneid'):
        continue
    parts = line.strip().split('\t')
    if len(parts) < max(sample_indices.values()) + 1:
        continue
    length = int(parts[5])
    for sample, idx in sample_indices.items():
        count = int(parts[idx])
        if count >= 1:
            sample_lengths[sample].append(length)

print("Sample\tTotal unigenes (expressed)\tTotal length\tN50 length\tMean length")
for sample in sorted(sample_lengths.keys()):
    lengths = sample_lengths[sample]
    total_unigenes = len(lengths)
    total_length = sum(lengths)
    n50 = calculate_n50(lengths)
    mean_len = total_length / total_unigenes if total_unigenes > 0 else 0
    print(f"{sample}\t{total_unigenes}\t{total_length}\t{n50}\t{mean_len:.2f}")

