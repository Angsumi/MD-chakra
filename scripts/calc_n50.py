import sys

lengths = []
total_length = 0
with open('Referance.fna') as f:
    current_len = 0
    for line in f:
        if line.startswith('>'):
            if current_len > 0:
                lengths.append(current_len)
                total_length += current_len
            current_len = 0
        else:
            current_len += len(line.strip())
    if current_len > 0:
        lengths.append(current_len)
        total_length += current_len

lengths.sort(reverse=True)
n50 = 0
cum_sum = 0
for length in lengths:
    cum_sum += length
    if cum_sum >= total_length / 2:
        n50 = length
        break

print(f"Total sequences: {len(lengths)}")
print(f"Total length: {total_length}")
print(f"N50 length: {n50}")
print(f"Mean length: {total_length / len(lengths):.2f}")
