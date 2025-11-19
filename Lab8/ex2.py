# 1. Make an artificial DNA sequence of 200-400 bases in length, in which to simulate 3-4 transposable elements.
# 2. Implement a software application to detect the positions of these transposable elements (start, end) within the created DNA sequence.
# 3. Downlaod from NCBI 3 bacterial genomes of your choosing. Try to find in these genomes possible transposibles for this one must detect posible inverted repeats without prior knowledge about their existance in the
# sequence. The inverted repeat should have a min length of 4 letters and a max of 6 letters
# 4. Make a rapot about the results from #3. Write this raport in a .txt or a .docx file. Upload this raport to moodle with the proj.

import os

def rev_comp(seq):
    comp = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    return "".join(comp.get(b, "N") for b in seq[::-1])


def load_fasta(filepath):
    seq = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


def find_inverted_repeats(sequence, min_len=4, max_len=6):
    results = []
    n = len(sequence)

    for L in range(min_len, max_len + 1):
        for i in range(n - 2*L + 1):
            left = sequence[i:i+L]
            right = sequence[i+L:i+2*L]

            if right == rev_comp(left):
                results.append({
                    "left_seq": left,
                    "left_start": i + 1,
                    "left_end": i + L,
                    "right_seq": right,
                    "right_start": i + L + 1,
                    "right_end": i + 2*L
                })
    return results

for x in range(1, 4):
    fasta_file = f"Bacteria{x}.fasta"
    print(f"\n=== Processing {fasta_file} ===")

    if not os.path.exists(fasta_file):
        print(f"File not found: {fasta_file}")
        continue

    genome = load_fasta(fasta_file)
    repeats = find_inverted_repeats(genome)

    print(f"Found {len(repeats)} inverted repeats (4–6 bp).")

    for r in repeats[:20]:  # print first 20 repeats
        print(f" IR: {r['left_seq']} ... {r['right_seq']}  "
              f"({r['left_start']}-{r['left_end']} , {r['right_start']}-{r['right_end']})")
