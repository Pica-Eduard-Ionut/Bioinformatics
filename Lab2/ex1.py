from itertools import product

S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

def percentages(sequence, k):
    alphabet = sorted(set(sequence))
    size = len(sequence) - k + 1
    results = {}

    for combo in map(''.join, product(alphabet, repeat=k)):
        count = sum(sequence[i:i+k] == combo for i in range(size))
        percentage = (count / size) * 100
        results[combo] = percentage

    return results

dinucleotides = percentages(S, 2)
trinucleotide = percentages(S, 3)

print("Dinucleotide percentages:")
for k, v in dinucleotides.items():
    print(f"{k}: {v:.2f}%")

print("\nTrinucleotide percentages:")
for k, v in trinucleotide.items():
    print(f"{k}: {v:.2f}%")
