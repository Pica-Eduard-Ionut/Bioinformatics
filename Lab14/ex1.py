# For M1 we have: S1 = “ATCGATTCGATATCATACACGTAT”, known to belong to a CpG island.
# For M2 we have: S2 = “CTCGACTAGTATGAAGTCCACGCTTG”, known to belong to other regions in the genome.
#
# Follow the steps below to implement a software application:
# 1. for the CpG+ model: Count the transition frequencies from the known sequence "S1" which does belong to a CpG island.
# 2. for the CpG- model: Also count the transition frequencies from the sequence "S2" which does not belong to a CpG island.
# 3. Make the log likelihood matrix:

import math
from collections import defaultdict

S1 = "ATCGATTCGATATCATACACGTAT"
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"

bases = ["A", "C", "G", "T"]

def count_transitions(sequence):
    counts = defaultdict(lambda: defaultdict(int))
    for i in range(len(sequence) - 1):
        counts[sequence[i]][sequence[i+1]] += 1
    return counts

counts_plus = count_transitions(S1)
counts_minus = count_transitions(S2)


def transition_probabilities(counts, pseudocount=1):
    probs = defaultdict(dict)
    for b1 in bases:
        total = sum(counts[b1][b2] for b2 in bases) + 4 * pseudocount
        for b2 in bases:
            probs[b1][b2] = (counts[b1][b2] + pseudocount) / total
    return probs

p_plus = transition_probabilities(counts_plus)
p_minus = transition_probabilities(counts_minus)

log_likelihood = defaultdict(dict)
for b1 in bases:
    for b2 in bases:
        log_likelihood[b1][b2] = math.log2(p_plus[b1][b2] / p_minus[b1][b2])

print("Transition counts (CpG+ model):")
for b1 in bases:
    print(b1, dict(counts_plus[b1]))

print("\nTransition counts (CpG- model):")
for b1 in bases:
    print(b1, dict(counts_minus[b1]))

print("\nLog-likelihood matrix:")
for b1 in bases:
    print(b1, {b2: round(log_likelihood[b1][b2], 3) for b2 in bases})
