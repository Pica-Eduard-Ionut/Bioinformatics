# Implement a software application that makes DNA patterns from DNA sequences of gene promoters. Compute the C+G% and the Kappa Index of Coincidence values from each sliding window.
# To ensure that the algorithm implementation works correctly, use the following steps:
# 1.Use the sequence: S=“CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG”.
# 2. Use a sliding window of length 30b.
# 3. Build a function to process the CpG content. The value that the function must return is: CG = 29.27
# 4. Build a function to process the Index of Coincidence. The value that the function must return is: IC = 27.53
# 5. Plot the pattern on a chart.
# 6. Calculate the center of weight of the pattern.
# 7. Take the center of each pattern and plot it on a second chart.
# 8. Take the DNA sequence of a promoter and generate a pattern. Open PromKappa and see if your pattern is the same.

import numpy as np
import matplotlib.pyplot as plt

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"

# promoter
# S = "CGGAACCGGGCGCCGCGGAGCAGGGGGCGGCGCTCGTCGGGGAGGCTCTGGCGTCGCAGGAGCCCACCGCGGGGGTTGGGAGGAGGCTCGGGCATGGCGGGCTGCAGTTCCCCAGCCCTGCGCCGCGGGGAGGCAGCTGAGGCCTGGCGAGAATTCGAGCGCAGTGCCGGCACTGCTGGGGAACCCGGCGCACCCTCCGCAGCTGCTAGCCCGGGTGCTAAGCCCCTCACTGCCCGGGCCGGCCGGCGGCTCCGAGTGCGGGGCCCGCGGGGCCCACGCCCACCCGGAACTCGCGCTGGCCCGCGAGCGCCGCGCGCAGCCCAGGTTCCCGCCCGCGCCACTCCCTCCACACCTCCCCACAAGCAGAGGGAGCTGGCTCCCGCATTGGCCAGCCCAGAGAGGGGCCCCCAGAGCGCAGTGTCGGGCTGAAGGGCTTCTCGAGCATGCCCAGAGCGGACGCCGAGGCCGAGGAGGCGCTGAGAGCGAGCAAGGGCTGCGAGGGCTGCCAGCATGCTGTCACCTCTCAATAGGAGATCAGTTAATGCATACTGAAGGAAGGCTTGTTGGAGAAGAATCCTCTCCTGAACCCTGTGGAGACT"

windowLength = 30

def calculateCgContent(seq):
    cgCount = seq.count('C') + seq.count('G')
    totalLength = len(seq)
    cgPercent = (cgCount / totalLength) * 100
    return round(cgPercent, 2)


def calculateIc(sequence):
    sequence = sequence.upper()
    N = len(sequence)

    if N < 2:
        return 0.0

    totalPercentage = 0.0

    # Shifts from 1 to N-1
    for shift in range(1, N):
        matches = 0
        overlap = N - shift

        for i in range(overlap):
            if sequence[i] == sequence[i + shift]:
                matches += 1

        pct = (matches / overlap) * 100.0
        totalPercentage += pct

    icValue = totalPercentage / (N - 1)
    return round(icValue, 2)


def processSequence(seq, windowLen):
    windows = [seq[i:i+windowLen] for i in range(len(seq) - windowLen + 1)]
    cgValues = [calculateCgContent(w) for w in windows]
    icValues = [calculateIc(w) for w in windows]
    return cgValues, icValues

cgValues, icValues = processSequence(S, windowLength)

def plotPattern(cgValues, icValues):
    plt.figure(figsize=(12, 6))

    # Plot both CpG% and IC on the same graph
    plt.plot(cgValues, label='CpG Content (C+G%)', color='blue')
    plt.plot(icValues, label='Index of Coincidence (IC%)', color='green')

    plt.title('CpG Content and Index of Coincidence')
    plt.xlabel('Sliding Window Position')
    plt.ylabel('Percentage (%)')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

plotPattern(cgValues, icValues)


def calculateCenterOfWeight(values):
    weightedPositions = np.array([i * value for i, value in enumerate(values)])
    totalWeight = sum(values)
    return weightedPositions.sum() / totalWeight if totalWeight != 0 else 0

centerCg = calculateCenterOfWeight(cgValues)
centerIc = calculateCenterOfWeight(icValues)

def plotCenters(centerCg, centerIc, seqLength):
    plt.figure(figsize=(6, 4))

    plt.axvline(centerCg, color='blue',
                label=f'CpG Center: {centerCg:.2f}')
    plt.axvline(centerIc, color='green', linestyle='--',
                label=f'IC Center: {centerIc:.2f}')

    plt.title('Centers of Weight for CpG and IC')
    plt.xlabel('Position')
    plt.ylabel('Center Value')
    plt.xlim(0, seqLength - windowLength)
    plt.legend()
    plt.show()

plotCenters(centerCg, centerIc, len(S))

averageIC = round(np.mean(icValues), 2)

print("Calculated CG: ", calculateCgContent(S))
print("Calculated IC: ", averageIC)

