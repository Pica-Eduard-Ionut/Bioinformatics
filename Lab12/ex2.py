# Download 10 influenza genomes. Adapt your application from the previous assignment in order to scan each genome for possible motifs.
# For each genome make a chart that shows the signal with the most likely locations of real functional motifs.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

motifs = [
    "GTCATTACTA",
    "ACACAATAGA",
    "GCGAGGGGTG",
    "GGGGGGGGGG",
    "TTTTTTTTTT",
    "AATCCAAAGA",
    "AAGAACATAA",
    "AGGGTTCAGG",
    "CTATTGTCTT"
]

nucleotides = ['A', 'C', 'G', 'T']
ntIndex = {nt: i for i, nt in enumerate(nucleotides)}

motifLength = len(motifs[0])

countMatrix = np.zeros((4, motifLength))
for motif in motifs:
    for position, nucleotide in enumerate(motif):
        countMatrix[ntIndex[nucleotide], position] += 1

pseudoCount = 1
countMatrixPc = countMatrix + pseudoCount
frequencyMatrix = countMatrixPc / countMatrixPc.sum(axis=0)
nullModel = 0.25
logLikelihoodMatrix = np.log(frequencyMatrix / nullModel)

def readFasta(filePath):
    seq = ""
    with open(filePath, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq += line.strip().upper()
    return seq

def scanSequence(sequence, pwmMatrix, motifLength):
    scores = []
    for i in range(len(sequence) - motifLength + 1):
        window = sequence[i:i + motifLength]
        score = 0
        for pos, nt in enumerate(window):
            if nt in ntIndex:
                score += pwmMatrix[ntIndex[nt], pos]
        scores.append(score)
    return scores

def getTopPeaks(scores, topN=5):
    scoresArr = np.array(scores)

    topIndices = scoresArr.argsort()[-topN:][::-1]
    return topIndices, scoresArr[topIndices]

fig, axes = plt.subplots(5, 2, figsize=(15, 20))
axes = axes.flatten()

for i in range(1, 11):
    fastaFile = f"Influenza{i}.fasta"
    if not os.path.exists(fastaFile):
        print(f"File {fastaFile} not found, skipping...")
        continue

    genomeSeq = readFasta(fastaFile)
    scores = scanSequence(genomeSeq, logLikelihoodMatrix, motifLength)
    topIndices, topScores = getTopPeaks(scores, topN=5)

    ax = axes[i - 1]
    ax.plot(scores, color='blue', label='PWM Score')
    ax.scatter(topIndices, topScores, color='red', s=50, label='Top Motifs')
    ax.set_title(f"{fastaFile}")
    ax.set_xlabel("Position")
    ax.set_ylabel("Log-Likelihood Score")
    ax.grid(True)
    ax.legend()

plt.tight_layout()
plt.show()
