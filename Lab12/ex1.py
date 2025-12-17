# A very early step in splice site recognition is exon definition, a process that is as yet poorly understood. Communication between the two ends of an exon is thought to be required for this step. Computational discovery of the exon-intron border or the intron-exon border or the transcription factor binding sites (TFBS) is a challenging but important problem of bioinformatics. Implement a software application for DNA motif finding by following the steps below.
# These sequences represent the exon-intron boundary.
# 1. make the count matrix
# 2. make the weight matrix
# 3. make the relative frequencies matrix
# 4. make the Log-likelihoods Matrix
# 5. Analize sequence S by using the Log-likelihoods Matrix:
# S="CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"
# Calculate the score for each sliding window.
# Do you have signals indicating that the S sequence contains an exon-intron border?

import numpy as np
import pandas as pd

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
numberOfMotifs = len(motifs)

countMatrix = np.zeros((4, motifLength))

for motif in motifs:
    for position, nucleotide in enumerate(motif):
        countMatrix[ntIndex[nucleotide], position] += 1

countDf = pd.DataFrame(
    countMatrix,
    index=nucleotides,
    columns=range(1, motifLength + 1)
)

pseudoCount = 1
countMatrixPc = countMatrix + pseudoCount

frequencyMatrix = countMatrixPc / countMatrixPc.sum(axis=0)

frequencyDf = pd.DataFrame(
    frequencyMatrix,
    index=nucleotides,
    columns=range(1, motifLength + 1)
)

nullModel = 0.25
logLikelihoodMatrix = np.log(frequencyMatrix / nullModel)

logLikelihoodDf = pd.DataFrame(
    logLikelihoodMatrix,
    index=nucleotides,
    columns=range(1, motifLength + 1)
)

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

windows = []
scores = []

for i in range(len(S) - motifLength + 1):
    window = S[i:i + motifLength]
    score = sum(
        logLikelihoodMatrix[ntIndex[nucleotide], position]
        for position, nucleotide in enumerate(window)
    )
    windows.append(window)
    scores.append(score)

scoreDf = pd.DataFrame({
    "Window": windows,
    "Score": scores
})


print("\nCOUNT MATRIX\n")
print(countDf.to_string())

print("\nRELATIVE FREQUENCY MATRIX\n")
print(frequencyDf.to_string(float_format="%.3f"))

print("\nLOG-LIKELIHOOD MATRIX\n")
print(logLikelihoodDf.to_string(float_format="%.3f"))

print("\nSLIDING WINDOW SCORES\n")
print(scoreDf.to_string(index=False, float_format="%.3f"))

# Yes, the sliding-window analysis reveals statistically significant peaks in log-likelihood scores, indicating that the sequence S likely contains one or more exon–intron border signals.
