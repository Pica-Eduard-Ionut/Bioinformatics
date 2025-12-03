# download from moodle the fasta file containing promoter sequences and use it as an input for an updated application which is able to generate and save the objective digital stain (ODS) inside the folder. The promotor file can be found inside the PromKappa package on moodle or Github inside folder "bin"

# OX is CG content and OY is IC

import os
import numpy as np
import matplotlib.pyplot as plt


def read_fasta(filepath):
    sequences = {}
    header = None
    seq_chunks = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(seq_chunks)

                header = line[1:]
                seq_chunks = []

            else:
                seq_chunks.append(line)

        if header is not None:
            sequences[header] = "".join(seq_chunks)

    return sequences

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

def calculateCenterOfWeight(values):
    weightedPositions = np.array([i * value for i, value in enumerate(values)])
    totalWeight = sum(values)
    return weightedPositions.sum() / totalWeight if totalWeight != 0 else 0

def plotPattern(cgValues, icValues, output_path, seq_name):
    plt.figure(figsize=(12, 6))
    plt.plot(cgValues, label='CpG Content (C+G%)', color='blue')
    plt.plot(icValues, label='Index of Coincidence (IC%)', color='green')

    plt.title(f'CpG Content and IC - {seq_name}')
    plt.xlabel('Sliding Window Position')
    plt.ylabel('Percentage (%)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()

def plotCenters(centerCg, centerIc, seqLength, windowLength, output_path, seq_name):
    plt.figure(figsize=(6, 4))
    plt.axvline(centerCg, color='blue', label=f'CpG Center: {centerCg:.2f}')
    plt.axvline(centerIc, color='green', linestyle='--', label=f'IC Center: {centerIc:.2f}')

    plt.title(f'Centers of Weight - {seq_name}')
    plt.xlabel('Position')
    plt.xlim(0, seqLength - windowLength)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def plotODS(cgValues, icValues, output_path, seq_name):
    plt.figure(figsize=(6, 6))
    plt.scatter(cgValues, icValues, s=10, alpha=0.7, color="purple")

    plt.xlabel("CG Content (%)")
    plt.ylabel("Index of Coincidence (%)")
    plt.title(f"ODS Plot (CG vs IC) - {seq_name}")
    plt.grid(True)
    plt.savefig(output_path, dpi=200)
    plt.close()


fasta_file = "promoters.fasta"
windowLength = 30

out_folder = "ODS_output"
os.makedirs(out_folder, exist_ok=True)

promoters = read_fasta(fasta_file)
if not promoters:
    print("No sequences found in promoters.fasta.")


for seq_name, seq in promoters.items():
    seq = seq.upper()
    print(f"Processing: {seq_name} (length {len(seq)})")

    cgValues, icValues = processSequence(seq, windowLength)

    centerCg = calculateCenterOfWeight(cgValues)
    centerIc = calculateCenterOfWeight(icValues)

    plotPattern(cgValues, icValues,
                f"{out_folder}/{seq_name}_CG_IC_pattern.png",
                seq_name)

    plotCenters(centerCg, centerIc, len(seq), windowLength,
                f"{out_folder}/{seq_name}_centers.png",
                seq_name)

    plotODS(cgValues, icValues,
            f"{out_folder}/{seq_name}_ODS.png",
            seq_name)

print("\nODS generation complete! Files saved in:", out_folder)

