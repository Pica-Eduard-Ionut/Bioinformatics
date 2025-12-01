# Download from NCBI or other sources 10 influenza genomes. Use the restriction enzyme ECO-R1 to
# digest each genome and plot the electrophoresis Gel simulation for each genome.

# 10 genomes from the same virus (but they differ a bit)

import os
import matplotlib.pyplot as plt
import numpy as np
import random

def loadFasta(path):
    with open(path, "r") as f:
        lines = f.readlines()

    return "".join([l.strip() for l in lines if not l.startswith(">")])

genomes = []
for i in range(1, 11):
    fileName = f"Influenza{i}.fasta"
    if os.path.exists(fileName):
        genomes.append(loadFasta(fileName))

    else:
        print(f"Warning: {fileName} not found")

ecoriSite = "GAATTC"

def digestSequence(seq, site=ecoriSite):
    parts = seq.split(site)
    fragments = []
    for i, p in enumerate(parts):
        length = len(p)
        if i < len(parts) - 1:
            length += len(site)
        fragments.append(length)
    return fragments

allFragmentSets = [digestSequence(g) for g in genomes]

def lanePositions(fragments, top=50, bottom=750):
    sortedFrags = sorted(fragments)
    minBp = sortedFrags[0]
    maxBp = sortedFrags[-1]
    span = maxBp - minBp if maxBp != minBp else 1

    positions = []

    for bp in fragments:
        norm = (bp - minBp) / span
        y = bottom - norm * (bottom - top)

        y += random.uniform(-3, 3)

        positions.append((bp, y))

    return positions

plt.figure(figsize=(16, 10))
ax = plt.gca()
ax.set_facecolor("black")

laneSpacing = 2.5
bandWidth = 0.9

for laneIndex, fragments in enumerate(allFragmentSets):
    laneX = 2 + laneIndex * laneSpacing

    yPositions = lanePositions(fragments)

    for bp, y in yPositions:
        plt.hlines(
            y,
            laneX - bandWidth / 2,
            laneX + bandWidth / 2,
            color="white",
            linewidth=4
        )

    plt.text(
        laneX,
        30,
        f"Influenza {laneIndex + 1}",
        color="white",
        ha="center",
        fontsize=12
    )

plt.xlim(0, laneSpacing * len(allFragmentSets) + 3)
plt.ylim(0, 800)
plt.gca().invert_yaxis()
plt.xticks([])
plt.yticks([])
plt.title(
    "EcoRI Digest — Influenza Genomes (Per-Genome Scaled Bands)",
    color="white",
    fontsize=18
)

plt.show()

