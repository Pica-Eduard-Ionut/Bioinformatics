# 1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).
# 2. Implement a software application that detects repetitions (between 6b and 10b) in this DNA sequence.
# 3. plot the frequencies of the repetitions found
# 4. download 10 influenza genomes; for each plot the frequencies of found repetitions

from collections import Counter
import matplotlib.pyplot as plt
import os

# https://www.ncbi.nlm.nih.gov/nuccore/JC598668.1?report=fasta
S = """GCTAGGGGTTCGACGAAATTTTTTTTTATACAGTGTGTGGCCGCGAGAGGGTTAGAGGGTTTGGCACAAG
TTTGGCGGCAGCTCCAGCATCCGTTTTTGGAACTCGATTCCGATTACCGCAGGTCTCTCGGACTAGTTGA
GAGGGTAAGTGTTTTGTTTTTGTTTTATTAATTTGTTCTGTGTTTCCTTCTCCCTATACAAACGTGTTTG
TTAATTTATAGGCTTGTATGCCCCGCTTCATTTCACTCATTATGTAGGCGTTACACTGACCTTTGCACTG
CGTCGGGGATCTATGGGTTTCTGTGTGCTTACTACGATTTTAAACTTGTTATTTTACGGTTTGAATTTTT
AATACTTGTCCTCGTTGAATAGATGATAGCCTGTTAATAGATTTGAGTCAGGCATTGCGGAGTTATATGT
TACGCGGTTCCCAGCCTATAAGAATCGTGGTGTTGGCGCCAAAAAATGCGCGCCAACAAAAAGGCGCCAA
AAAATTGCGCGCCATTGTTTGGCGCCTTTTTCTGCGCCGTTTTCAAAATCGCGCCATACCAATTTCAAAG
TTCCCGCCATTTGGCCAACACGCTATTATCCCTGCATGATCTTCTTTAATTGGACGACATTCCTCGATTC
CCGATCCACATATCCAGTGACAGGAGTTCGGAATAAACGTTGTGATACGCGATCGAGTTCGCGCCATACC
AATTTCAAAGTTCCCGCCATTTGGCCAACACGCTATTATCCCTGCATGATCTTCTTTAATTGGACGACAT
TCCTCGATTCCCGATCCACATATCCAGTGACAGGAGTTCGGAATAAACGTTGTGATACGCGATCGAGTTT
TCGTGGCATATTCCTACGGAAACCTATTGTTCTGTGGTTGGTTTCGATCTATCGTTCTCGTACTGCGTGA
CCTCTACGGAACAATAGTTTTCCAGGAGATTTCCCGGTTTCGACTGCCGAAGCATGGAAACGTCCTGGGA
AAATCTGTTGTTCCGTAGTGTTCTCGTGACACTAACTCGAGATCCCTGCGAAATGACAGTTTTCTCTGGG
AATTACATCGTCCTGATTGTCGCGACATGGAATGGAAGCCTCATAGGAAGAACTCGATGTGATGATGCTC
TCTAGCCAAGAGAGCCGCGAACGCTCCAAGGAGAACTGTTATCTCGGGGAGATCTCGATCTCTCCTACCA
GCAACTCGAGATCTCTACGAGATTACAGTTTTTGGGGGAAATGTGTCCTCAGAACTGCTTAATCGTAGAA
GCTTCCTAGTGGATGGCGTTGTCTCGTAGAGGTCCAGATCTCTCCTGTTGGCAACTCGAAATCTCTACGA
GATAACAGTTTGTCTAGGAAACTTTCCTCCCAACTAAAGAGCGATGACTTAGGAAGTAAACGTGCCCTCA
TCACCGCCCTTACACACTGCTAGTCATTCATGTACATTGCGATTGTGCCTTGGTGCGGGGCGGTTCCTAG
GCACCATTTATCTTGTATTCCTGTACATCCCCTCCTTAATACTTTAATTGGAGCCACATCGTTTTGCCAT
GTTGATCTGTTTCCTTTCCCTGTAACAATGGTTACATTCTAAACGAACTAACTTGCTTTGAAGAAATCTA
CGAATTGATAGACTGTATATGCGCCCCTTTTCTTTATTAAAAACGCTTAACAAGGTAACAATATACGTAA
CATTCGAACAAGTACCAAGATAACTATATAACGGCTTATATAACTATC"""

def getRepeats(seq, minLen=6, maxLen=10):
    seq = seq.replace("\n", "")
    repeated = {}
    for length in range(minLen, maxLen+1):
        subSeqs = [seq[i:i+length] for i in range(len(seq)-length+1)]
        counts = Counter(subSeqs)
        # keep only those occurring more than once
        repeated[length] = {sub: ct for sub, ct in counts.items() if ct > 1}
    return repeated

def read_fasta(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq


repeatCountsPerGenome = {}
genomeFiles = []

for i in range(1,11):
    fname = f"Influenza{i}.fasta"
    if not os.path.exists(fname):
        print(f"[Warning] File {fname} not found; skipping.")
        continue
    genomeFiles.append(fname)
    print(f"Processing {fname} ...")
    seq = read_fasta(fname)
    repeats = getRepeats(seq,6,10)
    repeatCountsPerGenome[fname] = repeats

for fname, repeats in repeatCountsPerGenome.items():
    allItems = []
    for length, subs in repeats.items():
        for subseq, ct in subs.items():
            allItems.append((subseq + f"({length}bp)", ct))
    allItemsSorted = sorted(allItems, key=lambda x: x[1], reverse=True)
    labels, values = zip(*allItemsSorted) if allItemsSorted else ([],[])
    plt.figure(figsize=(14,6))
    plt.bar(labels, values, color='skyblue', edgecolor='black')
    plt.title(f"Repeated subsequences (6-10 bp) in {fname}")
    plt.xlabel("Subsequence (with length)")
    plt.ylabel("Occurrences")
    plt.xticks(rotation=90, fontsize=6)
    plt.tight_layout()
    plt.show()


