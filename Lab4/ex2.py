# download from NCBI the FASTA files containing the COVID-19 genome and the Influenza gene. Use the AI to compare the codon frequencies between the two.
# a) make a chart that shows top 10 most frequent codons for COVID-19
# b) make a chart that shows top 10 most frequent codons for Influenza
# c) compare the two results and show the most frequent codons between the two
# d) show in the output of the console top 3 amino acids for each genome

import sys
from collections import Counter
import matplotlib.pyplot as plt


lookupTable = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",

    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",

    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",

    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}


def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"): continue
            seq.append(line.strip().upper())
    return "".join(seq)

def dna_to_rna(dna):
    return dna.replace("T","U")

def codon_freqs(rna_seq, frame=0):
    freqs = Counter()
    for i in range(frame, len(rna_seq)-2, 3):
        codon = rna_seq[i:i+3]
        if len(codon)==3:
            freqs[codon] += 1
    return freqs

def aa_freqs_from_codons(codon_counts):
    aa_counts = Counter()
    for codon, count in codon_counts.items():
        aa = lookupTable.get(codon)
        if aa and aa!="Stop":
            aa_counts[aa] += count
    return aa_counts

def plot_top_codons(freqs, title, top_n=10):
    top = freqs.most_common(top_n)
    codons, counts = zip(*top)
    plt.figure(figsize=(10,6))
    plt.bar(codons, counts)
    plt.title(title)
    plt.xlabel("Codon (RNA)")
    plt.ylabel("Frequency")
    plt.show()

def main():
    # provide your paths
    covid_path = "Covid19.fasta"
    flu_path = "Influenza.fasta"

    print("Loading COVID-19 genome...")
    covid_dna = read_fasta(covid_path)
    covid_rna = dna_to_rna(covid_dna)
    covid_codons = codon_freqs(covid_rna, frame=0)
    covid_aa = aa_freqs_from_codons(covid_codons)

    print("Loading Influenza genome...")
    flu_dna = read_fasta(flu_path)
    flu_rna = dna_to_rna(flu_dna)
    flu_codons = codon_freqs(flu_rna, frame=0)
    flu_aa = aa_freqs_from_codons(flu_codons)

    # a) Top10 COVID-19 codons
    print("\nTop 10 codons – SARS-CoV-2:")
    for codon, ct in covid_codons.most_common(10):
        print(f"{codon}: {ct}")
    plot_top_codons(covid_codons, "Top 10 Codons – SARS-CoV-2")

    # b) Top10 Influenza codons
    print("\nTop 10 codons – Influenza A:")
    for codon, ct in flu_codons.most_common(10):
        print(f"{codon}: {ct}")
    plot_top_codons(flu_codons, "Top 10 Codons – Influenza A")

    # c) Compare: find codons common in each top list
    covid_top10 = {c for c,_ in covid_codons.most_common(10)}
    flu_top10 = {c for c,_ in flu_codons.most_common(10)}
    common = covid_top10 & flu_top10
    print("\nCodons common in top10 of both viruses:")
    for c in sorted(common):
        print(c)

    # d) Top3 amino acids each
    print("\nTop 3 amino acids – SARS-CoV-2:")
    for aa, ct in covid_aa.most_common(3):
        print(f"{aa}: {ct}")
    print("\nTop 3 amino acids – Influenza A:")
    for aa, ct in flu_aa.most_common(3):
        print(f"{aa}: {ct}")

if __name__=="__main__":
    main()
