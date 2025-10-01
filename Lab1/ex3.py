# use the AI to adapt your current algorithm in order to make an app that takes a FASTA file and read the seq content from it and display the relational percentages for the symbols present in the
# alphabet of seq. NOTE: FASTA represents a file format that contains DNA, ARN or proteins seq thus, it contains the file information for your input

import tkinter as tk
from tkinter import filedialog


def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                sequence += line.upper()
    return sequence


def calculate_percentages(sequence):
    alphabet = sorted(set(sequence))
    freq = {symbol: 0 for symbol in alphabet}

    for char in sequence:
        freq[char] += 1

    percentages = {}
    for symbol in alphabet:
        percentages[symbol] = round(freq[symbol] / len(sequence) * 100, 2)

    return alphabet, freq, percentages


def open_file():
    fasta_file = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
    )

    if fasta_file:
        sequence = read_fasta(fasta_file)
        alphabet, freq, percentages = calculate_percentages(sequence)

        # print("\nFile chosen:", fasta_file)
        # print("\nAlphabet:", "".join(alphabet))
        # print("\nFrequency:", freq)
        # print("\nRelative frequencies:")
        # for symbol in alphabet:
        #     print(f"{symbol}: {percentages[symbol]} %")

        output_text.delete(1.0, tk.END)
        output_text.insert(tk.END, f"File chosen: {fasta_file}\n\n")
        output_text.insert(tk.END, f"Alphabet: {''.join(alphabet)}\n\n")
        output_text.insert(tk.END, f"Frequency: {freq}\n\n")
        output_text.insert(tk.END, "Relative frequencies:\n")
        for symbol in alphabet:
            output_text.insert(tk.END, f"{symbol}: {percentages[symbol]} %\n")


# UI
root = tk.Tk()
root.title("FASTA Analyzer")

btn = tk.Button(root, text="Choose FASTA File", command=open_file, width=30, height=2)
btn.pack(pady=10)

output_text = tk.Text(root, width=60, height=20, wrap="word")
output_text.pack(padx=10, pady=10)

root.mainloop()
