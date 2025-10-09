# Design an application using the AI, which contains a GUI that allows the user to choose a fasta file. The content of the file should be analyzed by using a sliding window of 30 positions. The content for each sliding window
# should be used in order to extract the relative frequencies (percentages) of the symbols found in the alphabet of the sequence thus your input will be the DNA seq from the fasta file and the output should be the values of
# the relative frequencies of each symbol from the alphabet of the sequence. Translate it as lines on a chart, thus your chart in the case of DNA should have 4 lines which reflect the values found over the sequence
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def read_fasta(filepath):
    sequence = ""
    try:
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line.startswith(">"):
                    sequence += line.upper()
        return sequence
    except Exception as e:
        messagebox.showerror("Error", f"Could not read file: {e}")
        return None

def sliding_window_frequencies(seq, window_size=30):
    if len(seq) < window_size:
        messagebox.showerror("Error", "Sequence is shorter than window size.")
        return None

    positions = []
    freqs = {'A': [], 'T': [], 'C': [], 'G': []}

    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        total = len(window)
        positions.append(i + 1)

        for base in freqs:
            count = window.count(base)
            freqs[base].append((count / total) * 100)

    return positions, freqs

def plot_frequencies(positions, freqs):
    plt.figure(figsize=(12, 6))

    colors = {'A':'green', 'T':'red', 'C':'blue', 'G':'orange'}

    def smooth(values, window=10):
        return np.convolve(values, np.ones(window)/window, mode='same')

    for base in freqs:
        y_smooth = smooth(freqs[base], window=10)
        plt.plot(positions, y_smooth, label=base, color=colors[base], linewidth=2)

    plt.title("Relative Nucleotide Frequencies (Sliding Window = 30)")
    plt.xlabel("Sequence Position (Window Start)")
    plt.ylabel("Frequency (%)")
    plt.ylim(0, 100)
    plt.legend(title="Nucleotide")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()

def open_fasta():
    filepath = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
    )
    if not filepath:
        return

    seq = read_fasta(filepath)
    if seq:
        result = sliding_window_frequencies(seq, 30)
        if result:
            positions, freqs = result
            plot_frequencies(positions, freqs)

root = tk.Tk()
root.title("DNA Frequency Analyzer")

label = tk.Label(root, text="Choose a FASTA file to analyze (Sliding Window = 30):")
label.pack(pady=10)

button = tk.Button(root, text="Select FASTA File", command=open_fasta)
button.pack(pady=5)

root.mainloop()
