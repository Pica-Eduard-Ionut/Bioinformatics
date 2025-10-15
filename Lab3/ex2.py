# Design an application that uses the sliding window method in order to read the melting temperature(tm) over the sequence S. Use a sliding window of 8 positions and choose a fasta file as input.
# Plot it in a chart as last time
import math
import matplotlib.pyplot as plt
from tkinter import filedialog, Tk

def calc_freq(seq):
    entries = "AGTC"
    freq = {i: 0 for i in entries}
    for i in seq:
        if i in freq:
            freq[i] += 1
    return freq

def calc_tm(freq):
    # Simple Wallace rule
    return 4 * (freq["G"] + freq["C"]) + 2 * (freq["A"] + freq["T"])

# Na+ = 0.001 M
def calc_tm2(seq, freq):
    Na = 0.001
    GC_percentage = (freq["G"] + freq["C"]) / len(seq) * 100
    tm = 81.5 + 16.6 * math.log10(Na) + 0.41 * GC_percentage - 600 / len(seq)
    return tm

def sliding_window_tm(seq, window_size=8):
    tm_values = []
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        freq = calc_freq(window)
        tm = calc_tm(freq)
        tm_values.append((i + 1, tm))
    return tm_values

def sliding_window_tm2(seq, window_size=8):
    tm_values = []
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        freq = calc_freq(window)
        tm = calc_tm2(window, freq)
        tm_values.append((i + 1, tm))
    return tm_values

def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                sequence += line.upper()
    return sequence

if __name__ == "__main__":
    root = Tk()
    root.withdraw()

    fasta_file = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
    )

    if fasta_file:
        sequence = read_fasta(fasta_file)
        tm_values1 = sliding_window_tm(sequence)
        tm_values2 = sliding_window_tm2(sequence)

        positions = [pos for pos, _ in tm_values1]
        tm1_list = [tm for _, tm in tm_values1]
        tm2_list = [tm for _, tm in tm_values2]

        plt.figure(figsize=(10, 5))
        plt.plot(positions, tm1_list, label="First formula", color="blue", linewidth=2)
        plt.plot(positions, tm2_list, label="Second formula", color="red", linewidth=2)
        plt.title("Melting Temperature (Tm) Over Sequence")
        plt.xlabel("Window Start Position")
        plt.ylabel("Melting Temperature (°C)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
