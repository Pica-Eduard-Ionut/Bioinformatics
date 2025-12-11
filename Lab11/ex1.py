import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"

match_score = 1
mismatch_score = -1
gap_penalty = -2


def needleman_wunsch(seq1, seq2):
    n, m = len(seq1), len(seq2)
    score = np.zeros((n+1, m+1), dtype=int)
    traceback = np.zeros((n+1, m+1), dtype=str)

    for i in range(1, n+1):
        score[i][0] = i * gap_penalty
        traceback[i][0] = "U"

    for j in range(1, m+1):
        score[0][j] = j * gap_penalty
        traceback[0][j] = "L"

    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            up = score[i-1][j] + gap_penalty
            left = score[i][j-1] + gap_penalty

            best = max(diag, up, left)
            score[i][j] = best

            traceback[i][j] = "D" if best == diag else "U" if best == up else "L"

    aligned1, aligned2 = "", ""
    path = []
    i, j = n, m
    path.append((i, j))

    while i > 0 or j > 0:
        if traceback[i][j] == "D":
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        elif traceback[i][j] == "U":
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

        path.append((i, j))

    return aligned1, aligned2, score, set(path)

aligned1, aligned2, score_matrix, path = needleman_wunsch(S1, S2)

matches = sum(a == b for a, b in zip(aligned1, aligned2))
length = len(aligned1)
similarity = matches / length * 100

print("Alignment:")
print(aligned1)
print(aligned2)
print()
print("Matches:", matches)
print("Length:", length)
print(f"Similarity: {similarity:.2f}%")
print("Traceback end at:", max(path))

purple_red_cmap = LinearSegmentedColormap.from_list(
    "purple_red",
    [
        "#000020",  # deep blue
        "#240040",  # purple
        "#550055",  # dark magenta
        "#770033",  # plum
        "#aa0022",  # red
        "#ff0033"   # bright red
    ]
)

norm_score = score_matrix.astype(float)
norm_score -= norm_score.min()
norm_score /= norm_score.max()

def create_traceback_grid(n, m, path):
    grid = np.zeros((n+1, m+1))
    for i in range(n+1):
        for j in range(m+1):
            if (i, j) in path:
                grid[i][j] = 1
    return grid


trace_grid = create_traceback_grid(len(S1), len(S2), path)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.imshow(norm_score, cmap=purple_red_cmap, interpolation="nearest")
ax1.set_title("Graphic representation of the alignment matrix")
ax1.set_xticks([])
ax1.set_yticks([])

ax2.imshow(trace_grid, cmap=LinearSegmentedColormap.from_list("yellow_red", ["#ffffaa", "#ff0000"]))
ax2.set_title("Traceback path deviation")
ax2.set_xticks([])
ax2.set_yticks([])

plt.tight_layout()
plt.show()
