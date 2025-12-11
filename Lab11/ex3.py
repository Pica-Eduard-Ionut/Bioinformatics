# 3. Formulate 3 scoring equations that are able to show the level of similarity between the 2 sequences aligned in the previous assignment. Implement each of the scoring equations in your current implementation

def load_fasta(filename):
    seq = ""
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq += line.strip()
    return seq

influenza_seq = load_fasta("Influenza.fasta")
covid_seq = load_fasta("Covid19.fasta")

def smith_waterman(s1, s2, match=2, mismatch=-1, gap=-2):
    m, n = len(s1), len(s2)
    H = [[0]*(n+1) for _ in range(m+1)]
    max_score = 0
    max_pos = (0,0)

    for i in range(1, m+1):
        for j in range(1, n+1):
            score_diag = H[i-1][j-1] + (match if s1[i-1]==s2[j-1] else mismatch)
            score_up = H[i-1][j] + gap
            score_left = H[i][j-1] + gap
            H[i][j] = max(0, score_diag, score_up, score_left)
            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = (i,j)

    # Backtracking
    i,j = max_pos
    align1, align2 = "", ""
    while i>0 and j>0 and H[i][j] != 0:
        current = H[i][j]
        if current == H[i-1][j-1] + (match if s1[i-1]==s2[j-1] else mismatch):
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i -= 1
            j -= 1
        elif current == H[i-1][j] + gap:
            align1 = s1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = s2[j-1] + align2
            j -= 1
    return align1, align2, max_score


def make_segments(seq, window=300, step=150):
    segments = []
    for start in range(0, len(seq)-window, step):
        end = start + window
        segments.append((start,end,seq[start:end]))
    return segments

influenza_segments = make_segments(influenza_seq, window=300, step=150)

def compute_similarity_scores(align1, align2):
    matches = 0
    mismatches = 0
    gaps = 0
    length = len(align1)

    for a, b in zip(align1, align2):
        if a == "-" or b == "-":
            gaps += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1

    # 1. Percent Identity
    pid = matches / length * 100

    # 2. Gap-Adjusted Score
    gas = (matches - gaps) / length * 100

    # 3. Weighted Similarity Score
    wss = (2*matches - mismatches - gaps) / length * 100

    return pid, gas, wss

matches = []

for seg_start, seg_end, subseq in influenza_segments:
    align1, align2, score = smith_waterman(subseq, covid_seq)
    if score > 50:
        clean_align1 = align1.replace("-", "")
        clean_align2 = align2.replace("-", "")

        i_start = influenza_seq.find(clean_align1, seg_start)
        i_end   = i_start + len(clean_align1)
        j_start = covid_seq.find(clean_align2)
        j_end   = j_start + len(clean_align2)

        if i_start != -1 and j_start != -1:
            pid, gas, wss = compute_similarity_scores(align1, align2)
            matches.append({
                "influenza_start": i_start,
                "influenza_end": i_end,
                "covid_start": j_start,
                "covid_end": j_end,
                "match_seq_inf": clean_align1,
                "match_seq_covid": clean_align2,
                "PID": pid,
                "GAS": gas,
                "WSS": wss
            })

import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(18,8))

y_covid = 0
y_base_influenza = 1
spacing = 0.3

# Draw COVID genome line
ax.hlines(y_covid, xmin=0, xmax=len(covid_seq), color="blue", linewidth=3)
ax.text(len(covid_seq)/2, y_covid-0.05, "COVID genome", ha="center", fontsize=12)

# Draw each influenza match as a separate line with scores
for idx, m in enumerate(matches):
    y = y_base_influenza + idx*spacing
    x1, x2 = m["influenza_start"], m["influenza_end"]

    # influenza segment line
    ax.hlines(y, xmin=x1, xmax=x2, color="red", linewidth=6)

    # full match sequence text
    seq_inf = m["match_seq_inf"]
    ax.text((x1+x2)/2, y+0.05, seq_inf, fontsize=8, rotation=45, ha="center")

    # similarity scores above the segment
    score_text = f"PID:{m['PID']:.1f} GAS:{m['GAS']:.1f} WSS:{m['WSS']:.1f}"
    ax.text((x1+x2)/2, y+0.12, score_text, fontsize=7, rotation=45, ha="center")

    # connection to COVID
    cx = (m["covid_start"] + m["covid_end"])/2
    ax.plot([(x1+x2)/2, cx], [y, y_covid], color="gray", linestyle="--")

    # COVID matching sequence
    seq_cov = m["match_seq_covid"]
    ax.text(cx, y_covid+0.05, seq_cov, fontsize=8, rotation=45, ha="center")

ax.set_ylim(-0.5, y_base_influenza + len(matches)*spacing + 0.5)
ax.set_xlim(0, max(len(influenza_seq), len(covid_seq)))
ax.set_yticks([y_covid])
ax.set_yticklabels(["COVID"])
ax.set_xlabel("Genome coordinate")
plt.title("Influenza → COVID Local Alignment Matches (Full Sequences + Similarity Scores)")
plt.tight_layout()
plt.show()
