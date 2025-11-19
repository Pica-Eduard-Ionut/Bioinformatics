# 1. Make an artificial DNA sequence of 200-400 bases in length, in which to simulate 3-4 transposable elements.
# 2. Implement a software application to detect the positions of these transposable elements (start, end) within the created DNA sequence.


S = (
    "GCTAGGATCGTACCGATGCTAGGCTAACGATGCGTTAACCTAGGCTAACGTAGGCTAGCT"
    "ATCGTAGGCTA"
    "GGTACCGATGCTAGGCTAATCGATCGGATCCGATGCTAAGGCTAACGTAGCTAGGCTAAC"
    "GGGATTAACCTT"
    "TAGGCTAACGATGCTAGGTACCGATGCTAGGCTAAGGCTAACGATGCTAGGCTAACGATG"
    "TTAACCGGAT"
    "CGAATCGATGCTAGGCTAACGATGCTAGGTACCGATGCTAGGCTAACGATGCTAGGCTAA"
    "CCTAGGTCGA"
    "GCTAGGCTAACGTAGGCTAGCTAACGATGCTAGGCTAATCG"
)

transposableElements = {
    "TE1": "ATCGTAGGCTA",
    "TE2": "GGGATTAACCTT",
    "TE3": "TTAACCGGAT",
    "TE4": "CCTAGGTCGA",
    "TE5": "TAGGTC" # overlapping transposeable element
}


def find_te_positions(S, transposableElements):
    results = []

    for name, motif in transposableElements.items():
        mLen = len(motif)
        for i in range(len(S) - mLen + 1):
            if S[i:i+mLen] == motif:
                results.append({
                    "TE": name,
                    "start": i + 1,
                    "end": i + mLen,
                    "motif": motif
                })
    return results


hits = find_te_positions(S, transposableElements)


print("Detected Transposable Elements:\n")
for hit in hits:
    print(f"{hit['TE']}  |  Start: {hit['start']}   End: {hit['end']}  |  Motif: {hit['motif']}")
