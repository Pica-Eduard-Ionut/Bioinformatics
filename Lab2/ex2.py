# Find in sequence S only the dinucleotides and trinucleotides that exists, without the use of the bruteforce engine
# in order to achieve the results one must verify these combinations starting from the beginning of the sequence
# example: S = "abaa" => dinucleotides(ab,ba,aa)
S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

dinucleotides = []
trinucleotides = []

for i in range(0,len(S)-1):
    di = S[i:i+2]
    if di not in dinucleotides:
        dinucleotides.append(di)

for i in range(0,len(S)-2):
    tri = S[i:i+3]
    if tri not in trinucleotides:
        trinucleotides.append(tri)

print(dinucleotides)

print(trinucleotides)
print(len(trinucleotides))
