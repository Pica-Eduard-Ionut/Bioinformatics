# The melting temperature (tm) is the temp at which one part of a particular DNA will disociate and become a single strand of DNA primerland and seq are of critical importance in designing the parameters of a successful
# amplification. The tm of a nucleic acid duplex increases both with its length and with increasing GC content
#    tm = 4(G+C) + 2(A+T)
#    tm = 81.5+16.6(log10([Na^+]))+.41*(%GC)-600/length
# Implement an app that computes the tm of a DNA seq, by using one of the formula or both of them. Input = a string of DNA ; Output = temp in Celsius
# Obs: sequence should be between 6-12
import math
S = "TTTTAACCGCTC"

def calc_freq(seq):
    entries = "AGTC"
    freq = {}
    for i in entries:
        freq[i] = 0

    for i in seq:
        freq[i] = freq[i] + 1
        #print(i)
    return freq

def calc_tm(freq):
    tm = 4*(freq["G"]+freq["C"]) + 2*(freq["A"]+freq["T"])
    return tm

# Na^+ as 0.001
def calc_tm2(seq,freq):
    Na = 0.001
    GC_percentage = (freq["G"] + freq["C"]) / len(seq) * 100
    tm = 81.5 + 16.6 * math.log10(Na) + 0.41 * GC_percentage - 600 / len(seq)
    return tm


freq = calc_freq(S)

print(str(calc_tm(freq)) + " degrees Celcius")
print(str(calc_tm2(S,freq)) + " degrees Celcius")
