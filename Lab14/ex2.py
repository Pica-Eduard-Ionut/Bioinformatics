# Take two poetries from the internet one from Mihai Eminescu and the other from Nichina Stanescu, make sure that the two poetries are of almost equal size.
# Use these 2 texts to make 2 transition matrices which capture the transition probabilities between the words inside these texts.
# Use a new text which combines other poetries from the two poets.
# Generate a log-likelihood matrix inbetween the two models and use this matrix in order to scan the third text which is a combination of other poetries of the two authors.
# Your scanner must be able to compute the scores of a sliding window and show the user on a chart which pieces of this text are from Mihai Eminescu and which pieces are from Nichita Stanescu.
# !!! Combine from the first two poems for the mixed one

text_eminescu = """
Vântu-o foaie vestejită
Mi-au adus mișcând fereasta -
Este moartea ce-mi trimite
Fără plic scrisoarea-aceasta.

Voi păstra-o, voi întinde-o
Între foile acele,
Ce le am din alte timpuri
De la mâna dragei mele.

Cum copacu-și uită foaia
Ce pe vânt mi-a fost trimisă,
Astfel ea uitat-au poate
Aste foi de dânsa scrise.

Vorbele iubirii moarte
Vinovate-mi stau de față,
Dovedite de minciună
Cer să sting a lor viață.

Dulcea lor zădărnicie
Nu mă-ndur s-o pun pe foc,
Deși-mi stau atât de triste
Că nu pot muri pe loc.

Voi păstra întreg amarul
Și norocul ăstor foi,
În durerea vechii pierderi
Recitindu-mă-napoi;

Numai vestea blând-a morții,
Foaia tristă le-am adaos:
Moartea vindec-orice rană,
Dând la patime repaos.
"""

text_stanescu = """
Somnul cu fierăstraie-n el
taie capetele cailor
și caii aleargă nechezind cu sânge,
ca niște mese roșii, fugite pe străzi,
de la cina cea de taină.

Și caii aleargă, în aburii roșii
clătinând umbre. În șăi, fantome.
Frunze se lipesc de gâturile lor
sau se prăbușesc de-a dreptul în ele,
cum se prăbușește umbra copacului în fântâni.

Aduceți găleți, aduceți căni mari de sticlă,
aduceți căni și pahare.
Aduceți căștile vechi rămase din război,
aduceți-i pe toți cărora le lipsește un ochi,
sau în loc de braț au un loc liber,
care poate fi umplut.

Peste tot sânge de cal decapitat
curge în voie, și
eu cel care-am văzut primul
acestea
vă vestesc că am băut din el ...
"""

# mixed by chatgpt and it resulted in uniform distribution
mixed_text = """
Vântu-o foaie vestejită
Mi-au adus mișcând fereasta -
Este moartea ce-mi trimite
Fără plic scrisoarea-aceasta.

Somnul cu fierăstraie-n el
taie capetele cailor
și caii aleargă nechezind cu sânge,
ca niște mese roșii, fugite pe străzi.

Voi păstra-o, voi întinde-o
Între foile acele,
Ce le am din alte timpuri
De la mâna dragei mele.

Frunze se lipesc de gâturile lor
sau se prăbușesc de-a dreptul în ele,
cum se prăbușește umbra copacului în fântâni.

Cum copacu-și uită foaia
Ce pe vânt mi-a fost trimisă,
Astfel ea uitat-au poate
Aste foi de dânsa scrise.

Aduceți găleți, aduceți căni mari de sticlă,
aduceți căni și pahare.
Aduceți căștile vechi rămase din război,
aduceți-i pe toți cărora le lipsește un ochi.

Vorbele iubirii moarte
Vinovate-mi stau de față,
Dovedite de minciună
Cer să sting a lor viață.

Peste tot sânge de cal decapitat
curge în voie, și
eu cel care-am văzut primul
acestea
vă vestesc că am băut din el ...
"""


import re
from collections import defaultdict
import math
import matplotlib.pyplot as plt

def tokenize(text):
    text = text.lower()
    words = re.findall(r"\b\w+\b", text)
    return words

def transition_counts(words):
    counts = defaultdict(lambda: defaultdict(int))
    for i in range(len(words) - 1):
        counts[words[i]][words[i+1]] += 1
    return counts

def transition_probabilities(counts, vocab, alpha=1):
    probs = defaultdict(dict)
    for w1 in vocab:
        total = sum(counts[w1][w2] for w2 in vocab) + alpha * len(vocab)
        for w2 in vocab:
            probs[w1][w2] = (counts[w1][w2] + alpha) / total
    return probs

words_E = tokenize(text_eminescu)
words_S = tokenize(text_stanescu)

vocab = set(words_E) | set(words_S)

counts_E = transition_counts(words_E)
counts_S = transition_counts(words_S)

P_E = transition_probabilities(counts_E, vocab)
P_S = transition_probabilities(counts_S, vocab)

LL = defaultdict(dict)

for w1 in vocab:
    for w2 in vocab:
        LL[w1][w2] = math.log2(P_E[w1][w2] / P_S[w1][w2])

def scan_text(words, LL, window_size=12, step=1):
    scores = []
    positions = []

    for i in range(0, len(words) - window_size + 1, step):
        score = 0
        transitions = 0

        for j in range(window_size - 1):
            w1 = words[i + j]
            w2 = words[i + j + 1]

            if w1 in LL and w2 in LL[w1]:
                score += LL[w1][w2]
                transitions += 1

        if transitions > 0:
            score /= transitions

        scores.append(score)
        positions.append(i)

    return positions, scores


mixed_words = tokenize(mixed_text)

positions, scores = scan_text(mixed_words, LL, window_size=15, step=1)


plt.figure(figsize=(12, 5))
plt.plot(positions, scores, label="Log-likelihood score")
plt.axhline(0, color="black", linestyle="--")

plt.title("Authorship Detection: Eminescu vs Stănescu")
plt.xlabel("Word position (sliding window)")
plt.ylabel("Log₂ likelihood score")

plt.text(0, max(scores)*0.9, "Eminescu", color="green")
plt.text(0, min(scores)*0.9, "Stănescu", color="red")

plt.show()
