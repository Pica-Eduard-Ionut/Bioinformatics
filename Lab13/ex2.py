# 1. Use a random DNA sequence of about 50 letters. Use the sequence to compute the transition probabilities between letter. Your output should be the transition matrix stored as a json file.
# 2. Use a random text (english) of about 300 letters ( that implies spaces punctuation and everything that composes the english language ) and compute the transition probabilities between words.
# Store the transition matrix as a json file. For ease of implementation you could represent each new word by using one symbol of your choosing (ascii format can be used entirely).

import random
import json
import string
from collections import defaultdict

def computeTransitionMatrix(sequence):
    transitions = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)

    for currentSymbol, nextSymbol in zip(sequence[:-1], sequence[1:]):
        transitions[currentSymbol][nextSymbol] += 1
        totals[currentSymbol] += 1

    probabilities = {}
    for symbol in transitions:
        probabilities[symbol] = {
            nextSymbol: count / totals[symbol]
            for nextSymbol, count in transitions[symbol].items()
        }

    return probabilities

dnaAlphabet = ["A", "C", "G", "T"]
dnaSequence = [random.choice(dnaAlphabet) for _ in range(50)]

dnaTransitionMatrix = computeTransitionMatrix(dnaSequence)

with open("dna_transition_matrix.json", "w") as f:
    json.dump(dnaTransitionMatrix, f, indent=2)

randomText = "In many circumstances, an operation on a data structure can be made more efficient by including additional information. This technique is known as using augmented data structures. For example, threaded trees use additional pointers so that it is possible to find the pre-order successor to a given node."

words = randomText.split()
encodedWords = {}
decodedWords = {}
encodedSequence = []

currentCode = 33

for word in words:
    if word not in encodedWords:
        encodedWords[word] = chr(currentCode)
        decodedWords[chr(currentCode)] = word
        currentCode += 1
    encodedSequence.append(encodedWords[word])

wordTransitionMatrix = computeTransitionMatrix(encodedSequence)

with open("word_transition_matrix.json", "w") as f:
    json.dump({
        "encoding": decodedWords,
        "transitionMatrix": wordTransitionMatrix
    }, f, indent=2)
