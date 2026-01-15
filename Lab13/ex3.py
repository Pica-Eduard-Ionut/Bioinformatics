# 4. Use the transition matrix from the json output in order to syntethize new sequences of text based on the transition matrix.

import json
import random

def weightedChoice(probabilities):
    symbols = list(probabilities.keys())
    weights = list(probabilities.values())
    return random.choices(symbols, weights=weights, k=1)[0]

def generateSequence(transitionMatrix, startSymbol, length):
    sequence = [startSymbol]
    currentSymbol = startSymbol

    for _ in range(length - 1):
        if currentSymbol not in transitionMatrix:
            currentSymbol = random.choice(list(transitionMatrix.keys()))
        else:
            currentSymbol = weightedChoice(transitionMatrix[currentSymbol])
        sequence.append(currentSymbol)

    return sequence

with open("word_transition_matrix.json", "r") as f:
    data = json.load(f)

decoder = data["encoding"]
transitionMatrix = data["transitionMatrix"]

startSymbol = random.choice(list(transitionMatrix.keys()))
generatedSymbols = generateSequence(transitionMatrix, startSymbol, 30)

generatedText = " ".join(decoder[symbol] for symbol in generatedSymbols)

print(generatedText)
