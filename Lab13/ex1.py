import numpy as np

def predictSteps(transitionMatrix, initialVector, numSteps=5):
    currentState = initialVector
    predictions = [currentState]

    for _ in range(numSteps):
        currentState = transitionMatrix @ currentState
        predictions.append(currentState)

    return predictions

transitionMatrix = np.array([
    [1, 0, 0, 0],
    [1, 1, 0, 0],
    [0, 1, 1, 0],
    [0, 0, 1, 1]
])

initialVector = np.array([1, 0, 0, 0])

results = predictSteps(transitionMatrix, initialVector)

for step, state in enumerate(results):
    print(step, state)
