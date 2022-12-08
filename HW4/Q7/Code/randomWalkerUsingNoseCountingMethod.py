import numpy as np
import matplotlib.pyplot as plt
from math import comb


def randomWalkerByNoseCounting(limit=0.9999, numOfSites=20, numOfTimeSteps=1000,
                                     pRightValue=0.5, initialPosition=1):
    #This function makes the array of probabilities then return it.
    def makeProbabilities(initialPosition=initialPosition):
        probabilities = np.zeros((numOfTimeSteps, numOfSites))
        x0 = initialPosition
        probabilities[0][x0] = 1
        p = pRightValue
        q = 1 - p
        for y in range(1, numOfTimeSteps):
            for x in range(numOfSites):
                if x == 0 or x == 1:
                    probabilities[y][x] = q * probabilities[y-1][x+1]
                elif x == numOfSites-1 or x == numOfSites-2:
                    probabilities[y][x] = p * probabilities[y-1][x-1]
                else:
                    probabilities[y][x] = p * probabilities[y-1][x-1] + q * probabilities[y-1][x+1]
        return probabilities
    #This function calculate the life expectancy for a given probabilities array.
    def calculateLifeExpectancy(initialPosition):
        probabilities = makeProbabilities(initialPosition)
        lifeExpectency = 0
        probabilityOfDeath = np.zeros(numOfTimeSteps)
        probabilityOfDeath[0] = probabilities[0][0] + probabilities[0][numOfSites-1]
        i = 0
        while np.sum(probabilityOfDeath) < limit and i < numOfTimeSteps:
            probabilityOfDeath[i] = probabilities[i][0] + probabilities[i][numOfSites-1]
            lifeExpectency += probabilityOfDeath[i] * i
            i += 1
        return lifeExpectency
    #This function calculate the life expectancy for each initial position.
    def plotLifeExpectancyPerInitialPosition():
        lifeExpectencies = np.zeros(numOfSites)
        for i in range(numOfSites):
            lifeExpectencies[i] = calculateLifeExpectancy(i)
        sites = [x for x in range(numOfSites)]
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.plot(sites, lifeExpectencies, '-o')
        ax.set_title(f'Bounds are at 0 and {numOfSites} and probability of going to right is {pRightValue}')
        ax.set_ylabel('life expectency')
        ax.set_xlabel('initial position')
        plt.show()

    plotLifeExpectancyPerInitialPosition()


randomWalkerByNoseCounting(pRightValue=0.7)







