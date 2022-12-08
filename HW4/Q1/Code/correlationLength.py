import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sp



def percolation(lengthOfLattice = 10, probability = 0.7, numOfRuns = 100):
    '''
    This function make a random vlaued by 0,1 lattice.
    '''
    def makeBinaryLattice(probability=probability, lengthOfLattice=lengthOfLattice):
        randomVluesMatrix = np.random.rand(lengthOfLattice, lengthOfLattice)
        binaryMatrix = np.zeros((lengthOfLattice+1, lengthOfLattice+1), dtype=int)
        for i in range(1, lengthOfLattice+1):
            for j in range(1, lengthOfLattice+1):
                if randomVluesMatrix[i-1][j-1] <= probability:
                    binaryMatrix[i][j] = 1
                else:
                    binaryMatrix[i][j] = 0
        return binaryMatrix

    '''
    This function make a labeled lattice using HK's algorithm and has 
    some other features.
    '''
    def makeLabeledLattice(probability=probability, lengthOfLattice=lengthOfLattice):
 
        def mergeClusterLabels(x, y):
            labelValues[findClusterLabel(x)] = findClusterLabel(y)

        def findClusterLabel(x):
            y = x
            while (labelValues[y] != y):
                y = labelValues[y]
            while (labelValues[x] != x):
                z = labelValues[x]
                labelValues[x] = y
                x = z
            return y

        blockedBinaryMatrix = makeBinaryLattice(probability, lengthOfLattice)
        largestLabel = 0
        labelGrid = np.zeros((lengthOfLattice + 1, lengthOfLattice + 1), dtype=int)
        labelValues = [x for x in range(lengthOfLattice * lengthOfLattice)]
        for x in range(1, lengthOfLattice + 1):
            for y in range(1, lengthOfLattice + 1):
                if blockedBinaryMatrix[x][y]:
                    leftLabel = labelGrid[x][y-1]
                    upLabel = labelGrid[x-1][y]
                    if (leftLabel == 0 and upLabel == 0):
                        largestLabel = largestLabel + 1
                        labelGrid[x][y] = largestLabel
                    elif (leftLabel != 0 and upLabel == 0):
                        labelGrid[x][y] = findClusterLabel(leftLabel)
                    elif (leftLabel == 0 and upLabel != 0):
                        labelGrid[x][y] = findClusterLabel(upLabel)
                    else:
                        mergeClusterLabels(leftLabel, upLabel)
                        labelGrid[x][y] = findClusterLabel(leftLabel)
        
        labeledLattice = np.zeros((lengthOfLattice, lengthOfLattice), dtype=int)
        for i in range(lengthOfLattice):
            for j in range(lengthOfLattice):
                labeledLattice[i][j] = findClusterLabel(labelGrid[i+1][j+1])
        isPercolate = False
        labelsOfInfiniteClusters = []
        for x in labeledLattice[0]:
            if x != 0:
                for y in labeledLattice[lengthOfLattice-1]:
                    if y != 0 and x == y:
                        isPercolate = True
                        if x not in labelsOfInfiniteClusters:
                            labelsOfInfiniteClusters.append(x)
        
        labelsOfNoneInfiniteClusters = set({})
        for x in range(lengthOfLattice):
            for y in range(lengthOfLattice):
                if labeledLattice[x][y] != 0 and labeledLattice[x][y] not in labelsOfInfiniteClusters:
                    labelsOfNoneInfiniteClusters.add(labeledLattice[x][y])
        if len(labelsOfNoneInfiniteClusters) != 0:
            areaOfEachNoneInfiniteCluster = np.zeros(len(labelsOfNoneInfiniteClusters))
            i = 0
            for label in labelsOfNoneInfiniteClusters:
                for x in range(lengthOfLattice):
                    for y in range(lengthOfLattice):
                        if labeledLattice[x][y] == label:
                            areaOfEachNoneInfiniteCluster[i] += 1
                i += 1
            correlationLength = np.sqrt(np.max(areaOfEachNoneInfiniteCluster))
            return correlationLength
        else:
            return 0


    '''
    This function plot the correlation length per probability of turning a cell on.
    '''
    def correlationLenghtPerProbability():
        pValues = [0,0.1,0.2,0.3,0.4,0.5,0.57,0.575,0.58,0.585,0.59,
        0.595,0.60,0.605,0.61,0.615,0.62,0.63,0.635,0.7,0.8,0.9,1]
        lengthValues = [10, 15, 23, 34, 51, 76, 114, 171]
        probabilityOfMaxCorrelationLength = np.zeros(len(lengthValues))
        fig, ax1 = plt.subplots(figsize=(15, 7))
        ax1.set_title('correlation length per probability')
        ax1.set_xlabel('probability of turning on a random cell')
        ax1.set_ylabel('correlation length')
        for l in range(len(lengthValues)):
            averageCorrelationLengths = np.zeros(len(pValues))
            for p in range(len(pValues)):
                averageCorrelationLength = 0
                for r in range(numOfRuns):
                    averageCorrelationLength += makeLabeledLattice(probability=pValues[p],
                     lengthOfLattice=lengthValues[l])
                averageCorrelationLength = averageCorrelationLength / numOfRuns 
                averageCorrelationLengths[p] = averageCorrelationLength
            probabilityOfMaxCorrelationLength[l] = pValues[np.argmax(averageCorrelationLengths)]
            np.save(f'averageCorrelationLengthsL={lengthValues[l]}.npy', averageCorrelationLengths)
            ax1.plot(pValues, averageCorrelationLengths, 'o-', label=f'L = {lengthValues[l]}')
        np.save('probabilityOfMaxCorrelationLength.npy', probabilityOfMaxCorrelationLength)
        ax1.legend()
        plt.show()

    def criticalPower():
        probabilityOfMaxCorrelationLength = np.load('probabilityOfMaxCorrelationLength.npy')
        lengthValues = [10, 15, 23, 34, 51, 76, 114, 171]
        print(probabilityOfMaxCorrelationLength)
        slope, intercept, r, pp, se = sp.linregress(np.log(lengthValues), 
                                                    np.log(np.abs(probabilityOfMaxCorrelationLength-0.593)))
        fig, ax2 = plt.subplots(figsize=(15, 7))
        ax2.set_title('critical power')
        ax2.set_xlabel('log(L)')
        ax2.set_ylabel('log(|Pc(L) - Pc_inf|)')
        x = np.log(lengthValues)
        ax2.scatter(np.log(lengthValues), np.log(np.abs(probabilityOfMaxCorrelationLength-0.593)))
        ax2.plot(x, slope * x + intercept, label='fit with equation: {:.3} * ln(L) + {:0.3}'.format(slope, intercept))
        ax2.legend()
        plt.show()
        
    correlationLenghtPerProbability()
    criticalPower()



percolation()


