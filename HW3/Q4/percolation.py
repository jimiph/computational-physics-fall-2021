import numpy as np
import matplotlib.pyplot as plt


def percolation(lengthOfLattice = 10, probability = 0.5, isShowLabeledLattice = False,
                isShowBinaryLattice = False, deltaP = 0.05,
                numOfRuns = 100, isShowPlotOfProbabilityOfPercolation = False):
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
    This function shows the made binary random valued lattice.
    '''
    def showBinaryLattice():
        fig , ax = plt.subplots(figsize=(10, 7))
        m = makeBinaryLattice()
        binaryLattice = np.zeros((lengthOfLattice, lengthOfLattice), dtype=int)
        for i in range(lengthOfLattice):
            for j in range(lengthOfLattice):
                binaryLattice[i][j] = m[i+1][j+1]
        ax.imshow(binaryLattice, interpolation = 'none')
        ax.invert_yaxis()
        ax.set_xlabel(f'The lattice is {lengthOfLattice}*{lengthOfLattice}')
        ax.set_title(f'The labeled lattice with 0 and 1 (p = {probability})')
        plt.show()
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
        
        return labeledLattice, isPercolate, labelsOfInfiniteClusters
    '''
    This function shows the made labeled lattice.
    '''
    def showLabeledLattice():
        labeledLattice, percolationStatus, c = makeLabeledLattice()
        fig , ax = plt.subplots(figsize=(10, 7))
        ax.pcolormesh(labeledLattice, cmap = plt.get_cmap('viridis'))
        ax.set_xlabel(f'The lattice is {lengthOfLattice}*{lengthOfLattice} and p={probability}')
        ax.set_title(f'Percolation Status:{percolationStatus}')
        plt.show()
    '''
    This function plot the probability distribution of making an infinite cluster
    and connection of a random chosed point to it.
    '''
    def infiniteClusterProbabilityDisturb():
        def calculatePercent(array):
            numOfTrues = 0
            for x in array:
                if x == True:
                    numOfTrues += 1
            return (numOfTrues / len(array))
        pValues = np.array([i * deltaP for i in range(1 + int(1/deltaP))])
        qConnectionValues = np.zeros(1 + int(1/deltaP))
        status = np.zeros(numOfRuns)
        lengthValues = [10, 100, 200]
        colors = ['r', 'b', 'g']
        fig , ax = plt.subplots(figsize=(10, 7))
        # ax.set_title('Probability of finding a cluster as a function of P')
        ax.set_title('Probability of finding a connection with the infinite cluster as a function of P')
        ax.set_xlabel('P: probability of turning a cell on')
        # ax.set_ylabel('Q: probability of percloation')
        ax.set_ylabel('Q: probability of connetion')
        for c in range(len(colors)):
            connectionValues = np.zeros(numOfRuns)
            qValues = np.zeros(1 + int(1/deltaP))
            for i in range(len(pValues)):
                for run in range(numOfRuns):
                    labeledLattice, status[run], labelsOfInfiniteClusters = makeLabeledLattice(probability=pValues[i],
                     lengthOfLattice=lengthValues[c])
                    x = np.random.randint(lengthValues[c])
                    y = np.random.randint(lengthValues[c])
                    if labeledLattice[x][y] in labelsOfInfiniteClusters:
                        connectionValues[run] = True
                qConnectionValues[i] = calculatePercent(connectionValues)
                qValues[i] = calculatePercent(status)
            # ax.plot(pValues, qValues, 'o-', color = f'{colors[c]}', label=f'L = {lengthValues[c]}')
            ax.plot(pValues, qConnectionValues, 'o-', color=f'{colors[c]}', label=f'{lengthValues[c]}')
        plt.legend()
        plt.show()

        pass
    if isShowLabeledLattice == True:
        showLabeledLattice()
    if isShowBinaryLattice == True:
        showBinaryLattice()
    if isShowPlotOfProbabilityOfPercolation == True:
        infiniteClusterProbabilityDisturb()

#To see the random values lattice, uncomment the below line:
# percolation(isShowBinaryLattice=True)


#To see the diagram of probability of percolation or probability of finding a
#connection to the infinite lattice, uncomment the below line:
# percolation(isShowPlotOfProbabilityOfPercolation=True)


#To see the labeled lattice using HK's algorithm, uncomment the below line:
# percolation(isShowLabeledLattice=True)
