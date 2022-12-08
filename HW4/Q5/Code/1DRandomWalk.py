import numpy as np
import matplotlib.pyplot as plt
from numpy.lib import average


# I use this function to check the truness of the equations (4) and (5) of the reference book.
def randomWalkTest(numOfTimeSteps=1000, numOfWalkers=10000,
                   pRightValues=[0.1, 0.4, 0.5, 0.9], isAnalyseAndPlot=False):

    #This function calculate the position of random walker per time and save the 
    #data as a .npy file.
    def makePositionsOfWalkerPerTime(pRightValue=0.1, fileName='test.npy'):
        positionsPerTime = np.zeros((numOfTimeSteps, numOfWalkers), dtype=int)
        positionsPerTime[0][:] = 0
        for walker in range(numOfWalkers):
            for t in range(1, numOfTimeSteps):
                positionsPerTime[t][walker] = positionsPerTime[t-1][walker] + np.random.choice(
                    [1, -1], p=[pRightValue, 1-pRightValue])

        np.save(fileName, positionsPerTime)


    #This function makes the data (position of random walker per time) for several 
    #values of probability.
    def makeEnsemblesForPValues():
        filesNames = [f'randomWalkerPRight={p}Q5.npy' for p in pRightValues]
        for x in range(len(pRightValues)):
            makePositionsOfWalkerPerTime(
                fileName=filesNames[x], pRightValue=pRightValues[x])
        return filesNames


    # This function uses the function makeEnsembles() to calculate <x(t)> and Var[(x(t))]
    def analyser():
        # If the data sheets aren't there, use this code first to make them.
        # filesNames = makeEnsemblesForPValues()
        # If you have data sheets use this code:
        filesNames = [f'randomWalkerPRight={p}Q5.npy' for p in pRightValues]
        averagePositions = []
        variancePositions = []
        time = [x for x in range(numOfTimeSteps)]
        for x in range(len(pRightValues)):
            file = filesNames[x]
            positionsPerTime = np.load(file)
            averagePositions.append(np.average(positionsPerTime, axis=1))
            variancePositions.append(np.var(positionsPerTime, axis=1))
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 7))
        ax1.set_title('average of position per time')
        ax1.set_xlabel(f'time: 0 - {numOfTimeSteps}')
        ax1.set_ylabel('<x(t)>')
        ax2.set_title('variance of position per time')
        ax2.set_xlabel(f'time: 0 - {numOfTimeSteps}')
        ax2.set_ylabel('Var[x(t)]')
        for x in range(len(pRightValues)):
            ax1.plot(time, averagePositions[:]
                     [x], label=f'p={pRightValues[x]}')
            ax2.plot(time, variancePositions[:]
                     [x], label=f'p={pRightValues[x]}')
        ax1.legend()
        ax2.legend()
        ax1.grid()
        ax2.grid()
        fig.tight_layout(pad=1.0)
        plt.show()
    if isAnalyseAndPlot == True:
        analyser()


# To show the plot of average postion per time and variance of position per time:
randomWalkTest(isAnalyseAndPlot=True)
