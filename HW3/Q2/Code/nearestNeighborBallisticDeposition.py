import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats as scs


def ballisticSideDeposition(numOfParticles=int(1.23*10**6), numOfTimeSteps=40,
                        numOfSites=50, numOfEnsembles=1, numOfSubLayers=4,
                        isShowLayer = False, isShowPlotsOfSTD = False):
    if isShowPlotsOfSTD == True:
        timeSteps = [int(1000 * 1.2 ** x) for x in range(numOfTimeSteps)]
        numOfTimeStepsInEachSubLayer = int(numOfTimeSteps / numOfSubLayers) 
    elif isShowLayer == True:
        timeSteps = [int(numOfParticles/numOfTimeSteps) for x in range(numOfTimeSteps)]
        numOfTimeStepsInEachSubLayer = int(numOfTimeSteps / numOfSubLayers) 
    def makeData(fileName='dataSet.csv'):
        points = np.zeros((numOfParticles, numOfSites))
        maxHeihgts = np.zeros((numOfTimeSteps, numOfSites))
        for sl in range(numOfSubLayers):
            if sl % 2 == 0:
                color = 1
            elif sl % 2 == 1:
                color = 2
            for t in range(sl * numOfTimeStepsInEachSubLayer, (sl+1) * numOfTimeStepsInEachSubLayer):
                if t != 0:
                    maxHeihgts[t] = maxHeihgts[t-1]
                for i in range(timeSteps[t]):
                    x = np.random.randint(numOfSites)
                    left = maxHeihgts[t][x-1]
                    right = maxHeihgts[t][(x+1) % numOfSites]
                    mid = maxHeihgts[t][x] + 1
                    newHeight = int(np.max([left, right, mid]))
                    maxHeihgts[t][x] = newHeight
                    points[newHeight][x] = color
        finalPoints = np.zeros((int(np.max(maxHeihgts))+1, numOfSites))
        for i in range(int(np.max(maxHeihgts))+1):
            finalPoints[i][:] = points[i][:]
        if isShowLayer == True:
            return finalPoints
        if isShowPlotsOfSTD == True:
            for i in range(1, numOfTimeSteps):
                maxHeihgts[i] = maxHeihgts[i] + maxHeihgts[i-1]
        np.savetxt(fileName, maxHeihgts, delimiter=',')
    def makeEnsembles():
        filesNames = []
        for i in range(numOfEnsembles):
            filesNames.append(f'dataSet{i}Q2L={numOfSites}.csv')
            makeData(fileName=filesNames[i])
        return filesNames
    def showLayer():
        finalPoints = makeData()
        fig , ax = plt.subplots(figsize=(10, 7))
        ax.imshow(finalPoints, interpolation = 'none', cmap = cm.hot)
        ax.invert_yaxis()
        ax.set_xlabel(f'The size of the lattice is {numOfSites}')
        ax.set_ylabel(f'Number of total particles is {numOfParticles}')
        ax.set_title(f'->Ballistic Deposition (KENAR NESHAST)<-')
        plt.show()
    if isShowLayer == True:
        showLayer()
    def dataAnalysers():
        # filesNames = makeEnsembles()
        filesNames = [f'dataSet{i}Q2L={numOfSites}.csv' for i in range(numOfEnsembles)]
        STDValuesOfEnsembles = np.zeros((numOfEnsembles, numOfTimeSteps))
        for i in range(numOfEnsembles):
            dataSet = np.loadtxt(filesNames[i], delimiter = ',')
            for k in range(numOfTimeSteps):
                STDValuesOfEnsembles[i][k] = np.std(dataSet[k])
        averageSTDValues = STDValuesOfEnsembles[0]
        for i in range(1, numOfEnsembles):
            averageSTDValues = averageSTDValues + STDValuesOfEnsembles[i]
        averageSTDValues = averageSTDValues/numOfEnsembles
        x1 = np.log(timeSteps)[10:30]
        y1 = np.log(averageSTDValues)[10:30]
        x2 = np.log(timeSteps)[35:40]
        y2 = np.log(averageSTDValues)[35:40]
        slope1, intercept1, r, p, se = scs.linregress(x1, y1)
        slope2, intercept2, r, p, se = scs.linregress(x2, y2)
        plt.plot(x1, slope1*x1+intercept1, color = 'r',
                label='fit with equation: {:.3f} * ln(x) + {:.3f}'.format(slope1, intercept1))
        plt.plot(x2, slope2*x2+intercept2, color='g',
                label='fit with equation: {:.3f} * ln(x) + {:.3f}'.format(slope2, intercept2))
        plt.scatter(np.log(timeSteps), np.log(averageSTDValues), s=1)
        plt.title(f'size of lattice is {numOfSites} and number of ensembles is {numOfEnsembles}')
        plt.xlabel(f'ln(t) for times steps 1000 * 1.2 ** x in range({numOfTimeSteps})')
        plt.ylabel('ln(STD) per ln(t)')
        plt.axis('equal')
        plt.legend()
        plt.show()
    if isShowPlotsOfSTD == True:
        dataAnalysers()





# ballisticSideDeposition(isShowLayer=True)
# ballisticSideDeposition(isShowPlotsOfSTD=True)

x = np.log([50, 100, 150, 200, 250])
y = [2.956, 1.749, 3.674, 3.502, 5.116]
slope, intercept, s, p, r = scs.linregress(x, y)
plt.scatter(x, y)
plt.plot(x, slope*x+intercept, label='fit with equation {:.3f} ln(L) + {:.3f}'.format(slope, intercept))
plt.xlabel('ln(L)')
plt.ylabel('ln(Ws)')
plt.legend()
plt.show()
