from types import LambdaType
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as scs
# np.random.seed(0)
def ballisticDeposition(numOfParticles=10**5, numOfTimeSteps=40,
                        numOfSites=50, numOfEnsembles=20, numOfSubLayers=4,
                        isShowLayer = False, isShowPlotsOfSTD = False,
                        colors = ['blue', 'orange']):
    if isShowPlotsOfSTD == True:
        #to have a better plot for STD per time we use logaritmic time steps
        timeSteps = [int(1000 * 1.2 ** x) for x in range(numOfTimeSteps)]
    elif isShowLayer == True:
        #to have a better shape of layer we use uniform time steps
        timeSteps = [int(numOfParticles/numOfTimeSteps) for x in range(numOfTimeSteps)]
    #This function makes heights per time stpes
    def makeData(fileName='dataSet.csv'):
        heights = np.zeros((numOfTimeSteps, numOfSites))
        for i in range(numOfTimeSteps):
            for j in range(timeSteps[i]):
                x = np.random.randint(numOfSites)
                left = heights[i][x-1]
                right = heights[i][(x+1) % numOfSites]
                mid = heights[i][x]
                minHeight = np.min([left, right, mid])
                if minHeight == mid:
                    heights[i][x] += 1
                elif left == right:
                    heights[i][np.random.choice([x-1, (x+1) % numOfSites])] += 1
                elif minHeight == left:
                    heights[i][x-1] += 1
                elif minHeight == right:
                    heights[i][(x+1) % numOfSites] += 1
        if isShowPlotsOfSTD == True:
            for i in range(numOfTimeSteps):
                if i != 0:
                    heights[i] = heights[i] + heights[i-1]
            np.savetxt(fileName, heights, delimiter = ',')

        elif isShowLayer == True:
            np.savetxt(fileName, heights, delimiter = ',')
            return fileName

        return fileName
    #This function makes ensembles
    def makeEnsembles():
        filesNames = []
        for i in range(numOfEnsembles):
            filesNames.append(f'dataSet{i}Q1L={numOfSites}.csv')
            makeData(fileName=filesNames[i])
        return filesNames
    #This function will show the layer shape
    def showLayer():
        fileName = makeData()
        dataSet = np.loadtxt(fileName, delimiter = ',')
        sites = [x for x in range(numOfSites)]
        numOfTimeStepsInEachSubLayer = int(numOfTimeSteps / numOfSubLayers)
        bottomSetter = np.zeros(numOfSites)
        for sl in range(numOfSubLayers):
            if sl % 2 == 0:
                color = colors[0]
            else:
                color = colors[1]
            for i in range(sl * numOfTimeStepsInEachSubLayer, (sl+1) * numOfTimeStepsInEachSubLayer):
                plt.subplots_adjust()
                plt.bar(sites, dataSet[i], bottom = bottomSetter, width = 1.0, color = color)
                bottomSetter = bottomSetter + dataSet[i]
        plt.title(f'The shape of layer (for {fileName})')
        plt.xlabel(f'Size of lattice is {numOfSites}')
        plt.ylabel(f'Nmuber of total particles is {numOfParticles}')
        plt.show()
    #This function calculate the STD values.
    def dataAnalysers():
        # filesNames = makeEnsembles()
        filesNames = [f'dataSet{i}Q1L={numOfSites}.csv' for i in range(numOfEnsembles)]
        STDValuesOfEnsembles = np.zeros((numOfEnsembles, numOfTimeSteps))
        for i in range(numOfEnsembles):
            dataSet = np.loadtxt(filesNames[i], delimiter = ',')
            for k in range(numOfTimeSteps):
                STDValuesOfEnsembles[i][k] = np.std(dataSet[k])
        averageSTDValues = STDValuesOfEnsembles[0]
        for i in range(1, numOfEnsembles):
            averageSTDValues = averageSTDValues + STDValuesOfEnsembles[i]
        averageSTDValues = averageSTDValues/numOfEnsembles
        '''
        for dataSet{i}Q1L=200.csv files:
        x1 = np.log(timeSteps)[14:39]
        y1 = np.log(averageSTDValues)[14:39]
        x2 = np.log(timeSteps)[40:45]
        y2 = np.log(averageSTDValues)[40:45]
        '''
        x1 = np.log(timeSteps)[10:35]
        y1 = np.log(averageSTDValues)[10:35]
        x2 = np.log(timeSteps)[37:40]
        y2 = np.log(averageSTDValues)[37:40]


        plt.scatter(np.log(timeSteps), np.log(averageSTDValues), s=1, color='b', label='data')
        slope1, intercept1, r, p, se = scs.linregress(x1, y1)
        slope2, intercept2, r, p, se = scs.linregress(x2, y2)
        plt.plot(x1, slope1*x1+intercept1, color = 'r',
                label='fit with equation: {:.3f} * ln(x) + {:.3f}'.format(slope1, intercept1))
        plt.plot(x2, slope2*x2+intercept2, color='g',
                label='fit with equation: {:.3f} * ln(x) + {:.3f}'.format(slope2, intercept2))

        
        plt.title(f'size of lattice is {numOfSites} and number of ensembles is {numOfEnsembles}')
        plt.xlabel(f'ln(t) for times steps 1000 * 1.2 ** x in range({numOfTimeSteps})')
        plt.ylabel('ln(STD) per ln(t)')
        plt.legend()
        plt.show()
    if isShowPlotsOfSTD == True:
        dataAnalysers()
    elif isShowLayer == True:
        showLayer()

# ballisticDeposition(isShowLayer=True)
# ballisticDeposition(isShowPlotsOfSTD=True)
