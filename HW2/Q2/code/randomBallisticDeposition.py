import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.function_base import average



def makeLayers(isShowLayer = False, numTotalParticles = 10**5, numSites = 200
                , numLayers = 4, isCalculateHeights = False, isSqueredHeightValues = False):
    #To change the color of the sublayers
    colors = ['blue', 'orange']
    #Discrete position
    sites = [x for x in range(numSites)]
    #Number of total prticles
    numTotalParticlesInLayer = int(numTotalParticles / numLayers)
    heights = np.zeros((numLayers, numSites))
    #To set the bottom of the bars; beacuse we want to make sublayers with
    #differen colors.
    bottomSetter = np.zeros(numSites)
    for layer in range(numLayers):
        
        for particle in range(numTotalParticlesInLayer):
            randomPosition = np.random.randint(numSites)
            #Randomly increas the hight of a site
            heights[layer][randomPosition] += 1                    

        if isShowLayer == True:
            if layer % 2 == 0:
                color = colors[0]
            else:
                color = colors[1]
            plt.subplots_adjust()
            plt.bar(sites, heights[layer], bottom = bottomSetter, width = 1.0, color = color)
            bottomSetter = bottomSetter + heights[layer]
    if isShowLayer == True:
        plt.title(f'Random Ballistic Deposition for {numTotalParticles} Particles')
        plt.show()
    totalHeights = np.zeros((numSites))
    if isCalculateHeights == True: 
        for i in range(numLayers):
            totalHeights = totalHeights + heights[i]
        if isSqueredHeightValues == False:
            return totalHeights
        elif isSqueredHeightValues == True:
            return np.square(totalHeights)

def dataCalculations(isShowAveragePlot = False, isShowSTDPlot = False):
    heightAverages = []
    squaredHeightAverages = []
    numParticlesInEachStep = []
    #I calculate the height averages for each 1000 particles
    numTotalParticles = 10 ** 5
    if isShowAveragePlot == True:
        for i in range(int(numTotalParticles / 1000)):
            numParticlesInEachStep.append(i * 1000)
            heightAverages.append(np.average(makeLayers(numTotalParticles=i * 1000, 
                                                        isCalculateHeights=True)))
        plt.title(f'This is average height for each 1000 particles.')
        plt.plot(numParticlesInEachStep, heightAverages)
        plt.show()
    elif isShowSTDPlot == True:
        for i in range(int(numTotalParticles / 1000)):
            numParticlesInEachStep.append(i * 1000)
            squaredHeightAverages.append(np.average(makeLayers(numTotalParticles=i * 1000, 
                                                        isCalculateHeights=True,
                                                        isSqueredHeightValues=True)))
        for i in range(int(numTotalParticles / 1000)):
            heightAverages.append(np.average(makeLayers(numTotalParticles=i * 1000, 
                                                        isCalculateHeights=True)))
        stdValues = np.sqrt(squaredHeightAverages - np.square(heightAverages))
        plt.title(f'This is std for each 1000 particles.Time is from 0 to {numTotalParticles}')
        plt.xlabel('log(t)')
        plt.ylabel('log(w(t))')
        plt.scatter(np.log(stdValues), np.log(numParticlesInEachStep), s=1)
        plt.show()


#to show the layer call this function:
#makeLayers(isShowLayer=True)

#to show the average heights plot, call this funcoin:
# dataCalculations(isShowAveragePlot=True)

#to show the log-log STD plot, call this funcoin:
# dataCalculations(isShowSTDPlot=True)
