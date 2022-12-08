import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats as sp

np.random.seed(2)
def ballisticDeposition(numOfParticles=10**6, numOfTimeSteps=100,
                        numOfSites=800, numOfSubLayers=4,
                        isShowLayer=False, isShowPlotOfWidths=False):
    timeSteps = [int(numOfParticles/numOfTimeSteps)
                 for x in range(numOfTimeSteps)]
    numOfTimeStepsInEachSubLayer = int(numOfTimeSteps / numOfSubLayers)

    def makeData(fileName='particles.csv'):
        particles = np.zeros((numOfParticles, numOfSites))
        particles[1][int(numOfSites/2)] = 1
        maxHeihgts = np.zeros((numOfTimeSteps, numOfSites))
        for sl in range(numOfSubLayers):
            if sl % 2 == 0 and isShowLayer == True:
                color = 1
            elif sl % 2 == 1 and isShowLayer == True:
                color = 2
            else:
                color = 1
            for t in range(sl * numOfTimeStepsInEachSubLayer, (sl+1) * numOfTimeStepsInEachSubLayer):
                if t != 0:
                    maxHeihgts[t] = maxHeihgts[t-1]
                for i in range(timeSteps[t]):
                    x = np.random.randint(numOfSites)
                    left = maxHeihgts[t][x-1]
                    right = maxHeihgts[t][(x+1) % numOfSites]
                    mid = maxHeihgts[t][x] + 1
                    newHeight = int(np.max([left, right, mid]))
                    if newHeight > 1 and (particles[newHeight-1][x] != 0 or particles[newHeight][x-1] != 0
                                          or particles[newHeight][(x+1) % numOfSites] != 0):
                        maxHeihgts[t][x] = newHeight
                        particles[newHeight][x] = color
                    elif (particles[newHeight][x-1] != 0 or particles[newHeight][(x+1) % numOfSites] != 0):
                        maxHeihgts[t][x] = newHeight
                        particles[newHeight][x] = color
        finalParticles = np.zeros((int(np.max(maxHeihgts))+1, numOfSites))
        for i in range(int(np.max(maxHeihgts))+1):
            finalParticles[i][:] = particles[i][:]
        np.savetxt(fileName, finalParticles, delimiter=',')
        return fileName

    def showLayer():
        fileName = makeData()
        particles = np.loadtxt(fileName, delimiter=',')
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.imshow(particles, interpolation='none', cmap=cm.hot)
        ax.invert_yaxis()
        ax.set_xlabel(f'The size of the lattice is {numOfSites}')
        ax.set_ylabel(f'Number of total particles is {numOfParticles}')
        ax.set_title(f'->Ballistic Deposition (KENAR NESHAST)<-')
        plt.show()

    def dataAnalysers():
        fileName = makeData()
        particles = np.loadtxt(fileName, delimiter=',')
        widths = np.zeros((len(particles)))
        particlesTillEachHeight = np.zeros((len(particles)))
        for h in range(1, len(particles)):
            numOfLeftZeros = 0
            numOfRightZeros = 0
            x = 0
            while particles[h][x] == 0 :
                numOfLeftZeros += 1
                x  = x + 1
            y = numOfSites - 1
            while particles[h][y] == 0 :
                numOfRightZeros += 1
                y  = y - 1
            widths[h] = numOfSites - (numOfRightZeros + numOfLeftZeros)
        for h in range(1, len(particles) + 1):
            x = np.zeros(h)
            for i in range(h):
                x[i] = widths[i]
            widths[h-1] = int(np.max(x))
        for h in range(len(particles)):
            particlesTillEachHeight[h] = np.sum(particles[h])
        for h in range(len(particles)):
            if h != 0:
                particlesTillEachHeight[h] = particlesTillEachHeight[h] + particlesTillEachHeight[h-1]
        
        '''
        In the following you can uncomment the commented code and comment other codes 
        to show you the linear plot; otherwise the log-log plot will be showen.
        '''
        slope, intercept, r, p, se = sp.linregress(np.log(particlesTillEachHeight)[566:628], 
                                                    np.log(widths)[566:628])
        # slope, intercept, r, p, se = sp.linregress((particlesTillEachHeight)[566:628], 
                                                    # (widths)[566:628])
        '''
        I used the saved data bottom to find best resion for fitting the regression line.
        The special numbers 566 and 628 are from here.
        '''
        # np.savetxt('widths.csv',np.log(widths[400:]), delimiter = ',' )
        # np.savetxt('particlesTillEachHeight.csv', np.log(particlesTillEachHeight[400:]), delimiter = ',')
        plt.scatter(np.log(particlesTillEachHeight[400:]), np.log(widths[400:]), label='particles', s=1)
        # plt.scatter((particlesTillEachHeight[400:700]), (widths[400:700]), label='particles', s=1)
        plt.plot(np.log(particlesTillEachHeight[566:628]), 
                slope * np.log(particlesTillEachHeight[566:628]) + intercept, color='r',
                label = 'fit with equation: {:.3} * ln(x) + {:0.3}'.format(slope, intercept))
        # plt.plot((particlesTillEachHeight[566:628]), 
                # slope * (particlesTillEachHeight[566:628]) + intercept, color='r',
                # label = 'fit with equation: {:.3} * x + {:0.3}'.format(slope, intercept))
        plt.legend()
        plt.xlabel(f'ln(time steps: {int(numOfParticles/numOfTimeSteps)})')
        # plt.xlabel(f'time steps: {int(numOfParticles/numOfTimeSteps)}')
        plt.ylabel('ln(correlation length)')
        # plt.ylabel('correlation length')
        plt.title('correlation length')
        plt.grid()
        plt.show()
    if isShowLayer == True:
        showLayer()
    if isShowPlotOfWidths == True:
        dataAnalysers()

        

#If you wnat to show the plot of the length of correlation per time use this:
# ballisticDeposition(isShowPlotOfWidths=True)
#If you want tot show the tree per lattice sites use this:
# ballisticDeposition(isShowLayer=True)

