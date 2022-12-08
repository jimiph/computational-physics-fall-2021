import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm

def showJulia(c):
    #This funcion map z to z^2 + c
    def fZ(z):
        return z**2 + c
    #Set the size of the figure
    plt.figure(figsize = (10, 7))
    #Set the Interval of x values and y values 
    xMin, xMax = -2, 2
    yMin, yMax = -2, 2
    width = xMax - xMin
    height = yMax - yMin
    #Set the number of points or pixels in each dimension
    numOfXs = 500
    numOfYs = 500
    #Set the maximum of iterations
    maxIteratoin = 1000
    #To check that a point is divegint or convergin  we define a maximum values
    #of norm or the radius of circle divergence. 
    maxNorm = 0.5 * (1 + np.sqrt(1 + 4 * abs(c)))
    julia = np.zeros((numOfXs, numOfYs))
    for nX in range(0, numOfXs):
        for nY in range(0, numOfYs):
            currentIteration = 0
            z = nX / numOfXs * width + xMin + 1j * (nY / numOfYs * height + yMin)
            #Mapping each point in complex plane again and again!
            while abs(z) < maxNorm and currentIteration < maxIteratoin:
                z = fZ(z)
                currentIteration += 1
            #We assign a values to each point in the complex plane. We use
            #the currentIteration to do this. Then we will map a color to each
            #point dependent to this value.
            julia[-nY, nX] = currentIteration + 1
    plt.title(f'Julia fractal for c={c} and maximum iteration={maxIteratoin},' + 
                f'\nalso number of pixels is {numOfXs} * {numOfYs}')
    plt.imshow(julia, interpolation='nearest', cmap=cm.hot,
                norm = LogNorm())
    plt.show()




showJulia(complex(-0.8, +0.16))

