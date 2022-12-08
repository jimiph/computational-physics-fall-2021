import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss


# np.random.seed(0)
def massDimensionOfClusterOfPercolation(probability=0.57, lengthOfLattice=100, numOfRuns=200,
                                        isShowCluster=False, isAnalyse=False):
 
    #This function return 1 or -1 with a uniform distubution.
    def randomTurningOn(p=probability):
        rand = np.random.uniform()
        # print(rand)
        if rand < p:
            return 1
        else:
            return -1
    
    #This function grow a cluster from a given point and return it as 2D numpy array.
    def makeCluster(startPoint=(int(lengthOfLattice/2), int(lengthOfLattice/2)),
     probability=probability, lengthOfLattice=lengthOfLattice):

        binaryMatrix = np.zeros((lengthOfLattice+2, lengthOfLattice+2), dtype=int)
        binaryMatrix[startPoint[0]][startPoint[1]] = 1
        binaryMatrix[0][:] = -1
        binaryMatrix[lengthOfLattice+1][:] = -1
        for x in range(lengthOfLattice+1):
            binaryMatrix[x][0] = -1
            binaryMatrix[x][lengthOfLattice+1] = -1

        #By This function we turn on or block the neighbors of a light cell.
        def checkNeighbors(x, y):
            cells = []
            if binaryMatrix[x][y] == 1:
                if binaryMatrix[x-1][y] == 0:
                    r = randomTurningOn(probability)
                    binaryMatrix[x-1][y] = r
                    if r == 1:
                        cells.append((x-1, y))
                if binaryMatrix[x+1][y] == 0:
                    r = randomTurningOn(probability)
                    binaryMatrix[x+1][y] = r
                    if r == 1:
                        cells.append((x+1, y))
                if binaryMatrix[x][y-1] == 0:
                    r = randomTurningOn(probability)
                    binaryMatrix[x][y-1] = r
                    if r == 1:
                        cells.append((x, y-1))
                if binaryMatrix[x][y+1] == 0:
                    r = randomTurningOn(probability)
                    binaryMatrix[x][y+1] = r
                    if r == 1:
                        cells.append((x, y+1))      
            return cells      


        boundCells = [(startPoint[0], startPoint[1])]
        while len(boundCells) != 0:
            currentBoundCells = []
            for cell in boundCells:
                x0 = cell[0]
                y0 = cell[1]
                currentBoundCells = currentBoundCells + checkNeighbors(x0, y0)
            boundCells = currentBoundCells

        finalCluster = np.zeros((lengthOfLattice, lengthOfLattice))
        for x in range(1, lengthOfLattice+1):
            for y in range(1, lengthOfLattice+1):
                if binaryMatrix[x][y] == -1:
                    finalCluster[x-1][y-1] = 0
                else:
                    finalCluster[x-1][y-1] = binaryMatrix[x][y]
        return finalCluster


    #This function make data and plot the log(s) - log(gr)
    def analyse():
        pValues = [0.5, 0.55, 0.59]
        #This function find the coordinates of all light cells and return it as a list.
        def findCoordinates(matrix):
            coordinates = []
            for x in range(lengthOfLattice):
                for y in range(lengthOfLattice):
                    if matrix[x][y] == 1:
                        coordinates.append((x, y))
            return coordinates

        def calculateArea(matrix):
            return len(findCoordinates(matrix))
            
        def calculateGyrationRadius(matrix):
            coordinates = findCoordinates(matrix)
            xs = np.zeros(len(coordinates))
            ys = np.zeros(len(coordinates))
            for i in range(len(coordinates)):
                xs[i], ys[i] = coordinates[i]
            xmid, ymid = np.average(coordinates, axis=0)
            r = np.sqrt(np.average(np.square(xs-xmid)+np.square(ys-ymid)))
            return r

        for p in pValues:
            areas = np.zeros(numOfRuns)
            gyrationRadiuses = np.zeros(numOfRuns)
            for i in range(numOfRuns):
                cluster = makeCluster(probability=p)
                areas[i] = calculateArea(cluster)
                gyrationRadiuses[i] = calculateGyrationRadius(cluster)
            np.save(f'areasForP={p}Q2.npy', areas)
            np.save(f'gyrationRadiusesForP={p}Q2.npy', gyrationRadiuses)
            
            xx = np.array([gyrationRadiuses[i]  for i in range(numOfRuns) if gyrationRadiuses[i] != 0])
            yy = np.array([areas[j] for j in  [np.where(gyrationRadiuses==i)[0][0] for i in xx]])
            x = np.log(xx)
            y = np.log(yy)
            slope, intercept, r, pp, se = ss.linregress(x, y)
            plt.scatter(x, y, s=1)
            plt.plot(x, slope*x+intercept, color='r', label='fit with equation: {:.3} * ln(gyration radius) + {:0.3}'.format(slope, intercept))
            plt.title(f'p = {p}')
            plt.xlabel('ln(gyration radius)')
            plt.ylabel('ln(areas)')
            plt.legend()
            plt.show()



    if isShowCluster == True:
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(makeCluster())
        ax.invert_yaxis()
        ax.set_title(f'p={probability} and L={lengthOfLattice}')
        plt.show()

    if isAnalyse == True:
        analyse()
    



#To show the cluster:
# massDimensionOfClusterOfPercolation(isShowCluster=True)

#To show the plots:
# massDimensionOfClusterOfPercolation(isAnalyse=True)


