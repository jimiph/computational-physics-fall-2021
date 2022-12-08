import numpy as np
import matplotlib.pyplot as plt



def randomWalkerWithTraps(numOfWalkers = 20000, numOfSites = 20, pRightValue = 0.3):
    def makeLifeTimePerInitialPosition(pRightValue = pRightValue):
        lifeTimePerInitianlPosition = np.zeros(numOfSites)
        for x0 in range(numOfSites):
            x = x0
            while True:
                if x == 0:
                    break
                elif x == numOfSites - 1:
                    break
                else:
                    x = x + np.random.choice([1, -1], p=[pRightValue, 1-pRightValue])
                    lifeTimePerInitianlPosition[x0] += 1
        return lifeTimePerInitianlPosition
                

    def analyser():
        average = np.zeros(numOfSites)
        for i in range(numOfWalkers):
            average = average + makeLifeTimePerInitialPosition()
        average = average / numOfWalkers
        np.save(f'meanLifeSpanPerInitialPositionPRight={pRightValue}Q6.npy', average)
        positions = [x for x in range(numOfSites)]
        averageLifeSpan = np.average(average)
        print(f'average life-span of the random walker over all initial positions is: {averageLifeSpan}')
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.plot(positions, average, 'o-')
        ax.set_title(f'mean life-span for {numOfWalkers} random walkers per initial position. Right going P={pRightValue}')
        ax.set_xlabel('x0')
        ax.set_ylabel('mean life-span <T>')
        plt.show()

    analyser()


randomWalkerWithTraps(numOfSites=10, pRightValue=0.5)



# numOfSites = 5
# pRightValue = 0.5

# lifeTimePerInitianlPosition = np.zeros(numOfSites)
# for x0 in range(numOfSites):
#     x = x0
#     while True:
#         if x == 0:
#             break
#         elif x == numOfSites - 1:
#             break
#         else:
#             x = x + np.random.choice([1, -1], p=[pRightValue, 1-pRightValue])
#             lifeTimePerInitianlPosition[x0] += 1
# print(lifeTimePerInitianlPosition)

