import numpy as np
import matplotlib.pyplot as plt


"""
f1 : to make the yellow part of the figure of 2.6
f2: to make the red part of the figure of 2.6
f3: to make the dark blue part of the figure of 2.6
f4: to make the light blue part of the figure of 2.6

"""

#This function make a scaling on y by 0.16
def f1(point):
    y = point[1]
    point[0] = 0
    point[1] = 0.16 * y
    return point

# This function makes a scaling and rotation and translation
def f2(point):
    x = point[0]
    y = point[1]
    point[0] = 0.85 * x + 0.04 * y
    point[1] = -0.04 * x + 0.85 * y + 1.6
    return point

# This function makes a scaling and rotation and translation
def f3(point):
    x = point[0]
    y = point[1]
    point[0] = 0.2 * x - 0.26 * y
    point[1] = 0.23 * x + 0.22 * y + 1.6
    return point

# This function makes a scaling and rotation and translation
def f4(point):
    x = point[0]
    y = point[1]
    point[0] = -0.15 * x + 0.28 * y
    point[1] = 0.26 * x + 0.24 * y + 0.44
    return point

#To make the fern more attractive and also the same the picture in the book, i used
# this probability distribution. I found this numbers by searching on the web.
# This funcion gets a random number and with a specific probability distribusion
# return another random number which belongs to the set {0, 1, 2, 3} 
def probabilityDisturb(randNum):
    if randNum < 1:
        newRand = 0
    elif 1 <= randNum < 86:
        newRand = 1
    elif 86 <= randNum < 93:
        newRand = 2
    else:
        newRand = 3
    return newRand

#This function operates on a given point to make a new random point.
def makeRandomTransformedPoint(oldPoint, levelOfFractal):
    functions = [f1, f2, f3, f4]
    for i in range(0, levelOfFractal):
        randomNum = probabilityDisturb(np.random.randint(0, 100))
        newPoint = functions[randomNum](oldPoint)
        oldPoint = newPoint
    return newPoint

#This function operates the algorithm on each point to make the fern.
def makeRandomFern(numOfPoints, levelOfFractal):
    points = []
    for i in range(0, numOfPoints):
        oldPoint = [-1 + 2 * np.random.random(), -1 + 2 * np.random.random()]
        newPoint = makeRandomTransformedPoint(oldPoint, levelOfFractal)
        points.append(newPoint)
    xValuesOfPoints = []
    yValuesOfPoints = []
    for i in range(0, len(points)):
        xValuesOfPoints.append(points[i][0])
        yValuesOfPoints.append(points[i][1])
    plt.figure(figsize=(10, 10))
    plt.title(f"number of points is {numOfPoints} and number of iteration is {levelOfFractal}")
    plt.scatter(xValuesOfPoints, yValuesOfPoints, s=0.1, color='green')
    plt.axis("equal")
    plt.show()


makeRandomFern(10**5, 18)
