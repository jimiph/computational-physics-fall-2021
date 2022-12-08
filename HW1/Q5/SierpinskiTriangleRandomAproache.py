import numpy as np
import matplotlib.pyplot as plt


#This function contains a scaling transformation by 0.5 factor
def f1(point):
    point[0] = 0.5 * point[0]
    point[1] = 0.5 * point[1]
    return point

#This function contains a scaling transformation by 0.5 factor and 
# a translation in x direction by 0.5
def f2(point):
    point[0] = 0.5 * point[0] + 0.5
    point[1] = 0.5 * point[1]
    return point

#This function contains a scaling transformation by 0.5 factor and 
# a translation in x direction by 0.25 and a translation in y direction
# by 0.25 * sqrt(3)
def f3(point):
    point[0] = 0.5 * point[0] + 0.25
    point[1] = 0.5 * point[1] + 0.25 * np.sqrt(3)
    return point

#This funcion gets a point and randomly operates one of the 3 functions
# f1, f2 or f3 on that point for levelOfFractal iterations; then return 
# the new point.
def makeRandomTransformedPoint(oldPoint, levelOfFractal):
    functions = [f1, f2, f3]
    for i in range(0, levelOfFractal):
        randomNum = np.random.randint(0, 3)
        newPoint = functions[randomNum](oldPoint)
        oldPoint = newPoint
    return newPoint

#This function operates the algorithm for a customized number of random points.
def makeRandomSierpinski(levelOfFractal, numOfPoints):
    points = []
    for i in range(0, numOfPoints):
        oldPoint = [np.random.rand(), np.random.rand()]
        newPoint = makeRandomTransformedPoint(oldPoint, levelOfFractal)
        points.append(newPoint)
    xValuesOfPoints = []
    yValuesOfPoints = []
    for i in range(0, len(points)):
        xValuesOfPoints.append(points[i][0])
        yValuesOfPoints.append(points[i][1])
    # plt.axes("equal")
    plt.title(f"number of points is {numOfPoints} and number of iteraions is {levelOfFractal}")
    plt.scatter(xValuesOfPoints, yValuesOfPoints, s = 0.2)
    plt.show()


makeRandomSierpinski(4, 10**4)

















# """
# size = 10**4

# points = np.random.rand(2, size)
# # print(points)
# complexPoints = points[0] + points[1] * 1j

# theta = np.deg2rad(23)
# trans = np.exp(theta * 1j)
# newComplexPoints = trans * complexPoints

# transMatrix = np.array([[np.cos(theta), -np.sin(theta)], 
#                         [np.sin(theta), np.cos(theta)]])

# newPoints = np.matmul(transMatrix, points)
# newPoints = transMatrix @ points

# plt.figure(figsize=(10, 10))
# # plt.scatter(newComplexPoints.real, newComplexPoints.imag, s=0.25)
# plt.scatter(newPoints[0], newPoints[1], s=0.25)

# plt.axis("equal")
# plt.show()
# """
