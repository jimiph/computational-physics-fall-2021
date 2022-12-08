from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt


"""
about vertices:
                *vertice(2)



    *vertice(0)              *vertice(1)

about newTriangles:
                *newTriangle(2)



    *newTriangle(0)              *newTriangle(1)

"""


# this funcion gets 3 vertices and return a newTriangle(0)
def newTriangleMaker0(vertices):
    side = vertices[1][0] - vertices[0][0]
    deltaX = side/2
    deltaY = (sqrt(3)/2) * side
    newVertex10 = vertices[1][0] - deltaX
    newVertex11 = vertices[1][1]
    newVertex20 = vertices[2][0] - deltaX/2
    newVertex21 = vertices[2][1] - deltaY/2
    newVertices = [vertices[0], [newVertex10,
                                 newVertex11], [newVertex20, newVertex21]]
    return newVertices


# this funcion gets 3 vertices and return a newTriangle(1)
def newTriangleMaker1(vertices):
    side = vertices[1][0] - vertices[0][0]
    deltaX = side/2
    deltaY = (sqrt(3)/2) * side
    newVertex20 = vertices[2][0] + deltaX/2
    newVertex21 = vertices[2][1] - deltaY/2
    newVertex00 = vertices[0][0] + deltaX
    newVertex01 = vertices[0][1]
    newVertices = [[newVertex00, newVertex01],
                   vertices[1], [newVertex20, newVertex21]]
    return newVertices


# this funcion gets 3 vertices and return a newTriangle(2)
def newTriangleMaker2(vertices):
    side = vertices[1][0] - vertices[0][0]
    deltaX = side/2
    deltaY = (sqrt(3)/2) * side
    newVertex10 = vertices[1][0] - deltaX/2
    newVertex11 = vertices[1][1] + deltaY/2
    newVertex00 = vertices[0][0] + deltaX/2
    newVertex01 = vertices[0][1] + deltaY/2
    newVertices = [[newVertex00, newVertex01], [
        newVertex10, newVertex11], vertices[2]]
    return newVertices


# this funcion gets a traingle and make 3 new triangles
def newTriangleMaker(triangle=np.array([[0, 0], [1, 0], [0.5, sqrt(3)/2]])):
    newTriangle0 = newTriangleMaker0(triangle)
    newTriangle1 = newTriangleMaker1(triangle)
    newTriangle2 = newTriangleMaker2(triangle)
    newTriangles = [newTriangle0, newTriangle1, newTriangle2]
    return newTriangles


# order means number of iterations
# finally this function will creat the fractal
def showFractal(order):
    currentTriangle = newTriangleMaker()
    for i in range(0, order-1):
        new = []
        for triangle in currentTriangle:
            new = new + newTriangleMaker(triangle)
        currentTriangle = new
    currentTriangle = np.array(currentTriangle)
    ax = plt.axes()
    for triangle in currentTriangle:
        ax.add_patch(Polygon(triangle))
    plt.axis("equal")
    plt.title(f"Sierpinski's Triangle with {order} iteration")
    plt.show()


showFractal(10)
