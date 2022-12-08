from genericpath import samefile
import imp
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt


"""
fucntional algorithem
"""

# scaling by factor 1/3
def f1(vertices):
    newVertices = []
    for vertex in vertices:
        x = vertex[0]
        y = vertex[1]
        xNew = (1/3) * x
        yNew = (1/3) * y
        newVertex = [xNew, yNew]
        newVertices.append(newVertex)
    return newVertices

# scaling by factor 1/3 and clockwise rotaion with 60 degree. Then translation in 
# +x direction by 1/3
def f2(vertices):
    newVertices = []
    for vertex in vertices:
        x = vertex[0]
        y = vertex[1]
        xNew = (1/3) + (1/6)*(x - sqrt(3)*y)
        yNew = (1/6)*(sqrt(3)*x + y)
        newVertex = [xNew, yNew]
        newVertices.append(newVertex)
    return newVertices

# scaling by factor 1/3 and clockwise rotaion with 60 degree. Then translation in 
# +x direction by 1/2 and in direction +y by sqrt(3)/6
def f3(vertices):
    newVertices = []
    for vertex in vertices:
        x = vertex[0]
        y = vertex[1]
        xNew = (1/2) + (1/6) * (x+sqrt(3)*y)
        yNew = (sqrt(3)/6) + (1/6) * (-sqrt(3)*x+y)
        newVertex = [xNew, yNew]
        newVertices.append(newVertex)
    return newVertices

# scaling by factor 1/3 and translatoin by 2/3 in +x direction
def f4(vertices):
    newVertices = []
    for vertex in vertices:
        x = vertex[0]
        y = vertex[1]
        xNew = (2/3) + (1/3) * x
        yNew = (1/3) * y
        newVertex = [xNew, yNew]
        newVertices.append(newVertex)
    return newVertices

# n in number of iterations
def makeKokhFractal(n):
    initialVertices = [[0, 0], [1, 0]]
    functions = [f1, f2, f3, f4]

    functions = functions
    vertices = initialVertices
    for i in range(0, n):
        newVertices = []
        for func in functions:
            newVertices = newVertices + func(vertices)
        vertices.clear()
        vertices = newVertices
    for x in range(1, len(vertices)):
        xValues = [vertices[x-1][0], vertices[x][0]]
        yValues = [vertices[x-1][1], vertices[x][1]]
        plt.plot(xValues, yValues, color="blue")
    plt.title(f"Kokh's Fractal with {n} iterations")
    plt.axis('equal')
    plt.show()
# makeKokhFractal(4)

"""
objective algorithm
"""



