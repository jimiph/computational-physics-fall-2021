from numpy import random
import matplotlib.pyplot as plt
import numpy as np



# this function makes Khayyam's Triangle
def makeKhayyamNumbers(numOfRows):
    khayyamNumbers = [[1]]
    for i in range(0, numOfRows-1):
        old = khayyamNumbers[i]
        # print(old)
        new = []
        for j in range(0, len(old)-1):
            newNum = old[j] + old[j+1]
            new = new + [newNum]
        new.append(1)
        # new.append(1)
        new.insert(0, 1)
        # new.insert(0, 1)
        khayyamNumbers =  khayyamNumbers + [new]
    return khayyamNumbers

# This function make an numOfRows * numOfRows array. Because I want to use the
# imshow() function I need to make this array. Also I dublicate each array because
# I want to make a equilateral triangle as my output. This idea is shown below.
"""
for numOfRows = 3:
    0 0 1 1 0 0
    0 0 1 1 0 0
    0 1 1 1 1 0
    0 1 1 1 1 0
    0 1 2 2 1 0
    0 1 2 2 1 0
"""
def makeArrayOfKhayyamNumbers(numOfRows):

    khayyamNumbers = makeKhayyamNumbers(numOfRows)
    newKhayyamNumbers = []
    n = len(khayyamNumbers)
    for i in range(0, n):
        new = []
        for j in range(0, len(khayyamNumbers[i])):
            new.append(khayyamNumbers[i][j])
            new.append(khayyamNumbers[i][j])
        newKhayyamNumbers.append(new)

    for i in range(0, numOfRows):
        for j in range(0, (numOfRows - i -1)):
                newKhayyamNumbers[i].append(0)
                newKhayyamNumbers[i].insert(0, 0)
    m = len(newKhayyamNumbers)
    for i in range(0, 2*m, 2):
        newKhayyamNumbers.insert(i, newKhayyamNumbers[i])
    return newKhayyamNumbers

# This function: if a number is odd I put 1 instead of it and
# if it is even I put 0 instead of it. Then I use the function 
# imshow() to display this triangle. 
def showKhayyamTriangle(numOfRows):
    ArrayOfKhayyamNumbers = makeArrayOfKhayyamNumbers(numOfRows)
    for i in range(0, 2*numOfRows):
        for j in range(0, 2*numOfRows):
            ArrayOfKhayyamNumbers[i][j] = ArrayOfKhayyamNumbers[i][j] % 2
    plt.title(f"Sierpinski's Triangle using Khayyam's Triangle with {numOfRows}Rows")
    plt.imshow(ArrayOfKhayyamNumbers)
    plt.show()


showKhayyamTriangle(100)