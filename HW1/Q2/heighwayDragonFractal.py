from turtle import right, left, forward, speed, exitonclick, hideturtle
from numpy import sqrt

#This function gets numer of iteration and size of the dragon.
def dragon(iteration=4, size=200, direction1=right, direction2=left):
    if iteration <= 0:
        forward(size)
        return
 
    size /= sqrt(2) # to scale to call the funcion again
    direction1(45) # rotaion by 45 degree
    dragon(iteration-1, size, right, left)
    direction2(90) # rotaion by 90 degree
    dragon(iteration-1, size, left, right)
    direction1(45) # rotaion by 45 degree
 
speed(0)
hideturtle()
dragon(8)
exitonclick() # click to exit