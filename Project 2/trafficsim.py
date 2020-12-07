import random as rnd
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl

def getFlowrate(carVelocities, roadLength):
    """Takes array carVelocities and float or int roadLength. Returns the
    flowrate as the sum of all values in carVelocities devided by roadLength."""
    sum = 0
    for velocity in carVelocities:
        sum += velocity
    flowrate = sum/roadLength
    return flowrate

def randomizeCarPositions(nCars, roadLength):
    """This function takes integers nCars and roadLength. It then returns a
    sorted array of unique random integers from 0 to roadLength-1. The array is
    of length nCars."""
    carPositions = [x for x in range(roadLength)]
    while len(carPositions) > nCars:
        randomIndex = int(len(carPositions)*rnd.random()) % len(carPositions)
        del carPositions[randomIndex]
    return carPositions

def randomizeCarVelocities(nCars, vMax):
    """This function takes integers nCars, and vMax. It then returns an array
    of of random integers from 0 to vMax. The array is of length nCars.
    No integers in the returned array are the same."""
    carVelocities = []
    for i in range(nCars):
        randomVelocity = int((vMax+1)*rnd.random()) % (vMax+1)
        carVelocities.append(randomVelocity)
    return carVelocities

def traffic(nCars, carVelocities, vMax, carPositions, roadLength, pBrake):
    """Simulating the cellular automaton traffic model
    with periodic boundary conditions."""
    # accelerate on step towards max speed
    for i in range(nCars):
        if carVelocities[i] < vMax:
            carVelocities[i] += 1
    # brake if close to next car
    for i in range(nCars):
        distance = (carPositions[(i+1)%nCars] - carPositions[i]) % roadLength
        if carVelocities[i] >= distance:
            carVelocities[i] = distance - 1
    # random braking
    for i in range(nCars):
        if rnd.random() <= pBrake and carVelocities[i] > 0 :
            if carVelocities[i] > 0:
                carVelocities[i] -= 1
    # move cars forwards, the road is a circle
    for i in range(nCars):
        carPositions[i] += carVelocities[i]
        carPositions[i] %= roadLength


