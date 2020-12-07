
# -------------------- Import libraries

import random as rnd
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
from matplotlib import animation
from trafficsim import *

# -------------------- Calculation

def addToArrays(t):
    while t < tStop:
        t += 1
        traffic(nCars, carVelocities, vMax,
            carPositions, roadLength, pBrake)
        # store data
        tArray.append(t)
        carPositionsDict[t] = carPositions[:]

# iteration for positions vs time, nCars = 10
roadLength = 50
vMax = 2
pBrake = 0.5
tStop = 100      # stop calculating at this time

nCars = 10
carPositions = [x for x in range(nCars)]
carVelocities = [0 for x in range(nCars)]
t = 0
tArray = []
carPositionsDict = {}       # keys are t in tArray

addToArrays(t)
# calculate flow rate for fundamental diagram

## variables
nCarsArray = [0, 2, 5, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50]
nSimulations = 1000      # number of simulations for each nCars
roadLength = 50
vMax = 2
pBrake = 0.5
tStop = 100             # calculate flowrate at this time, after initial transition

## storage
flowrateArray = []      # will contain (nCars number of) arrays (of length nSimulations) corresponding to densityArray
densityArray =  []

for nCars in nCarsArray:
    flowrateArray_nCars = []
    for simulation in range(nSimulations):
        t = 0
        carPositions = randomizeCarPositions(nCars, roadLength)
        carVelocities = randomizeCarVelocities(nCars, vMax)
        while t < tStop:
            t += 1
            traffic(nCars, carVelocities, vMax,
                    carPositions, roadLength, pBrake)
        # store data
        flowrateArray_nCars.append(getFlowrate(carVelocities, roadLength))
    densityArray.append(nCars/roadLength)
    flowrateArray.append(flowrateArray_nCars[:])

# -------------------- Plot

def plotTraffic():
    plt.figure()
    colors = pl.cm.jet(np.linspace(0.1,0.9,10))
    plt.title('Traffic model - 10 cars, Road length - %d'%(roadLength))
    plt.xlabel('Time', fontsize = 12)
    plt.ylabel('Position', fontsize = 12)

    for i in range(len(tArray)):
        t = tArray[i]
        carPositions = carPositionsDict[t]
        for j in range(len(carPositions)):
            plt.plot(t, carPositions[j],
                marker = '.',markersize = 5,
                color = colors[j]
            )
    plt.show()
    
def plotDiagram():
    def fFitted(x,a,b,c):
        return (a * (x**b)) * np.exp(-c*x)

    meanFlowrateArray = [np.mean(array) for array in flowrateArray]

    popt, pcov = curve_fit(fFitted, densityArray, meanFlowrateArray)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    aRounded = round(a,2)
    bRounded = round(b,2)
    cRounded = round(c,2)

    xArray = np.linspace(min(densityArray), max(densityArray), 1001)

    plt.figure()
    plt.plot(xArray, [fFitted(x, a,b,c) for x in xArray],'r:', 
             label = r'$f(x) =  %s \, x^{%s} \, e^{- %s \, x }$'%(aRounded,bRounded,cRounded))
    plt.scatter(densityArray, meanFlowrateArray, color='tab:red', marker='o', label = 'Mean flowrates')
    plt.legend()
    plt.title('Fundamental diagram')
    plt.ylabel('Flow rate - $\Sigma V_{cars}/ L$', fontsize = 12)
    plt.xlabel('Density - $N_{cars} / L$', fontsize = 12)
    plt.show()

def output():
    plotTraffic()
    plotDiagram()

output()
