# -------------------- Import libraries

import random as rnd
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
from trafficsim import*


# -------------------- Calculation

# iteration for positions vs time, nCars = 10

roadLength = 50
vMaxArray = [1,2,5]
pBrake = 0.5
tStop = 100      # stop calculating at this time
nCars = 10
tArray = [t for t in range(1,tStop+1)]

carPositionsDictArray = []  # storing lists in dicts in a list

for vMax in vMaxArray:

    carPositions = [x for x in range(nCars)]
    carVelocities = [0 for x in range(nCars)]
    t = 0
    carPositionsDict = {}       # keys are t in tArray

    while t < tStop:
        t += 1
        traffic(nCars, carVelocities, vMax,
                carPositions, roadLength, pBrake)

        # store data
        carPositionsDict[t] = carPositions[:]
    carPositionsDictArray.append(carPositionsDict)


# calculate flow rate for fundamental diagram

# variables

roadLength = 50
vMaxArray = [1,2,5]
densityArray = [0, 0.04, 0.1, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] 
pBrake = 0.5
tStop = 100      # stop calculating at this time
nCars = 10
tArray = [t for t in range(1,tStop+1)]
nSimulations = 1000      # number of simulations for each nCars

## storage

flowrateDict = {}       # keys are items in roadLengthArray

for vMax in vMaxArray:
    print(vMax)
    flowrateArray = []      # will contain (nCars number of) arrays (of length nSimulations) corresponding to densityArray

    for density in densityArray:
        print(density)
        nCars = int(density * roadLength)
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
        flowrateArray.append(flowrateArray_nCars[:])
    flowrateDict[vMax] = flowrateArray[:]


# -------------------- Plot

def plotTraffic():
    fig = plt.figure()

    colors = pl.cm.jet(np.linspace(0.1,0.9,20))

    for i in range(len(vMaxArray)):
        vMax = vMaxArray[i]
        carPositionsDict = carPositionsDictArray[i]
        plt.clf()
        plt.title('Traffic simulation for $v_{max} = %d$'%(vMax),fontsize = 14)
        plt.xlabel('Time', fontsize = 14)
        plt.ylabel('Position', fontsize = 14)
        
        for j in range(len(tArray)):
            t = tArray[j]
            carPositions = carPositionsDict[t]
            for k in range(len(carPositions)):
                plt.plot(t, carPositions[k],
                    marker = '.',
                    markersize = 2,
                    color = colors[k]
                )
        plt.show()

def plotDiagram(i):
    vMax = vMaxArray[i]
    flowrateArray = flowrateDict[vMax]

    # regression
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

    plt.plot(xArray, [fFitted(x, a,b,c) for x in xArray], 
             label = '$V_{max}$ = ' + str(vMax))
    plt.legend()
    plt.title('Fundamental diagrams for $v_{max} = 1,2,5$', fontsize = 14)
    plt.ylabel('Flow rate - $\Sigma V_{cars}/ L$', fontsize = 14)
    plt.xlabel('Density - $N_{cars} / L$', fontsize = 14)

def output():
    plotTraffic()
    plt.figure()
    plotDiagram(0)
    plotDiagram(1)
    plotDiagram(2)
    plt.show()


output()

