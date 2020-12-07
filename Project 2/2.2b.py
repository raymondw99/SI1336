
# -------------------- Import libraries

import random as rnd
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
from trafficsim import*
# -------------------- Calculation

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

nCars = 25
roadLength = 50
vMax = 2
pBrake = 0.5
tStop = 100

# calculate flowrate at this time, after initial transition

def standardError(N):
    flowrateArray = []      # will contain (nCars number of) arrays (of length N) corresponding to densityArray
    densityArray =  []
    flowrateArray_nCars = []
    for n in range(N):
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
    error = np.std(flowrateArray)/np.sqrt(N)
    return error

def output():
    standardErrorArray = []
    simulationArray = [2100,2200,2300,2400,2500,2600,2700,2800,2900,3000]
    for i in simulationArray:
        standardErrorArray.append(round(standardError(i),5))
        print((round(standardError(i),5)))
        
    plt.ticklabel_format(style='plain', axis='x', useOffset=False)
    plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(100))
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.5f'))
    plt.plot(simulationArray,standardErrorArray)
    plt.title("Standard error for flow rate")
    plt.xlabel("N - Number of simulations", fontsize = 12)
    plt.ylabel("SE - Standard error", fontsize = 12)
    #plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.show()

output()

