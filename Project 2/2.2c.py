
# -------------------- Import libraries

import random as rnd
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit
from matplotlib import animation
from trafficsim import*

## variables

densityArray =  [0, 0.04, 0.1, 0.14, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] 
nSimulations = 1000      # number of simulations for each nCars
#roadLength = 50
vMax = 2
pBrake = 0.5
tStop = 100             # calculate flowrate at this time, after initial transition

def meanFlowrateArray(roadLength):
    flowrateArray = []
    nCarsArray = []
    for n in densityArray:
        nCarsArray.append(int(n*roadLength))

    ## storage

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
        flowrateArray.append(flowrateArray_nCars[:])
        meanFlowrateArray = [np.mean(array) for array in flowrateArray]
    print(meanFlowrateArray)
    return meanFlowrateArray

# -------------------- Plot

def plot(roadLength, color):
    # regression

    def fFitted(x,a,b,c):
        return (a * (x**b)) * np.exp(-c*x)

    popt, pcov = curve_fit(fFitted, densityArray, meanFlowrateArray(roadLength))
    a = popt[0]
    b = popt[1]
    c = popt[2]
    aRounded = round(a,2)
    bRounded = round(b,2)
    cRounded = round(c,2)
    
    xArray = np.linspace(min(densityArray), max(densityArray), 1001)
    plt.plot(xArray, [fFitted(x, a,b,c) for x in xArray],
             color,label = 'L = ' + str(roadLength))
    plt.legend()
    plt.title('Length dependence for fundamental diagrams')
    plt.ylabel('Flow rate - $\Sigma V_{cars}/ L$', fontsize = 12)
    plt.xlabel('Density - $N_{cars} / L$', fontsize = 12)

def output():
    plt.figure()
    plot(5,'tab:blue')
    plot(10,'tab:orange')
    plot(25,'tab:green')
    plot(50,'tab:red')
    plot(100,'tab:purple')
    plot(150,'tab:brown')
    plt.show()

output()
