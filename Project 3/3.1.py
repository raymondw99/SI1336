# Computes the volume of a 10-dimensional sphere using midpoint integration

import math as math
import numpy as np
from scipy.optimize import curve_fit
from time import process_time
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator

# Returns a list of n equidistant points between -bound and bound
def discretization(bound, n):
	return np.linspace(-float(bound) + float(bound) / n, float(bound) + float(bound) / n, n, False)

# Recursively computes an integral of a dim-dimensonal sphere
# with radius sqrt(radius). Pass volume2 = 1 at start.
# n is the number of points used for the midpoint method.
def recursiveIntegral(radius2, volume2, dim, n):
    volume2 *= radius2 * 2 * 2  / (n * n)
    if (dim > 1):
        partIntegral = 0
        for x in discretization(math.sqrt(radius2), n):
            partIntegral += recursiveIntegral(radius2 - x * x, volume2, dim - 1, n)
    else:
        partIntegral = math.sqrt(volume2) * n

    return partIntegral


# Monte Carlo approximation of the volume of a dim-dimensonal sphere
# with radius sqrt(radius)
def montecarlo(dim, radius, N):
    count_in_sphere = 0

    for count_loops in range(N):
        point = np.random.uniform(-1.0, 1.0, dim)
        distance = np.linalg.norm(point)
        if distance < 1.0:
            count_in_sphere += 1

    return np.power(2.0, dim) * (count_in_sphere / N)

# The number of dimensions
numDims = 10

timelist = []
errorlist = []

errorlist1 = [0.9682663823341708, 0.1660809555237166, 0.15717487652157347,
             0.1251325819655884, 0.09438377397356268] 
timelist1 = [0.10271, 0.429084, 3.925878, 23.295996, 111.834063] 

errorlist2 = [1.1533329367039888, 0.5112814632959886, 0.31594107517154857,
              0.10376734000000003, 0.10376734000000003, 0.057765464637436814,
              0.03379026373290017]
timelist2 = [0.02383748, 0.04916397, 0.07693159, 0.10376734000000003,
             0.2895373, 0.5990030300000001, 0.9335681499999999]
def calculate(t, error):
    print('time =', t)
    print('error   =', error)
    print('')
    timelist.append(t)
    errorlist.append(error)

def graph(errorlist, timelist, col, title):
    fig, ax = plt.subplots()
    plt.title(title, fontsize = 18)
    ax.plot(errorlist, timelist, color = col)
    ax.plot(errorlist, timelist, 'o', color = col)
    plt.xlabel("Integration error", fontsize = 18)
    plt.ylabel("Computational time (s)", fontsize = 18)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.show()

# The number of points in the midpoint method along one dimension

pointlist = [2, 3, 4, 5]
analytical = math.pi**(numDims / 2) / math.factorial(numDims / 2)

iterationlist = [750,1000,2500,5000,7500,
                10000, 25000, 50000, 75000, 100000]

##for point in pointlist:
##    t = process_time()
##    integral = recursiveIntegral(1, 1, numDims, point)
##    t = process_time() - t
##    error = abs(integral- analytical)
##    calculate(t, error)

for iteration in iterationlist:
    f2 = []
    t = process_time()
    for i in range(100):
        f2.append(montecarlo(numDims, 1, iteration))
    t = (process_time() - t)/100
    f1 = [i**2 for i in f2]
    error = sqrt((sum(f1)/len(f1) - (sum(f2)/len(f2))**2)/100)
    calculate(t, error)

graph(errorlist, timelist, 'tab:red', 'Monte Carlo estimation')                           
#graph(errorlist1, timelist1, 'tab:blue', 'Midpoint method')  
#graph(errorlist1, timelist1, 'tab:blue','Midpoint method')

