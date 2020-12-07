import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import curve_fit
from pylab import cm as plcm            # colors for plotting
import pickle                           # for storing calculations and reading again while editing plot details

plt.rcParams.update({'font.size': 14})

# ------------------------- Set global values

nToExamine = [5,10]                 # assignmet
nArray = [n for n in range(2,30)]   # range to find dependence on n
tolerance = 0.01                    # max relative error from exact solution

# ------------------------- Define functions

def initiateVMatrixes():
    """Initiates potential matrixes.
    - v has boundary values and initial guess of 9 everywhere else
    - vNew is a copy of v
    - vExact is the exact analytical solution, 10 everywhere"""
    global v, vNew, vExact
    # Initialize the grid to 0
    v = np.zeros((n+1, n+1))        # matrix of v, index are i: row, j:column
    # Set the boundary conditions
    for i in range(1,n):
        v[0,i] = 10
        v[n,i] = 10
        v[i,0] = 10
        v[i,n] = 10
    # Exact solution
    vExact = np.copy(v)
    for i in range(1,n):
        for j in range(1,n):
            vExact[i,j] = 10
    # Initial guess
    for i in range(1,n):
        for j in range(1,n):
            v[i,j] = 0.9*vExact[i,j]
    vNew = np.copy(v)

def relax():
    """One relax iteration. v[i,j] is set as the avarage of its neighbours."""
    global v, vNew, n
    for x in range(1,n):
        for y in range(1,n):
            vNew[x,y] = (v[x-1][y] + v[x+1][y] + v[x][y-1] + v[x][y+1])*0.25
    for x in range(1,n):
        for y in range(1,n):
            v[x,y] = vNew[x,y]

def calculate():
    """Main calculation function that first initalizes with initiateVMatrixes()
    and then uses relax() until v is within tolerance.
    1. Iterate for n = 5, 10
    2. Iterate for the range of n in nArray"""
    global v, vNew, n, stepsToToleranceArray
    # n = 5, 10
    for n in nToExamine:
        initiateVMatrixes()
        step = 0
        toleranceAcqurired = False
        while not toleranceAcqurired:
            step+=1
            relax()
            # Controll accuracy
            toleranceAcqurired = True   # run through v and set false if not acquired
            for i in range(1,n):
                for j in range(1,n):
                    if np.abs( (v[i,j]-vExact[i,j])/vExact[i,j] ) > tolerance:
                        toleranceAcqurired = False
            if toleranceAcqurired:
                print('Tolerance for n =', n, 'was met after', step, 'steps.')

    # range of n
    stepsToToleranceArray = []
    for n in nArray:
        print('Currently working with n = ', n)
        initiateVMatrixes()
        step = 0
        toleranceAcqurired = False
        while not toleranceAcqurired:
            step+=1
            relax()
            # Controll accuracy
            toleranceAcqurired = True   # run through v and set false if not acquired
            for i in range(1,n):
                for j in range(1,n):
                    if np.abs( (v[i,j]-vExact[i,j])/vExact[i,j] ) > tolerance:
                        toleranceAcqurired = False
            if toleranceAcqurired:
                stepsToToleranceArray.append(step)
                
# ----------------------- Plot

calculate()
# common plot values
commonLinewidth = 2

#plt.subplot(2,1,1) # linear plot
plt.figure()
print(stepsToToleranceArray)
print(nArray)
pol = np.poly1d(np.polyfit(nArray, stepsToToleranceArray, 2))
plt.plot(nArray, pol(nArray),
        linewidth = commonLinewidth,
        marker = '.',
        markersize = 10,
        color = 'tab:red', label = 'Computational data')
plt.plot(nArray, pol(nArray),
        linewidth = commonLinewidth,
         linestyle = '--',
        color = 'r', label = '$f(n) = 0.56n^2$')
plt.legend()
plt.title('Relaxations for tolerance ($1 \%$ accuracy)', fontsize = 18)
plt.ylabel('Relaxations', fontsize = 18)
plt.xlabel('Grid size (n)', fontsize = 18 )
plt.show()

plt.figure()
# calculate regression from n = 10 and up
log_nArray = [np.log(n) for n in nArray]
log_stepsToToleranceArray = [np.log(steps) for steps in stepsToToleranceArray]
def f(x,k,m):
    return k*x+m
[k,m], _ = curve_fit(f, log_nArray[10:], log_stepsToToleranceArray[10:])
# plot loglog with regression
plt.title('Relaxations for tolerance (log-log)', fontsize = 18)
plt.loglog(nArray, stepsToToleranceArray,
        linewidth = commonLinewidth,
        color = 'tab:blue',
        marker = '',
        markersize = 10, label = 'Computational data')
plt.loglog(nArray, [np.e**f(x,k,m) for x in log_nArray],
        linewidth = commonLinewidth,
        linestyle = '--',
        zorder = 3,
        color = 'b', label = '$log f(n) = 2.007log(n) - 0.597$')
plt.legend()
plt.ylabel('Relaxations', fontsize = 18)
plt.xlabel('Grid size (n)', fontsize = 18)
plt.show()

