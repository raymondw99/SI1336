import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import curve_fit
from pylab import cm as plcm            # colors for plotting
import pickle                           # for storing calculations and reading again while editing plot details

plt.rcParams.update({'font.size': 14})

# ------------------------- Set global values

nToExamine = [5,10]                 # assignment
nArray = [n for n in range(2,30)]   # range to find dependence on n
tolerance = 0.01                    # max relatice error from exact solution


stepsToToleranceArray1 = [1, 3, 5, 8, 10, 14, 18, 23, 29,
                          35, 41, 48, 56, 64, 73, 82, 92,
                          103, 113, 125, 137, 150, 163,
                          177, 191, 206, 222, 238]
stepsToToleranceArray2 = [1, 4, 8, 13, 20, 26, 35, 44,
                          56, 67, 81, 94, 110, 126, 144,
                          162, 182, 203, 225, 248, 272,
                          297, 324, 352, 381, 410, 442, 473]

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

def relax_checker():
    """One checker-relax iteration. v[i,j] is set as the avarage of its neighbours."""
    checker = 2
    global v, vNew, n
    for check in range(0,2):
        for x in range(1,n):
            for y in range(1,n):
                if (x*(n+1) + y) % 2 == check:
                    v[x,y] = (v[x-1][y] + v[x+1][y] + v[x][y-1] + v[x][y+1])*0.25

def calculate():
    """Main calculation function that first initalizes with initiateVMatrixes()
    and then uses relax() until v is within tolerance.
    1. Iterate for n = 5, 10
    2. Iterate for the range of n in nArray"""
    global v, vNew, n, stepsToToleranceArray
    stepsToToleranceArray = []
    for n in nArray:
        print('Currently working with n = ', n)
        initiateVMatrixes()
        step = 0
        toleranceAcqurired = False
        while not toleranceAcqurired:
            step+=1
            relax_checker()
            # Controll accuracy
            toleranceAcqurired = True   # run through v and set false if not acquired
            for i in range(1,n):
                for j in range(1,n):
                    if np.abs( (v[i,j]-vExact[i,j])/vExact[i,j] ) > tolerance:
                        toleranceAcqurired = False
            if toleranceAcqurired:
                stepsToToleranceArray.append(step)
        if n in [5,10]: print('n =', n, 'steps =', step)


# ----------------------- Main

calculate()         # comment to obly load data

# ----------------------- Plot

# common plot values
commonLinewidth = 2

plt.figure() # linear plot
plt.plot(nArray, stepsToToleranceArray,
        linewidth = commonLinewidth,
        linestyle = '--',
        marker = '.',
        markersize = 8,
        color = 'tab:cyan',
        zorder = 2)
plt.plot(nArray, stepsToToleranceArray1,
        linewidth = commonLinewidth,
        linestyle = '-',
        marker = 'x',
        markersize = 8,
        color = 'tab:blue',
        zorder = 1)
plt.plot(nArray, stepsToToleranceArray2,
        linewidth = commonLinewidth,
        marker = '.',
        markersize = 8,
        color = 'b',
        zorder = 0)
plt.title(r'Steps required to reach tolerance', fontsize = 18)
plt.legend(['Checker relaxation', 'Gauss-Seidel', 'Relaxation method'])
plt.ylabel('Relaxations', fontsize = 18)
plt.xlabel('Grid size (n)', fontsize = 18 )
plt.show()

plt.figure() # loglog plot
# calculate regression from n = 10 and up
log_nArray = [np.log(n) for n in nArray]
log_stepsToToleranceArray = [np.log(steps) for steps in stepsToToleranceArray]
def f(x,k,m):
    return k*x+m
[k,m], _ = curve_fit(f, log_nArray[10:], log_stepsToToleranceArray[10:])
# plot loglog with regression
plt.loglog(nArray, [np.e**f(x,k,m) for x in log_nArray],
        linewidth = commonLinewidth,
        linestyle = '--',
        zorder = 3,
        color = 'r')
plt.loglog(nArray, stepsToToleranceArray,
        linewidth = commonLinewidth,
        color = 'tab:red',
        marker = '.',
        markersize = 8)
plt.legend(['Computarional data', r'log(f(n))$ = %s \, \log(n) - %s$'%(round(k,3), round(np.abs(m),3))],
        fancybox = True)
plt.title('Relaxations for tolerance (log-log)', fontsize = 18)
plt.ylabel('Relaxations', fontsize = 18)
plt.xlabel('Grid size (n)', fontsize = 18)
plt.show()

