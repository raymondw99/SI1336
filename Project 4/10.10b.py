
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import curve_fit
from pylab import cm as plcm            # colors for plotting
import pickle                           # for storing calculations and reading again while editing plot details

plt.rcParams.update({'font.size': 14})

# ------------------------- Set global values

nToExamine = [5,10]
nArray = [n for n in range(2,30)]
tolerance = 0.01

stepsToToleranceArray1 = [1, 4, 8, 13, 20, 26, 35, 44,
                          56, 67, 81, 94, 110, 126, 144,
                          162, 182, 203, 225, 248, 272,
                          297, 324, 352, 381, 410, 442, 473]


# ------------------------- Define functions

def initiateVMatrixes():
    """Initiates potential matrixes.
    - v has boundary values and initial guess of 4 in the center
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
    if n % 2 == 0:
        v[n//2,n//2] = 4
    else:
        v[n//2, n//2] = 4
        v[n//2+1, n//2] = 4
        v[n//2, n//2+1] = 4
        v[n//2+1, n//2+1] = 4
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

calculate()         # comment to obly load data
commonLinewidth = 2

plt.figure() # linear plot
plt.plot(nArray, stepsToToleranceArray1,
        linewidth = commonLinewidth,
        marker = '.',
        markersize = 8,
        color = 'tab:blue')
plt.plot(nArray, stepsToToleranceArray,
        linewidth = commonLinewidth,
        marker = '.',
        markersize = 8,
        color = 'b')
plt.legend(['Good initial guess', 'Poor initial guess'])
plt.title('Relaxations for tolerance ($1 \%$ accuracy)', fontsize = 18)
plt.ylabel('Relaxations', fontsize = 18)
plt.xlabel('Grid size (n)', fontsize = 18 )
plt.show()

plt.figure() # loglog plot
# calculate regression
log_nArray = [np.log(n) for n in nArray]
log_stepsToToleranceArray = [np.log(steps) for steps in stepsToToleranceArray]
def f(x,k,m):
    return k*x+m
[k,m], _ = curve_fit(f, log_nArray[10:], log_stepsToToleranceArray[10:])
# plot loglog plot
plt.loglog(nArray, stepsToToleranceArray,
        linewidth = commonLinewidth,
        color = 'tab:red',
        marker = 'x',
        markersize = 6)
plt.loglog(nArray, [np.e**f(x,k,m) for x in log_nArray],
        linewidth = commonLinewidth,
        linestyle = '--',
        zorder = 3,
        color = 'r')
plt.legend(['Computational data',r'log(f(n))$ = %s \, \log(n) - %s$'%(round(k,3), round(np.abs(m),3))])
plt.title('Relaxations for tolerance (log-log)', fontsize = 18)
plt.ylabel('Relaxations', fontsize = 18)
plt.xlabel('Grid size (n)', fontsize = 18)
plt.show()
