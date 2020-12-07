import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 14})

# ------------------------- Set global values

n = 10
tolerance = 0.000001

# ------------------------- Define functions

def initiateVMatrixes_1():
    """Initiates potential matrixes with first boundary value: V = 10 on
    two opposing sides and 5 on the other two sides.
    - v has boundary values and initial guess of 9 everywhere else
    - vNew is a copy of v
    - vExact is the exact analytical solution, 10 everywhere"""
    global v, vNew
    # Initialize the grid to 0
    v = np.zeros((n+1, n+1))        # matrix of v, index are i: row, j:column
    # Set the boundary conditions
    for i in range(1,n):
        v[0,i] = 10
        v[n,i] = 10
        v[i,0] = 5
        v[i,n] = 5
    # Initial guess
    for i in range(1,n):
        for j in range(1,n):
            v[i,j] = 7.5
    vNew = np.copy(v)



def initiateVMatrixes_2():
    """Initiates potential matrixes with second boundary conditions: V = 10 on
    three sides and 0 on the fourth.
    - v has boundary values and initial guess of 9 everywhere else
    - vNew is a copy of v
    - vExact is the exact analytical solution, 10 everywhere"""
    global v, vNew
    # Initialize the grid to 0
    v = np.zeros((n+1, n+1))        # matrix of v, index are i: row, j:column
    # Set the boundary conditions
    for i in range(1,n):
        v[0,i] = 10
        v[n,i] = 10
        v[i,0] = 10
        v[i,n] = 0
    v[0,0]=10; v[n,0]=10
    # Initial guess
    for i in range(1,n):
        for j in range(1,n):
            v[i,j] = 10-(j/(n-1))*10
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

def calculate1():
    """Main calculation function that first initalizes with initiateVMatrixes()
    and then uses relax() until v is within tolerance.
    1. First boundary coditions
    2. Second boundary conditions"""
    global v, vNew, n, v1

    # First bondary conditions
    initiateVMatrixes_1()
    step = 0
    toleranceAcqurired = False
    while not toleranceAcqurired:
        if step%100==0:print('v =',v,'\nstep =',step)
        step+=1
        vOld = np.copy(v)
        relax()
        # Controll accuracy
        toleranceAcqurired = True   # run through v and set false if not acquired
        for i in range(1,n):
            for j in range(1,n):
                if np.abs( (v[i,j]-vOld[i,j])/vOld[i,j] ) > tolerance:
                    toleranceAcqurired = False
    print('Tolerance for n =', n, 'was met after', step, 'steps with first boundary conditions.')
    v1 = np.copy(v)
    return v1

def calculate2():
    global v, vNew, n, v2
    # Second boundary conditions
    initiateVMatrixes_2()
    step = 0
    toleranceAcqurired = False
    while not toleranceAcqurired:
        if step%100==0:print('v =',v,'\nstep =',step)
        step+=1
        vOld = np.copy(v)
        relax()
        # Controll accuracy
        toleranceAcqurired = True   # run through v and set false if not acquired
        for i in range(1,n):
            for j in range(1,n):
                if np.abs( (v[i,j]-vOld[i,j])/vOld[i,j] ) > tolerance:
                    toleranceAcqurired = False
    print('Tolerance for n =', n, 'was met after', step, 'steps with second boundary conditions.')
    v2 = np.copy(v)
    return v2   

# ----------------------- Plot

commonColormap = 'viridis'         # alternatives: inferno, plasma, Greys, Blues, BuPu, bone, afmhot
commonInterpolation = 'bicubic' # alternatives: nearest

v1 = calculate1()
v2 = calculate2()

plt.figure()
# 1 plot potential with colormap
plt.title('Initial boundary condition', fontsize = 18)
im1 = plt.imshow(v1,            # plot first result
        cmap=commonColormap,
        interpolation=commonInterpolation)
plt.colorbar(im1, orientation='vertical',    # add colorbar below plots
                  fraction=.05)
plt.xlabel('x', fontsize = 18)
plt.ylabel('y', fontsize = 18)
plt.gca().invert_yaxis()
plt.show()

plt.figure()
# 2 plot potential with colormap
plt.title('Second boundary condition',fontsize = 18)
im1 = plt.imshow(v2,            # plot first result
        cmap=commonColormap,
        interpolation=commonInterpolation)
plt.colorbar(im1, orientation='vertical',    # add colorbar below plots
                  fraction=.05)
plt.xlabel('x', fontsize = 18)
plt.ylabel('y', fontsize = 18)
plt.gca().invert_yaxis()
plt.show()


# 1 plot potential with countor
plt.figure()
cs1 = plt.contour(v1,levels=np.arange(0, 15, 1.3), colors='black', linestyles='dashed')
plt.clabel(cs1, inline=1, fontsize=10)
plt.title('Initial boundary condition')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# 2 plot potential with countor
plt.figure()
cs1 = plt.contour(v2,levels=np.arange(0, 15, 1.3), colors='black', linestyles='dashed')
plt.clabel(cs1, inline=1, fontsize=10)
plt.title('Second boundary condition')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

