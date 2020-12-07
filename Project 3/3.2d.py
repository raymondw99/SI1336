import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
import random as rnd
from data import*

matplotlib.rcParams.update({'font.size': 20})
#matplotlib.use('TkAgg')

# The potential
def V(r):
    epsilon = 1
    sigma = 1
    V = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return V

# Calculate the shortest periodic distance, unit cell [0,Lx],[0,Ly]
# This code assumes all particles are within [0,Lx],[0,Ly]
def pbc_dist(x1,y1,x2,y2,Lx,Ly):
    dx = x1 - x2
    dy = y1 - y2
    while dx < -0.5*Lx:
        dx += Lx
    while dx > 0.5*Lx:
        dx -= Lx
    while dy < -0.5*Ly:
        dy += Ly
    while dy > 0.5*Ly:
        dy -= Ly
    return math.sqrt(dx*dx + dy*dy)

# Number of particles
n = 20

# Unit cell size
Lx = 5.6
Ly = 5.6

# Boltzmann constant
kB = 1
# Temperature
T = 1.0

# kb*T
kBT = kB * T

start    = 1000
numsteps = 50000

numstepsperframe = 100
numframes        = int(numsteps/numstepsperframe)
rnd.seed()

outt = []
epot = []
v_avlist = []
# The maximum step size
# FOR OPTIMAL PERFORMANCE, TUNE THIS FOR THE GIVEN SYSTEM AND DENSITY
delta = 1

# Initialize the particle position to a nearly hexagonal lattice
x = []
y = []
for i in range (0,n):
    x.append(Lx/5*((i % 5) + 0.5*(int(i/5))))
    y.append(Lx/5*0.87*(int(i/5)))

# Determine the initial value of the potential
v = 0
for i in range(0,n):
    for j in range(i+1,n):
        r = pbc_dist(x[i],y[i],x[j],y[j],Lx,Ly)
        v += V(r)

fig = plt.figure()
ax  = plt.subplot(xlim=(0, Lx), ylim=(0, Ly))

step = 0

# The sum of the acceptance ratio over start to numsteps
sar = 0

# The sum of the potential and the potential squared
sv  = 0
svv = 0

# Perform one MC sweep over all particles
def mc_sweep():
    global n, start, step, x, y, v, sar, sv, svv
    
    for i in range(0,n):
        xn = x[i] + delta*(2*rnd.random() - 1)
        yn = y[i] + delta*(2*rnd.random() - 1)
        dV = 0
        for j in range(0,n):
            if i != j:
                r = pbc_dist(x[i],y[i],x[j],y[j],Lx,Ly)
                rn = pbc_dist(xn, yn, x[j], y[j], Lx, Ly)
                dV  += V(rn) - V(r)

        ratio = math.exp(-dV/kBT)

        sar += min(1,ratio)

        if ratio > rnd.random():
            # Accept the new configuration
            x[i] = xn
            y[i] = yn
            v += dV

        # Put the particle back in the box
        if x[i] < 0:
            x[i] += Lx
        if x[i] >= Lx:
            x[i] -= Lx
        if y[i] < 0:
            y[i] += Ly
        if y[i] >= Ly:
            y[i] -= Ly

    if step >= start:
        sv  += v
        svv += v*v

    if step % numstepsperframe == 0:
        outt.append(step)
        epot.append(v)
        if step > start:
            v_av  = sv/(step + 1 - start)
            vv_av = svv/(step + 1 - start)
            v_avlist.append(v_av)
            cVlist.append((vv_av - v_av*v_av)/(kB*T*T))
            print("<V>", v_av, "cV", (vv_av - v_av*v_av)/(kB*T*T), "ratio", sar/((step+1)*n))

    step += 1



def animate(framenr):
    global numstepsperframe, n, step, x, y, v, sar, sv, svv
    for it in range(numstepsperframe):
        mc_sweep()


    
    return ax.scatter(x, y, s=1500, marker='o', c="r"),

# Call the animator, blit=True means only re-draw parts that have changed
#anim = animation.FuncAnimation(fig, animate,
#                               frames=numframes, interval=50, blit=True, repeat=False)

# Depending on how you run python you might want to remove this plt.show()
# When running from Anaconda Spyder the show() needed to avoid a long delay.
# When running from the command line this is often not needed and leads to
# having to close the animation window and click on an empty window that appears

#plt.show()  # show the animation
#plt.waitforbuttonpress(timeout=30)

##plt.figure
##plt.plot(outt, epot)
##plt.draw()


fig, ax = plt.subplots(1,1)
plt.xlabel('Temperature', fontsize = 14)
plt.ylabel('V (J) - $C_{V}$ (J/K)', fontsize = 14)
plt.title('Average energy (V) and heat capacity ($C_{V}$)', fontsize = 14)
ax.plot(templist, cVlist, 'tab:red', label = 'Heat capacity')
ax.plot(templist, cVlist, 'ro')
ax.plot(templist, Vlist, 'tab:orange', label = 'Average energy')
ax.plot(templist, Vlist,'o', color="tab:orange")
ax.tick_params(axis='both', which='major', labelsize=14)
plt.legend(loc='upper right', prop={'size': 14})
plt.show()

