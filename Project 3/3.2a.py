# Performs a molecular dynamics simulation of particles in 2 dimensions

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
import random as rnd
from matplotlib import ticker
from data import*

matplotlib.rcParams.update({'font.size': 20})
#matplotlib.use('TkAgg')

# Generate two Gaussian random numbers with standard deviation sigma, mean 0
def gaussianRandomNumbers(sigma):
    w = 2
    while (w >= 1):
        rx1 = 2 * rnd.random() - 1
        rx2 = 2 * rnd.random() - 1
        w = rx1 * rx1 + rx2 * rx2
    
    w = math.sqrt(-2 * math.log(w) / w);
    return sigma * rx1 * w, sigma * rx2 * w


# Assigns Gaussian distributed velocities given an energy per particle
def thermalize(vx, vy, kineticEnergyPerParticle):
    for i in range(0, len(vx)):
        vx[i], vy[i] = gaussianRandomNumbers(kineticEnergyPerParticle)

# The pair potential
def pairEnergy(r):
    # Implement the pair energy calculation here
    epsilon = 1
    sigma = 1
    V = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return V


# The pair force              
def pairForce(r):
    # Implement the pair force calculation here
    sigma = 1
    epsilon = 1
    F = -4*epsilon*((-12/r)*((sigma/r)**12) + (6/r)*((sigma/r)**6))
    return F

# Calculate the shortest periodic distance, unit cell [0,Lx],[0,Ly]
# Returns the difference along x, along y and the distance
# This code assumes all particles are within [0,Lx],[0,Ly]
def pbc_dist(x1, y1, x2, y2, Lx, Ly):
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
    return dx, dy, math.sqrt(dx*dx + dy*dy)


# Number of particles
n = 24
# Mass of a particle
mass    = 1.0
invmass = 1.0/mass

# The number of particles per row in the initial grid
numPerRow = 6

# Unit cell size
Lx = numPerRow*1.12
Ly = numPerRow*1.12

# Boltzmann constant
kB = 1.0
# Temperature
T = 1

# kb*T
kBT = kB * T

dt = 0.1
nsteps = 20000
# Use a value of 1 if you want to see in detail how the particles interact
numStepsPerFrame = 100
numFrames = nsteps//numStepsPerFrame

#rnd.seed(1)

x = []
y = []
vx = []
vy = []
fx = []
fy = []

# Initialize the particle position to a nearly hexagonal lattice
x = []
y = []
for i in range (0,n):
    x.append(Lx*0.95/numPerRow*((i % numPerRow) + 0.5*(i/numPerRow)))
    y.append(Lx*0.95/numPerRow*0.87*(i/numPerRow))

    # Make the velocity lists
    vx.append(0)
    vy.append(0)

    # Make the force lists
    fx.append(0)
    fy.append(0)
    
# Assign velocities from a Gaussian distribution
thermalize(vx, vy, kBT)

fig = plt.figure()
ax  = plt.subplot(xlim=(0, Lx), ylim=(0, Ly))


# Start recording for Ekin and the heat capacity at time 5
startTimeForAveraging = 100
startStepForAveraging = startTimeForAveraging/0.01

step = 0

# Initialize observables
sumEpot  = 0
sumEpot2 = 0

outt = []
ekinList = []
epotList = []
etotList = []


# Perform one MD integration step
def integrate():
    global step, n, x, y, v, Ekin, Epot, sumEpot, sumEpot2
    
    # Clear the energy and potential
    Epot = 0
    for i in range(0, n):
        fx[i] = 0
        fy[i] = 0

    # Compute the pair potentials and forces
    for i in range(0,n):
        for j in range(i+1,n):
            dx, dy, r = pbc_dist(x[i],y[i],x[j],y[j],Lx,Ly)
            Epot  += pairEnergy(r)
            fij = pairForce(r)
            fx[i] += fij * dx / r
            fy[i] += fij * dy / r
            fx[j] -= fij * dx / r
            fy[j] -= fij * dy / r

    if step > startStepForAveraging:
        sumEpot  += Epot
        sumEpot2 += Epot * Epot

    Ekin = 0
    for i in range(0,n):

        # At the first step we alread have the "full step" velocity
        if step > 0:
            # Update the velocities with a half step
            vx[i] += fx[i]*invmass*0.5*dt
            vy[i] += fy[i]*invmass*0.5*dt

        # When temperature coupling, modify the velocity of one particle here


        # Add the kinetic energy of particle i to the total
        Ekin += 0.5*mass*(vx[i]*vx[i] + vy[i]*vy[i])

        # Update the velocities with a half step
        vx[i] += fx[i] *invmass*0.5*dt
        vy[i] += fy[i]*invmass*0.5*dt

        # Update the coordinates
        x[i]  += vx[i] * dt
        y[i]  += vy[i] * dt

        # Put particles back in the unit cell, useful for visualization
        x[i] = x[i] % Lx
        y[i] = y[i] % Ly

    # Increase the time step
    step += 1

# Perform several MD integration steps and record observables
def integrate_some_steps(framenr):
    global numStepsPerFrame, start, step, kBT, x, y, v, Ekin, Epot, sv, svv

    for i in range(0, numStepsPerFrame):
        integrate()
    
    t = step*0.01
    outt.append(t)
    ekinList.append(Ekin)
    epotList.append(Epot)
    etotList.append(Epot + Ekin)

    if step >= startStepForAveraging and step % 1000 == 0:
        EpotAv  = sumEpot/(step + 1 - startStepForAveraging)
        Epot2Av = sumEpot2/(step + 1 - startStepForAveraging)
        Cv = (Epot2Av - EpotAv * EpotAv) / (kBT * T)
        print('time', t, 'Cv =', Cv)

    return ax.scatter(x, y, s=1500, marker='o', c="r"),


# Call the animator, blit=True means only re-draw parts that have changed
##anim = animation.FuncAnimation(fig, integrate_some_steps,
##                               frames=numFrames, interval=50, blit=True, repeat=False)

# Depending on how you run python you might want to remove this plt.show()
# When running from Anaconda Spyder the show() needed to avoid a long delay.
# When running from the command line this is often not needed and leads to
# having to close the animation window and click on an empty window that appears

plt.show()  # show the animation
#plt.waitforbuttonpress(timeout=20)

##plt.figure()
##fig, ax = plt.subplots(1,1)
##ax.plot(outt, etotList)
##ax.tick_params(axis='both', which='major', labelsize=12)
##plt.title('Energy drift', fontsize = 12)
##plt.xlabel('Time', fontsize = 12)
##plt.ylabel('Total energy', fontsize = 12)
##ax.yaxis.get_offset_text().set_fontsize(12)
##print(etotList)

def mean(energies):
    mean = sum(energies)/len(energies)
    return mean

fig, ax = plt.subplots(1,1)
ax.plot(timelist, [x-mean(etotlist000001) for x in etotlist000001],
                   'tab:orange',label = '$dt = 10^{-5}$')
ax.plot(timelist, [x-mean(etotlist00001) for x in etotlist00001],
                   'g',label = '$dt = 10^{-4}$')
ax.plot(timelist, [x-mean(etotlist0001) for x in etotlist0001],
                   'r',label = '$dt = 10^{-3}$')
ax.plot(timelist, [x-mean(etotlist001) for x in etotlist001],
                   label = '$dt = 10^{-2}$')
ax.tick_params(axis='both', which='major', labelsize=14)
plt.title('Deviation from total energy', fontsize = 14)
plt.xlabel('Time', fontsize = 14)
plt.ylabel('Deviation', fontsize = 14)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax.yaxis.get_offset_text().set_fontsize(14)
plt.legend(loc='upper right', prop={'size': 14})
#plt.legend( ('Ekin','Epot','Etot') )

plt.show()
