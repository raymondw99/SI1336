# Python simulation of a double pendulum with real time animation
# BH 2020-10-27

from matplotlib import animation
from pylab import *

# Global constants
G = 9.8  # gravitational acceleration


# Kinetic energy
def Ekin(osc):
    return 1 / (2.0 * osc.m * osc.L * osc.L) * (
            osc.p1 * osc.p1 + 2.0 * osc.p2 * osc.p2 - 2.0 * osc.p1 * osc.p2 * cos(osc.q1 - osc.q2)) / (
                   1 + (sin(osc.q1 - osc.q2)) ** 2)


# Potential energy
def Epot(osc):
    return osc.m * G * osc.L * (3 - 2 * math.cos(osc.q1) - math.cos(osc.q2))


# Class that holds the parameter and state of a double pendulum
class Oscillator:
    m = 1  # mass of the pendulum bob
    L = 1  # arm length

    t = 0  # the time

    def p2squaredFromH(self):
        return (self.E - Epot(self)) * (1 + (sin(self.q1 - self.q2)) ** 2) * self.m * self.L * self.L

    E = 15
    p2 = -1
    while (p2 < 0):
        #q1 = math.pi * (2 * np.random.random() - 1)
        #q2 = math.pi * (2 * np.random.random() - 1)
        q1 = 1.1
        q2 = 0
        p1 = 0
    
        ########################################
        # I cannot figure out how to pass the class instance itself to p2SquaredFromH
        # I should probably restructure this, but I don't know how

        p2squared = 1
        if (p2squared >= 0):
            p2 = math.sqrt(p2squared)


# Class for storing observables for an oscillator
class Observables:

    def __init__(self):
        self.time = []  # list to store time
        self.q1list = []  # list to store q1
        self.q2list = []  # list to store q2
        self.p1list = []  # list to store p1
        self.p2list = []  # list to store p2
        self.epot = []  # list to store potential energy
        self.ekin = []  # list to store kinetic energy
        self.etot = []  # list to store total energy
        self.poincare_q1 = []  # list to store q1 for Poincare plot
        self.poincare_p1 = []  # list to store p1 for Poincare plot


# Differentiate of H with respect to p1
def dHdp1(q1, q2, p1, p2, m, L):
    return 1 / (2.0 * m * L * L) * (2*p1 - 2.0 * p2 * cos(q1 - q2))/ (1 + (sin(q1 - q2)) ** 2)


# Differentiate of H with respect to p2
def dHdp2(q1, q2, p1, p2, m, L):
    return 1 / (2.0 * m * L * L) * (4.0 * p2  - 2.0 * p1 * cos(q1 - q2)) / (1 + (sin(q1 - q2)) ** 2)


# Differentiate of H with respect to q1
def dHdq1(q1, q2, p1, p2, m, L):
    return 1 / (2.0 * m * L * L) * (
            -2 * (p1 * p1 + 2 * p2 * p2) * cos(q1 - q2) + p1 * p2 * (4 + 2 * (cos(q1 - q2)) ** 2)) * sin(
        q1 - q2) / (1 + (sin(q1 - q2)) ** 2) ** 2 + m * G * L * 2.0 * sin(q1)


# Differentiate of H with respect to q2
def dHdq2(q1, q2, p1, p2, m, L):
    return 1 / (2.0 * m * L * L) * (
            2 * (p1 * p1 + 2 * p2 * p2) * cos(q1 - q2) - p1 * p2 * (4 + 2 * (cos(q1 - q2)) ** 2)) * sin(q1 - q2) / (
                   1 + (sin(q1 - q2)) ** 2) ** 2 + m * G * L * sin(q2)


class BaseIntegrator:
    dt = 0.003  # time step

    def integrate(self,
                  osc,
                  obs,
                  ):
        """ Perform a single integration step """
        self.timestep(osc, obs)

        """ Append observables to their lists """
        obs.time.append(osc.t)
        obs.q1list.append(osc.q1)
        obs.q2list.append(osc.q2)
        obs.p1list.append(osc.p1)
        obs.p2list.append(osc.p2)
        obs.epot.append(Epot(osc))
        obs.ekin.append(Ekin(osc))
        obs.etot.append(Epot(osc) + Ekin(osc))
        
        obs.poincare_q1.append(osc.q1)
        obs.poincare_p1.append(osc.p1)
        
    def timestep(self, osc, obs):
        """ Implemented by the child classes """
        pass

# Runge-Kutta 4 integrator
class RK4Integrator(BaseIntegrator):
    def timestep(self, osc, obs):
        dt = self.dt

        # TODO: Add integration here
        # Derivatives at time t
        dp1dt = -dHdq1(osc.q1, osc.q2, osc.p1, osc.p2, osc.m, osc.L) 
        dp2dt = -dHdq2(osc.q1, osc.q2, osc.p1, osc.p2, osc.m, osc.L)
        dq1dt = dHdp1(osc.q1, osc.q2, osc.p1, osc.p2, osc.m, osc.L) 
        dq2dt = dHdp2(osc.q1, osc.q2, osc.p1, osc.p2, osc.m, osc.L)

        # Constants in RK4
        k1p1 = dt * dp1dt 
        k1p2 = dt * dp2dt 
        k1q1 = dt * dq1dt 
        k1q2 = dt * dq2dt

        k2p1 = dt * (-dHdq1(osc.q1+k1q1*0.5, osc.q2+k1q2*0.5, osc.p1+k1p1*0.5, osc.p2+k1p2*0.5, osc.m, osc.L)) 
        k2p2 = dt * (-dHdq2(osc.q1+k1q1*0.5, osc.q2+k1q2*0.5, osc.p1+k1p1*0.5, osc.p2+k1p2*0.5, osc.m, osc.L))
        k2q1 = dt * dHdp1(osc.q1+k1q1*0.5, osc.q2+k1q2*0.5, osc.p1+k1p1*0.5, osc.p2+k1p2*0.5, osc.m, osc.L)
        k2q2 = dt * dHdp2(osc.q1+k1q1*0.5, osc.q2+k1q2*0.5, osc.p1+k1p1*0.5, osc.p2+k1p2*0.5, osc.m, osc.L)

        k3p1 = dt * (-dHdq1(osc.q1+k2q1*0.5, osc.q2+k2q2*0.5, osc.p1+k2p1*0.5, osc.p2+k2p2*0.5, osc.m, osc.L))
        k3p2 = dt * (-dHdq2(osc.q1+k2q1*0.5, osc.q2+k2q2*0.5, osc.p1+k2p1*0.5, osc.p2+k2p2*0.5, osc.m, osc.L))
        k3q1 = dt * dHdp1(osc.q1+k2q1*0.5, osc.q2+k2q2*0.5, osc.p1+k2p1*0.5, osc.p2+k2p2*0.5, osc.m, osc.L)
        k3q2 = dt * dHdp2(osc.q1+k2q1*0.5, osc.q2+k2q2*0.5, osc.p1+k2p1*0.5, osc.p2+k2p2*0.5, osc.m, osc.L)

        k4p1 = dt * (-dHdq1(osc.q1+k3q1, osc.q2+k3q2, osc.p1+k3p1, osc.p2+k3p2, osc.m, osc.L))
        k4p2 = dt * (-dHdq2(osc.q1+k3q1, osc.q2+k3q2, osc.p1+k3p1, osc.p2+k3p2, osc.m, osc.L))
        k4q1 = dt * dHdp1(osc.q1+k3q1, osc.q2+k3q2, osc.p1+k3p1, osc.p2+k3p2, osc.m, osc.L)
        k4q2 = dt * dHdp2(osc.q1+k3q1, osc.q2+k3q2, osc.p1+k3p1, osc.p2+k3p2, osc.m, osc.L)

        osc.p1 = osc.p1 + (1/6)*(k1p1 + 2*k2p1 + 2*k3p1 + k4p1)
        osc.p2 = osc.p2 + (1/6)*(k1p2 + 2*k2p2 + 2*k3p2 + k4p2)
        osc.q1 = osc.q1 + (1/6)*(k1q1 + 2*k2q1 + 2*k3q1 + k4q1)
        osc.q2 = osc.q2 + (1/6)*(k1q2 + 2*k2q2 + 2*k3q2 + k4q2)

        osc.t += dt

        if osc.p2 > 0 and osc.q2**2 < 0.005**2:     # due to q2=0 makes the pendelum too damped
            obs.poincare_q1.append(osc.q1)
            obs.poincare_p1.append(osc.p1)
        


# Animation function which integrates a few steps and return a line for the pendulum
def animate(framenr, osc, obs, integrator, pendulum_lines, stepsperframe):
    for it in range(stepsperframe):
        integrator.integrate(osc, obs)

    x1 = math.sin(osc.q1);
    y1 = -math.cos(osc.q1);
    x2 = x1 + math.sin(osc.q2);
    y2 = y1 - math.cos(osc.q2)
    pendulum_lines.set_data([0, x1, x2], [0, y1, y2])
    return pendulum_lines,


class Simulation:

    def run(self,
            integrator,
            tmax=50.,  # final time
            stepsperframe=50,  # how many integration steps between visualising frames
            outfile='energy1.pdf'
            ):
        numframes = int(tmax / (stepsperframe * integrator.dt))

        plt.clf()
        fig = plt.figure()
        ax = plt.subplot(xlim=(-2.2, 2.2), ylim=(-2.2, 2.2))
        plt.axhline(y=0)  # draw a default hline at y=1 that spans the xrange
        plt.axvline(x=0)  # draw a default vline at x=1 that spans the yrange
        pendulum_lines, = ax.plot([], [], lw=5)

        oscillator = Oscillator()
        obs = Observables()
        # Call the animator, blit=True means only re-draw parts that have changed
        anim = animation.FuncAnimation(fig, animate,  # init_func=init,
                                       fargs=[oscillator, obs, integrator, pendulum_lines, stepsperframe],
                                       frames=numframes, interval=25, blit=True, repeat=False)
        plt.show()

        # If you experience problems visualizing the animation and/or
        # the following figures comment out the next line 
        plt.waitforbuttonpress(30)

##        plt.figure()
##        plt.xlabel('q1')
##        plt.ylabel('p1')
##        plt.plot(obs.q1list, obs.p1list)
##        plt.tight_layout()  # adapt the plot area tot the text with larger fonts 
##        plt.show()
##    
##        plt.figure()
##        plt.xlabel('q2')
##        plt.ylabel('p2')
##        plt.plot(obs.q2list, obs.p2list)
##        plt.tight_layout()  # adapt the plot area tot the text with larger fonts 
##        plt.show()
        
        plt.figure()
        plt.xlabel('q1')
        plt.ylabel('p1')
        plt.plot(obs.poincare_q1, obs.poincare_p1, 'ro')
        plt.tight_layout()  # adapt the plot area tot the text with larger fonts
        plt.show()

##        plt.figure()
##        plt.xlabel('time')
##        plt.ylabel('energy')
##        plt.plot(obs.time, obs.epot, obs.time, obs.ekin, obs.time, obs.etot)
##        plt.legend(('Epot', 'Ekin', 'Etot'))
##        plt.tight_layout()  # adapt the plot area tot the text with larger fonts 
##        plt.show()
## 

simulation = Simulation()
simulation.run(integrator=RK4Integrator())
