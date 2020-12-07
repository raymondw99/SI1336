# Python simulation of simple planar pendulum with real time animation
# BH, OF 2020-10-20, latest verson 2020-10-28

from matplotlib import animation
from pylab import *

# Global constants
G = 9.8  # gravitational acceleration


class Oscillator:
    """Class for an general, simple oscillator"""
    m = 1  # mass of the pendulum bob
    c = 4  # c=g/L=4
    L = G / c  # string length
    omega_0 = 2
    gamma = 0.1

    t = 0  # the time
    theta = np.pi*0.1  # the position/angle
    dtheta = 0  # the velocity


# Class for storing observables for an oscillator
class Observables:
    def __init__(self):
        self.time = []  # list to store time
        self.pos = []  # list to store positions
        self.vel = []  # list to store velocities
        self.energy = []  # list to store energy


class BaseSystem:
    def force(self, osc):
        """ Implemented by the childclasses  """
        pass


class Harmonic(BaseSystem):
    def force(self, osc):
        
        return -osc.c * osc.theta


class Pendulum(BaseSystem):
    def force(self, osc):
        return -osc.c * np.sin(osc.theta)


class BaseIntegrator:
    dt = 0.01 # time step

    def integrate(self, simsystem, osc, obs):
        """ Perform a single integration step """
        self.timestep(simsystem, osc, obs)

        # Append observables to their lists
        obs.time.append(osc.t)
        obs.pos.append(osc.theta)
        obs.vel.append(osc.dtheta)
        obs.energy.append(
            0.5 * osc.m * osc.L ** 2 * osc.dtheta ** 2 + 0.5 * osc.m * G * osc.L * osc.theta ** 2)  # harmonic oscillator        

    def timestep(self, simsystem, osc, obs):
        """ Implemented by the childclasses """
        pass


class EulerIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        osc.t += self.dt

        accel = simsystem.force(osc) / osc.m
        
        # TODO: Implement the integration here, updating osc.theta and osc.dtheta
        c = osc.theta #Temporarily storing the value of theta at time t
        osc.theta += osc.dtheta*self.dt
        osc.dtheta -= (G/osc.L)*c*self.dt

        

class EulerCromerIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        osc.t += self.dt

        accel = simsystem.force(osc) / osc.m

        # TODO: Implement the integration here, updating osc.theta and osc.dtheta
        c = osc.dtheta
        osc.dtheta -= (G/osc.L)*osc.theta*self.dt
        osc.theta += osc.dtheta*self.dt

class VerletIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        osc.t += self.dt

        accel = simsystem.force(osc) / osc.m

        # TODO: Implement the integration here, updating osc.theta and osc.dtheta
        osc.theta = osc.theta + osc.dtheta*self.dt + 0.5*accel*(self.dt**2)
        osc.dtheta = osc.dtheta + 0.5*self.dt*(simsystem.force(osc)/osc.m + accel) 
        
class RungeKutta4(BaseIntegrator):
    def timestep(self,simsystem,osc,obs):
        
        accel1 = -osc.c * np.sin(osc.theta)
        a1 = accel1*self.dt
        b1 = osc.dtheta*self.dt
        osc.t += self.dt

        accel2 = -osc.c * np.sin(osc.theta + 0.5*b1)
        a2 = accel2*self.dt
        b2 = (osc.dtheta + 0.5*a1)*self.dt
        osc.t += self.dt + 0.5*self.dt
        
        accel3 = -osc.c * np.sin(osc.theta + 0.5*b1 + 0.5*b2)
        a3 = accel3*self.dt
        b3 = (osc.dtheta + 0.5*a1 + 0.5*a2)*self.dt
        osc.t += self.dt + 0.5*self.dt + 0.5*self.dt
        
        accel4 = -osc.c * np.sin(osc.theta + 0.5*b1 + 0.5*b2 + b3)
        a4 = accel4*self.dt
        b4 = (osc.dtheta + 0.5*a1 + 0.5*a2 + a3)*self.dt
        osc.t += self.dt + 0.5*self.dt + 0.5*self.dt + self.dt

        osc.dtheta += (1/6)*(a1+2*a2+2*a3+a4)
        osc.theta += (1/6)*(b1+2*b2+2*b3+b4)


# Animation function which integrates a few steps and return a line for the pendulum
def animate(framenr, simsystem, oscillator, obs, integrator, pendulum_line, stepsperframe):
    for it in range(stepsperframe):
        integrator.integrate(simsystem, oscillator, obs)

    x = np.array([0, np.sin(oscillator.theta)])
    y = np.array([0, -np.cos(oscillator.theta)])
    pendulum_line.set_data(x, y)
    return pendulum_line,


class Simulation:

    def run(self,
            simsystem,
            integrator,
            tmax= 100,  # final time
            stepsperframe=100,  # how many integration steps between visualising frames
            title="simulation",  # Name of output file and title shown at the top
            ):
        oscillator = Oscillator()
        obs = Observables()

        numframes = int(tmax / (stepsperframe * integrator.dt))

        plt.clf()
        # fig = plt.figure()
        ax = plt.subplot(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
        plt.axhline(y=0)  # draw a default hline at y=1 that spans the xrange
        plt.axvline(x=0)  # draw a default vline at x=1 that spans the yrange
        pendulum_line, = ax.plot([], [], lw=5)
        plt.title(title)
        # Call the animator, blit=True means only re-draw parts that have changed
        anim = animation.FuncAnimation(plt.gcf(), animate,  # init_func=init,
                                       fargs=[simsystem, oscillator, obs, integrator, pendulum_line, stepsperframe],
                                       frames=numframes, interval=25, blit=True, repeat=False)
        #plt.show()

        # If you experience problems visualizing the animation and/or
        # the following figures comment out the next line 
        plt.waitforbuttonpress(30)

        plt.clf()
        #plt.title("Harmonic-RungeKutta4 - (0)/π = 0.1, dt = 0.001", fontsize = 14)
        plt.title('Analytical versus numerical')
        analytical = np.multiply(np.cos(np.multiply(obs.time,2)),np.pi*0.1)
        plt.plot(obs.time, analytical - obs.pos)
        #plt.plot(obs.time, obs.pos, label="Position")
        #plt.plot(obs.time, obs.vel, label="Velocity")
        #plt.plot(obs.time, obs.energy, label="Energy")
        plt.xlabel('time', fontsize = 14)
        plt.ylabel('difference', fontsize = 14)
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        #plt.ylabel('energy', fontsize = 14)
        plt.legend()
        plt.show()


simulation = Simulation()
#simulation.run(simsystem=Harmonic(), integrator=VerletIntegrator(), title="Analytical-Numerical")
#simulation.run(simsystem=Harmonic(), integrator=RungeKutta4(), title="Harmonic-RungeKutta4")
#simulation.run(simsystem=Harmonic(), integrator=EulerIntegrator(), title="Harmonic-Euler")
#simulation.run(simsystem=Harmonic(), integrator=EulerCromerIntegrator(), title="Harmonic-EulerCromer")
simulation.run(simsystem=Harmonic(), integrator=VerletIntegrator(), title="Harmonic-Verlet")
#simulation.run(simsystem=Pendulum(), integrator=RungeKutta4(), title="Pendulum-RungeKutta4")
#simulation.run(simsystem=Pendulum(), integrator=EulerCromerIntegrator(), title="Pendulum-EulerCromer")
#simulation.run(simsystem=Pendulum(), integrator=VerletIntegrator(), title="Pendulum-Verlet")

