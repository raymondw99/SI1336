import numpy as np
import matplotlib.pyplot as plt
from pylab import *

g = 9.81
m = 1
k = 9
l = g /k

dt = 0.001

def acc_harmonic(theta):
    return -k * theta

def acc_pendulum(theta):
    return -k * np.sin(theta)

def energy_harmonic(theta, dtheta):
    return 0.5 * m * l**2 * dtheta**2 + 0.5 * m * g * l * theta**2

def energy_pendulum(theta, dtheta):
    return 0.5 * m * l**2 * dtheta**2 + m * g * l * (1. - np.cos(theta))


def velocity_verlet(A, E, theta, dtheta):
    a = A(theta)
    theta += dtheta * dt + 0.5 * a * dt * dt
    a_next = A(theta)
    dtheta += 0.5 * (a_next + a) * dt

    energy = E(theta, dtheta)
    return [theta, dtheta, energy]
          
def period(A, E, theta_i):
    theta = theta_i
    dtheta = 0 
    time = 0 
    
    def iterate():
        nonlocal theta, dtheta, time
        current = sign(theta)
        while sign(theta) == current:
            time += dt
            [theta, dtheta, energy] = velocity_verlet(A, E, theta, dtheta)
    
    iterate()
    start = time
    iterate()
    iterate()
    end = time
    period = end - start
    
    return period

def output():
  theta = np.linspace(0.01, np.pi/2, 20)
  period_pendulum, period_harmonic, period_perturbation = [], [], []
  for i in theta:
      period_harmonic.append(period(acc_harmonic, energy_harmonic, i))
      period_pendulum.append(period(acc_pendulum,energy_pendulum, i))
      period_perturbation.append(2*np.pi * np.sqrt(1/k) * (1 + i**2/16 + 
                                 11* i **4/3072) + 173* i**6/737820)

  plt.figure()
  plt.plot(theta/np.pi, period_pendulum)
  plt.plot(theta/np.pi, period_harmonic)
  plt.plot(theta/np.pi, period_perturbation)
  plt.title("Comparison of Pendulum/Harmonic \n against Perturbation series", x = 0.3)
  plt.figlegend(('Pendulum', 'Harmonic', 'Perturbation series O(n^6)'))
  plt.xlabel("θ(0)/pi (Position)", fontsize = 12)
  plt.ylabel("T (Seconds)", fontsize = 12)
  plt.show()

output()


