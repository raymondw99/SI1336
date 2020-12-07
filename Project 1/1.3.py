import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp
from scipy.optimize import curve_fit

dt = 0.001
omega_0 = 2
tmax = 20

def plot_labels(x, y):
    plt.xlabel(x)
    plt.ylabel(y)
    plt.show()

def acc(x, dx, gamma):
    return -omega_0**2 * x - gamma * dx

def energy(x, dx):
    return 0.5 * omega_0**2 * x**2 + 0.5 * dx**2

def velocity_verlet(acc, x, dx, gamma):
    a = acc(x, dx, gamma)
    x += dx * dt + 0.5 * a * dt * dt
    a_next = acc(x, dx, gamma)
    dx += 0.5 * (a_next + a) * dt

    return [x, dx]

def curve(gamma):
    x = 1
    dx = 0
  
    x_values = []
    dx_values = []
    energies = []
    times = []

    t = 0
    
    while t < tmax:
        x_values.append(x)
        dx_values.append(dx)
        energies.append(energy(x, dx))
        times.append(t)
        
        t += dt
        [x, dx] = velocity_verlet(acc, x, dx, gamma)
        
    return [times, gamma, x_values, dx_values, energies]

def graph(n, y, title):
    
    values = []
    gamma_values = [0.5,1,2,3]
    for gamma in gamma_values:
      values.append(curve(gamma))
  
    plt.figure()
    labels = []
    for (times, gamma, x_values, dx_values, energies) in values:
        labels.append("γ = " + str(gamma))

    for value in values:
        plt.plot(np.array(value[0]), np.array(value[n]))
    
    plt.xlim(0, tmax)
    plt.title(title)
    plt.figlegend(labels)
    plt.xlabel("t", fontsize = 14)
    plt.ylabel(y, fontsize = 14)
    plt.show()

def envelope(times, signal):
    try:
        vertex, cov = sp.find_peaks(signal)
        vertex_times = np.insert(times[vertex], 0, 0)
        vertex_points = np.insert(signal[vertex], 0, signal[0])

        def func(x, A, T):
            return A*np.exp(-T*x)

        (A, T), cov = curve_fit(func, vertex_times, vertex_points, [1, 0.1])

        return [func(times, A, T), 1/T]
    except:
        return [0, 0]

def analyze():
    relaxation_times = []
    gamma_values = np.linspace(0.1, 12, 40)

    for gamma in gamma_values:
        (times, g, x_values, dx_values, energies) = curve(gamma)
        enveloped_times, relaxation_time = envelope(np.array(times), 
                                           np.array(x_values))
        relaxation_times.append(np.amin(relaxation_time))

    plt.figure()
    plt.title("Relaxation time", fontsize = 12)
    plt.plot(gamma_values, relaxation_times)
    plt.xlabel("γ", fontsize = 12)
    plt.ylabel("Relaxation time", fontsize = 12)
    plt.show()

graph(2, "θ(t)", "Position over time")
graph(3, "θ'(t)","Velocity over time")
graph(4, "E(t)", "Energy over time")
analyze()
