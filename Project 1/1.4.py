import numpy as np
import matplotlib.pyplot as plt

dt = 0.001
k = 3
g = 9.81
gamma = 1

def acc(theta, dtheta):
    return -k**2 * np.sin(theta) - gamma * dtheta
    
def velocity_verlet(acc, theta, dtheta):
    a = acc(theta, dtheta)
    theta += dtheta * dt + 0.5 * a * dt * dt
    a_next = acc(theta, dtheta)
    dtheta += 0.5 * (a_next + a) * dt

    return [theta, dtheta]

def output():
    theta = np.pi / 2
    dtheta = 0
    
    times = []
    theta_vals = []
    dtheta_vals = []
    
    t = 0
    
    while t < 15:
        theta_vals.append(theta)
        dtheta_vals.append(dtheta)
        times.append(t)
        t += dt
        [theta, dtheta] = velocity_verlet(acc, theta, dtheta)
      
    plt.figure()
    plt.plot(dtheta_vals, theta_vals)
    plt.title("Space phase portrait: γ = 1, √g/l = 2, θ(0) = π/2, θ'(0) = 0 ",
              fontsize = 12)
    plt.xlabel("θ'(t)", fontsize = 12)
    plt.ylabel("θ(t)", fontsize = 12)
    plt.show()  

output()
