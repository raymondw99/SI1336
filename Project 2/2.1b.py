import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import math



def run_random(steps, r0, a, c, m):


    def random_generator():
        nonlocal r0 

        r0 = (a * r0 + c) % m
        return r0 / (m - 1)

    x_vals = []
    y_vals = []

    x = 0
    y = 0

    for i in range(steps):
        r = (a * r0 + c) % m
        direction = r%4

        if direction == 0:
            x += 1

        elif direction == 1:
            y += 1

        elif direction == 2:
            x -= 1

        elif direction == 3:
            y -= 1

        x_vals.append(x)
        y_vals.append(y)
    return x_vals, y_vals


def plot_walk(r0, a, c, m):
    plt.figure()
    plt.title("1000 steps (r = {0}, a = {1}, c = {2}, m = {3})".format(r0, a, c, m))
    x_vals, y_vals = run_random(1000, r0, a, c, m)
    plt.plot(x_vals, y_vals)
    plt.show()
    
plot_walk(1, 3, 4, 128)
plot_walk(1, 3, 4, 129)
plot_walk(1, 3, 4, 130)

plot_walk(1, 4, 4, 128)
plot_walk(1, 5, 4, 128)
plot_walk(1, 6, 4, 128)

plot_walk(1, 3, 5, 128)
plot_walk(1, 3, 6, 128)
plot_walk(1, 3, 7, 128)

plot_walk(1, 3, 4, 128)
plot_walk(2, 3, 4, 128)
plot_walk(3, 3, 4, 128)
