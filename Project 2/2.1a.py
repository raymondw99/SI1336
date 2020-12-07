import numpy as np
import matplotlib.pyplot as plt
import random 
import math

x_vals = []
y_vals = []

x = 0
y = 0
steps = 100

for i in range(steps):
    val = math.floor(random.uniform(0, 4))

    if val == 0:
        x += 1

    if val == 1:
        y += 1

    if val == 2:
        x -= 1

    if val == 3:
        y -= 1

    x_vals.append(x)
    y_vals.append(y)

plt.figure()
plt.title(str(steps) + " steps")
plt.plot(x_vals, y_vals)
plt.rc('xtick',labelsize=1)
plt.rc('ytick',labelsize=1)
plt.show()


