import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import math


def runRandomWalk(steps, not_backfire):
    x_vals = []
    y_vals = []

    x = 0
    y = 0
    indices = []
    
    direction = math.floor(rnd.uniform(0, 4))

    def intersectsSelf():
        for i in range(len(x_vals)):
            if x_vals[i] == x and y_vals[i] == y:
                return True
        return False

    def step():
        nonlocal x, y

        if direction == 0:
            x += 1

        if direction == 1:
            y += 1

        if direction == 2:
            x -= 1

        if direction == 3:
            y -= 1

    for i in range(steps):
        step()
        if intersectsSelf():
            return None

        if not_backfire:
            newDir = math.ceil(rnd.uniform(0, 3))
            direction = (direction - 2 + newDir) % 4
        else:
            direction = math.floor(rnd.uniform(0, 4))

        x_vals.append(x)
        y_vals.append(y)
        indices.append(i)
    return x_vals, y_vals, np.array(indices)


def runUntilWorks(steps, not_backfire):
    stepsRun = 1
    while True:
        val = runRandomWalk(steps, not_backfire)
        if val != None:
            return val, stepsRun
        stepsRun += 1

for i in range(5, 55, 5):
    plt.figure()
    plt.title("{0} non-intersecting steps".format(i))
    (X, Y, indices), stepsRun = runUntilWorks(i, True)
    plt.plot(X, Y)
    plt.show()


def means(start, end, step, not_backfire):
    indices = []
    means = []
    attempts = []
    for i in range(start, end, step):
        N = 50
        runWalks = 0
        for n in range(N):
            _, stepsRun = runUntilWorks(i, not_backfire)
            runWalks += stepsRun
        indices.append(i)
        means.append(N / runWalks)
        attempts.append(runWalks / N)

    return indices, means, attempts


indices1, means1, counts1 = means(1, 50, 1, True)
indices2, means2, counts2 = means(1, 25, 1, False)

plt.figure()
plt.title("Success rate")
plt.plot(indices1, means1, 'tab:red')
plt.plot(indices2, means2, 'tab:orange')
plt.figlegend(("Improved variant", "Original"))
plt.xlabel("N (Steps)", fontsize = 12)
plt.ylabel("Success ratio", fontsize = 12)
plt.show()

plt.figure()
plt.title("Attempts per success")
plt.plot(indices1, counts1, 'tab:blue')
plt.plot(indices2, counts2, 'tab:cyan')
plt.figlegend(("Improved variant", "Original"))
plt.xlabel("N (Steps)", fontsize = 12)
plt.ylabel("Attempts/success", fontsize = 12)
plt.show()
