import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import random as rnd
import math


def runRandomWalk(steps, canIntersect, preventBackfire):
    X = []
    Y = []
    indices = []

    x = 0
    y = 0

    direction = math.floor(rnd.uniform(0, 4))

    def intersectsSelf():
        for i in range(len(X)):
            if X[i] == x and Y[i] == y:
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
        if not canIntersect and intersectsSelf():
            return None

        if preventBackfire:
            newDir = math.ceil(rnd.uniform(0, 3))
            direction = (direction - 2 + newDir) % 4
        else:
            direction = math.floor(rnd.uniform(0, 4))

        X.append(x)
        Y.append(y)
        indices.append(i)
    return (np.array(X), np.array(Y), np.array(indices))


def runUntilWorks(steps, canIntersect, preventBackfire):
    stepsRun = 1
    while True:
        val = runRandomWalk(steps, canIntersect, preventBackfire)
        if preventBackfire == False:
            return val, 1

        if val != None:
            return val, stepsRun
        stepsRun += 1


def runWithChoice(simLength, canIntersect, preventBackfire):
    x = np.int_(np.linspace(5, simLength, 5))
    y = np.int_(np.linspace(10, 1000, 10))
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    for i in range(0, len(x)):
        xval = x[i]
        for j in range(0, len(y)):
            yval = y[j]

            summedSquares = 0
            for sims in range(0, yval):
                (Xv, Yv, indices), stepsRun = runUntilWorks(
                    xval, canIntersect, preventBackfire)
                r = (Xv[0]-Xv[-1])**2 + (Yv[0]-Yv[-1])**2
                summedSquares += r / yval

            Z[j, i] = summedSquares
    print()

def visualize():
    R2_intersecting = [math.sqrt(i) for i in [7.043999999999888,47.05199999999999,
                          105.06000000000003,173.5080000000005,
                          262.29200000000066]]
    R2_non_intersecting = [math.sqrt(i) for i in [4.15399999999994, 14.172000000000061,
                             25.417999999999935,35.272000000000055,
                             49.76400000000004]]

    length = [5,15,25,35,50]



    
    plt.figure()
    plt.title("Distance comparison")
    plt.plot(length, R2_intersecting,
             color = 'tab:blue',label = "Intersecting")
    plt.plot(length, R2_non_intersecting,
             color = 'tab:cyan',label = "Non-intersecting")
    plt.xlabel('N (Steps)', fontsize = 12)
    plt.ylabel('$\sqrt{⟨R^2⟩}$', fontsize = 12)
    plt.legend()
    plt.show()

    plt.figure()
    plt.title("Log-log distance comparison")
    plt.plot([math.log1p(i) for i in length],
             [math.log1p(i) for i in R2_intersecting],
             color = 'tab:red',label = "Intersecting")
    plt.plot([math.log1p(i) for i in length],
             [math.log1p(i) for i in R2_non_intersecting],
             color = 'tab:orange',label = "Non-intersecting")
    plt.xlabel('ln(N) (Steps)', fontsize = 12)
    plt.ylabel('$ln(\sqrt{⟨R^2⟩})$', fontsize = 12)
    plt.legend() 
    plt.show()


visualize() 
#runWithChoice(50, False, True)
#runWithChoice(50, True, False)
