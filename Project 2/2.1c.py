import random as rnd
from pylab import *
import numpy as np
from scipy.stats import sem

def atask(N=1000):
    x, y = [0], [0] #Initial position
    for nstep in range(1,N):
        step = int(rnd.random()*4)
        x.append(x[nstep-1])
        y.append(y[nstep-1])
        if step == 0:
            #print(step)
            x[nstep]+=1
        elif step == 1:
            x[nstep]-=1
        elif step == 2:
            y[nstep]+=1
        elif step == 3:
            y[nstep]-=1
            #print(step)
        else:
            raise Exception('Something is wrong')
    return x,y

meansqrt = []
meansqrtfluct = []
mean = []
stderror = []
iterations=5000
N = arange(1,1002,50)
for n in N:
    distances = []
    sumquad = 0
    sumsqrt = 0
    for i in range(iterations):
        x,y = atask(n)
        R2 = x[-1]**2+y[-1]**2
        sumquad += R2
        sumsqrt += sqrt(R2)
        distances.append(sqrt(R2))
    meansqrt.append(sqrt(sumquad/iterations))
    mean.append(sumquad/iterations)
    meansqrtfluct.append(sqrt((sumquad/iterations-(sumsqrt/iterations)**2)*iterations/(iterations-1)))
    stderror.append(sem(distances))
print(N[-1])

plt.figure()
plt.title('Comparison')
plt.plot(N,meansqrt,label=r'$\sqrt{<R^2>}$')
plt.plot(N,meansqrtfluct,label="RMSF")
plt.plot(N,stderror, label="STDE")
plt.xlabel('N (Steps)', fontsize = 12)
#plt.ylabel(fontsize=15)
plt.legend()
plt.show()


plt.figure()
plt.title('Length dependence on N')
plt.plot(N,mean)
plt.xlabel('N (Steps)', fontsize = 12)
plt.ylabel('Length $⟨R^2⟩$', fontsize = 12)
plt.show()

plt.figure()
plt.title('Standard error')
plt.plot(N,stderror)
plt.xlabel('N (Steps)', fontsize = 12)
plt.ylabel('Standard error', fontsize = 12)
plt.show()

