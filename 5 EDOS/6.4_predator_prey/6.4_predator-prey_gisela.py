"""
# PROBLEMA 6.4 - PREY-PREDATOR MODEL
# Gisela Martí Guerrero
"""

import numpy as np
import matplotlib.pyplot as plt

x_o = 1000    # conills
y_o = 40      # llops
a = 0.04      # alfa
b = 0.0005    # beta
g = 0.3       # gamma
d = 5e-5      # delta
dt = 0.01
time = np.arange(0,1000,dt)

def fx(x,y):
    return a*x - b*x*y

def fy(y,x):
    return -g*y + d*x*y

def RK4(x,y,dt,f):
    f0 = f(x,y)
    f1 = f(x+f0*dt/2, y)
    f2 = f(x+f1*dt/2, y)
    f3 = f(x+f2*dt, y)
    xt = x+dt/6*(f0+2*f1+2*f2+f3)
    return xt

# Main integration loop
X = [x_o]
Y = [y_o]
for t in time:
    x_i = X[-1]
    y_i = Y[-1]
    x_next = RK4(x_i,y_i,dt,fx)
    y_next = RK4(y_i,x_i,dt,fy)
    X.append(x_next)
    Y.append(y_next)

output_data = np.column_stack((time, X[:-1], Y[:-1]))
header = "Time,    X (rabbits), Y (wolves)"
np.savetxt("6.4_predator_prey.txt", output_data, header=header, delimiter='\t', fmt='%.6f')

plt.figure(figsize=(10, 6))
plt.plot(time,X[:-1], label="Rabbits", c="black")
plt.plot(time,Y[:-1], label="Wolves", c="red")
plt.title('Time Evolution of the predator-prey model')
plt.xlabel('Time')
plt.ylabel('Population')
plt.xlim([0,1000])
plt.legend(loc="upper left")
plt.tight_layout()
plt.show()

plt.figure(figsize=(7, 7))
plt.plot(X,Y, c="black")
plt.title('Parametric relationship for the predator-prey model')
plt.xlabel('X (rabbits)')
plt.ylabel('Y (wolves)')
plt.tight_layout()
plt.show()