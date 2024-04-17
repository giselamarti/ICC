"""
# PROBLEMA 6.2 - AUTOCATALISI
# Gisela Martí Guerrero
"""

import numpy as np
import matplotlib.pyplot as plt

k1 = 1e-4
k2 = 8
xo = 1e-3
yo = 0
dt = 1          
time = np.arange(0,600,dt)
methods = ["Simple", "Modified","Improved"]

def fx(x,y): 
    return -2*k1*x-2*k2*x*y

def fy(y,x):
    return 2*k1*x+2*k2*x*y

def euler(y,x,f,dt,mode):
    if mode.lower() == "simple": 
        a,b,d,g = 1,0,0,0
    elif mode.lower() == "modified": 
        a,b,d,g = 0,1,0.5,0.5
    elif mode.lower() == "improved":
        a,b,d,g = 0.5,0.5,1,1

    yt = y + dt*(a*f(y,x) + b*f(y+f(y,x)*d*dt,x))
    return yt

conc_x = [[xo] for _ in range(3)]
conc_y = [[yo] for _ in range(3)]
for t in time:
    for i,method in enumerate(methods):
        xi = conc_x[i][-1]
        yi = conc_y[i][-1]
        x_next = euler(xi,yi,fx,dt,mode=method)
        y_next = euler(yi,xi,fy,dt,mode=method)
        conc_x[i].append(x_next)
        conc_y[i].append(y_next)
        
conc_x = [cx[:-1] for cx in conc_x]
conc_y = [cy[:-1] for cy in conc_y]

output_data = np.column_stack((time, np.round(conc_x[0], 6), np.round(conc_y[0], 6),
                               np.round(conc_x[1], 6), np.round(conc_y[1], 6),
                               np.round(conc_x[2], 6), np.round(conc_y[2], 6)))
header = "Time, Simple_[MnO4-], Simple_[Mn2+], Modified_[MnO4-], Modified_[Mn2+], Improved_[MnO4-], Improved_[Mn2+]"
np.savetxt("6.2_output_concentrations.txt", output_data, header=header, delimiter='\t', fmt='%.6f')

plt.figure(figsize=(10, 6))
for i, method in enumerate(methods):
    plt.plot(time, conc_x[i], label=f'{method} [MnO4-]')
    plt.plot(time, conc_y[i], label=f'{method} [Mn2+]')
plt.title('Concentration Evolution Over Time')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (M)')
plt.xlim([0,600])
plt.ylim([-0.0001,0.0011])
plt.legend()
plt.grid(True)
plt.show()