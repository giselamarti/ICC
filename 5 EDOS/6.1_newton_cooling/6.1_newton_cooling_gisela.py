"""
# PROBLEMA 6.1 - NEWTON COOLING
# Gisela Martí Guerrero
"""

import numpy as np
import matplotlib.pyplot as plt

K = 0.05   # K/min
dt = 0.1  
time = np.arange(1,100,dt)
methods = ["Simple", "Modified","Improved"]

Tf = int(input("fixed temperature?"))
To = int(input("initial temperature?"))

def euler(x,f,dt,mode):
    if mode.lower() == "simple": 
        a,b,d,g = 1,0,0,0
    elif mode.lower() == "modified": 
        a,b,d,g = 0,1,0.5,0.5
    elif mode.lower() == "improved":
        a,b,d,g = 0.5,0.5,1,1

    xt = x + dt*(a*f(x) + b*f(x+f(x)*d*dt))
    return xt

def f(T): 
    return K*(Tf-T)

temperatures = [[To] for _ in range(3)]
for t in time:
    for i, method in enumerate(methods):
        Ti = temperatures[i][-1]
        T_next = euler(Ti,f,dt,method)
        temperatures[i].append(T_next)
        
temperatures = [temp[:-1] for temp in temperatures]
T_exact = Tf + (To-Tf)*np.exp(-K*time)

output_data = np.column_stack((time, temperatures[0], temperatures[1], temperatures[2], T_exact))
header = "Time,     Simple,     Modified,   Improved,  T_exact"
np.savetxt("6.1_output_temperatures.txt", output_data, header=header, delimiter='\t', fmt='%.6f')

plt.figure(figsize=(10, 6))
for i, method in enumerate(methods):
    plt.plot(time, temperatures[i], label=f"{method} Method")

plt.plot(time, T_exact, label="Exact Solution", linestyle='--', color='black')

plt.title('Temperature vs Time')
plt.xlabel('Time (min)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.grid(True)
plt.show()