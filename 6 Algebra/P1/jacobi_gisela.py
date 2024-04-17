import numpy as np
from pandas import *

def maximum_element(a):
    max_el = 0.0
    n = len(a)
    for i in range(n - 1):
        for j in range(i + 1, n):
            if abs(a[i,j]) >= max_el:
                max_el = abs(a[i,j])
                ip = i
                iq = j
    return max_el, ip, iq 

def rotation(a, v, ip, iq):
    n = len(a)
    h = a[iq,iq] - a[ip,ip]
    if abs(a[ip,iq]) < abs(h)*1.0e-36: # Machine epsilon ~ 1e-36
        t = a[ip,iq]/h
    else:
        theta = h/(2.0*a[ip,iq])
        t = 1.0/(abs(theta) + np.sqrt(1.0 + theta**2))
        if theta < 0.0: 
            t = -t
    c = 1.0/np.sqrt(1.0 + t**2)
    s = t*c                     
    tau = s/(1.0 + c)
    temp = a[ip,iq] #h
    a[ip,iq] = 0.0
    a[ip,ip] = a[ip,ip] - t*temp
    a[iq,iq] = a[iq,iq] + t*temp
    for i in range(ip):       # Case of rotations 1 < i < p
        temp = a[i,ip]
        a[i,ip] = temp - s*(a[i,iq] + tau*temp)
        a[i,iq] = a[i,iq] + s*(temp - tau*a[i,iq])
    for i in range(ip+1,iq):  # Case of rotations p < i < q
        temp = a[ip,i]
        a[ip,i] = temp - s*(a[i,iq] + tau*a[ip,i])
        a[i,iq] = a[i,iq] + s*(temp - tau*a[i,iq])
    for i in range(iq+1,n):   # Case of n > i > q
        temp = a[ip,i]
        a[ip,i] = temp - s*(a[iq,i] + tau*temp)
        a[iq,i] = a[iq,i] + s*(temp - tau*a[iq,i])
    for i in range(n):
        temp = v[i,ip]
        v[i,ip] = temp - s*(v[i,iq] + tau*v[i,ip])
        v[i,iq] = v[i,iq] + s*(temp - tau*v[i,iq])
        
def jacobi(a):
    n = len(a)
    max_rot = 10*(n**2)
    v = np.identity(n)*1.0
    for i in range(max_rot):
        max_a,ip,iq = maximum_element(a)
        if max_a < tolerance: 
            return np.diagonal(a),v
        rotation(a,v,ip,iq)
    
def hilbert(n):
    H_matrix = np.array([[1.0/(i + j - 1) for j in range(1, n+1)] for i in range(1, n+1)])
    return H_matrix

# Define the tolerance
tolerance = 1.0e-9

# Hilbert matrix
N = 10
HM = hilbert(N)

# PAP and QAQ quadrants
PAP = HM[:int(N/2),:int(N/2)]
QAQ = HM[int(N/2):,int(N/2):]

cols = ['c{}'.format(i) for i in range(1, 11)]
rows = ['r{}'.format(i) for i in range(1, 11)]
print("Input matrix:")
print(DataFrame(HM, columns=cols, index=rows))
print("")
print("PAP:")
print(PAP)
print("")
print("QAQ:")
print(QAQ)

# Obtain the eigenvalues and the transformation matrix
eigenval, v = jacobi(HM.copy()) 

# Check if the eigenvalue problem is satisfied
print("Is the eigenvalue problem satisfied?: ",np.all(np.isclose((HM@v),(eigenval*v))))
print("")
print("Eigenvalues:")
print(eigenval)
print("")
print("Transformation matrix:")
print(DataFrame(v, columns=cols, index=rows))

# Effective Hamiltonian matrix
eff_H = v.T@HM@v
print("Effective Hamiltonian:")
print(DataFrame(eff_H, columns=cols, index=rows))