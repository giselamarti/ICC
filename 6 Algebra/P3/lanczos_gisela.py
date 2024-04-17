import numpy as np
from pandas import *
import numpy.linalg as la

def lanczos(H,vg):
    Lv=np.zeros((len(vg),len(vg)), dtype=complex)
    Hk=np.zeros((len(vg),len(vg)), dtype=complex)
    Lv[0]=vg/la.norm(vg)
    w=np.dot(H,Lv[0]) 
    a=np.dot(np.conj(w),Lv[0])
    w=w-a*Lv[0]
    Hk[0,0]=a
    for j in range(1,len(vg)):
        b=(np.dot(np.conj(w),np.transpose(w)))**0.5
        Lv[j]=w/b
        w=np.dot(H,Lv[j])
        a=np.dot(np.conj(w),Lv[j])
        w=w-a*Lv[j]-b*Lv[j-1]
        Hk[j,j]=a
        Hk[j-1,j]=b
        Hk[j,j-1]=np.conj(b)
        
    return (Hk,Lv)

def lanczos(H,vg):
    Lv=np.zeros((len(vg),len(vg)), dtype=complex)
    Hk=np.zeros((len(vg),len(vg)), dtype=complex)
    Lv[0]=vg/la.norm(vg)
    w=np.dot(H,Lv[0]) 
    a=np.dot(np.conj(w),Lv[0])
    w=w-a*Lv[0]
    Hk[0,0]=a
    for j in range(1,len(vg)):
        b=(np.dot(np.conj(w),np.transpose(w)))**0.5
        Lv[j]=w/b
        w=np.dot(H,Lv[j])
        a=np.dot(np.conj(w),Lv[j])
        w=w-a*Lv[j]-b*Lv[j-1]
        Hk[j,j]=a
        Hk[j-1,j]=b
        Hk[j,j-1]=np.conj(b)
        
    return (Hk,Lv)

def QR(A):
    n, m = A.shape
    Q = np.empty((n, n))
    u = np.empty((n, n))
    u[:, 0] = A[:, 0]
    Q[:, 0] = u[:, 0] / np.linalg.norm(u[:, 0])

    for i in range(1, n):
        u[:, i] = A[:, i]
        for j in range(i):
            u[:, i] -= (A[:, i] @ Q[:, j]) * Q[:, j]
        Q[:, i] = u[:, i] / np.linalg.norm(u[:, i])
    R = np.zeros((n, m))
    for i in range(n):
        for j in range(i, m):
            R[i, j] = A[:, j] @ Q[:, i]

    return Q, R

def eigenvalues(A):
    A_old = np.copy(A)
    A_new = np.copy(A)
    diff = np.inf
    i = 0
    while (diff > tolerance) and (i < max_iter):
        A_old[:, :] = A_new
        Q, R = QR(A_old)
        A_new[:, :] = R @ Q
        diff = np.abs(A_new - A_old).max()
        i += 1
    eigvals = np.diag(A_new)
    return eigvals

def hilbert(n):
    H_matrix = np.array([[1.0 / (i+j-1) for j in range(1, n+1)] for i in range(1, n+1)])
    return H_matrix

# Define parameters
tolerance = 1.0e-9
max_iter = 1000

# Hilbert matrix
N = 100
HM = hilbert(N)

print("Input matrix:")
cols = ['c{}'.format(i) for i in range(1, 101)]
rows = ['r{}'.format(i) for i in range(1, 101)]
print(DataFrame(HM, columns=cols, index=rows))

v_rand = np.random.random(N)
T,V = lanczos(HM.copy(),v_rand)
eigenval = eigenvalues(T.real)

print("Tridiagonal matrix:")
print(DataFrame(T.real, columns=cols, index=rows))
print("")
print("Eigenvalues:")
print(eigenval)
