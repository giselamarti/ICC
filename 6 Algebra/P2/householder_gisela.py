import numpy as np
from pandas import *

def householder(A):
    n = A.shape[0]
    v = np.zeros(n)
    u = np.zeros(n)
    z = np.zeros(n)

    for k in range(0, n-2):
        if np.isclose(A[k+1,k], 0.0):
            a = -np.sqrt(np.sum(A[(k+1):,k]**2))
        else:
            a = -np.sign(A[k+1,k]) * np.sqrt(np.sum(A[(k+1):,k]**2))
        r2 = a**2 - a*A[k+1,k]
        v[k] = 0.0
        v[k+1] = A[k+1,k] - a
        v[(k+2):] = A[(k+2):,k]
        u[k:] = 1.0 / r2 * np.dot(A[k:,(k+1):], v[(k+1):])
        z[k:] = u[k:] - np.dot(u,v) / (2.0*r2) * v[k:]

        for l in range(k+1, n-1):
            A[(l+1):,l] = (A[(l+1):,l] - v[l] * z[(l+1):] - v[(l+1):] * z[l])
            A[l,(l+1):] = A[(l+1):,l]
            A[l,l] = A[l,l] - 2*v[l]*z[l]

        A[-1,-1] = A[-1,-1] - 2*v[-1]*z[-1]
        A[k,(k+2):] = 0.0
        A[(k+2):,k] = 0.0
        A[k+1,k] = A[k+1,k] - v[k+1]*z[k]
        A[k,k+1] = A[k+1,k]

    return A

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
N = 10
HM = hilbert(N)

print("Input matrix:")
cols = ['c{}'.format(i) for i in range(1, 11)]
rows = ['r{}'.format(i) for i in range(1, 11)]
print(DataFrame(HM, columns=cols, index=rows))

# Compute the Householder tridiagonal matrix of the Hilbert matrix
HM_TD = householder(HM.copy())

# Compute the eigenvalues
eig = eigenvalues(HM_TD.copy())

print("Householder Tridiagonal matrix:")
print(DataFrame(HM_TD, columns=cols, index=rows))
print("")
print("Eigenvalues: ")
print(eig)
