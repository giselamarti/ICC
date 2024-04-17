import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from matplotlib.animation import FuncAnimation, PillowWriter

# Define parameters
L = 1.0                      # Length
dL = 100                     # Grid points
dx = L/dL                    # Subintervals
xrange = np.arange(0,L+dx,dx)
t = 10000                    # Total time
tp = 1000                    # Time points
dt = t/(tp-1)                # Subintervals
trange = np.arange(0,t+dt,dt)
D = 1e-5                     # Diffusion coefficient
k = 1e-3                     # Velocity constant
animation_frames = 1000
animation_name = "Diffusion.gif"

# Create tridiagonal matrix with scipy
dim = dL+1
a = D * dt / dx**2 
main_diagonal = 1 + 2 * a
off_diagonal = -a
matrix = diags([off_diagonal, main_diagonal, off_diagonal], [-1, 0, 1], shape=(dim, dim), format='csr')
# Set boundary conditions
matrix[0, 0] = 1
matrix[-1, -1] = 1
matrix[0, 1] = 0
matrix[-1, -2] = 0

# Concentration arrays
conc = np.zeros(shape=(2, dL+1))
_, x_points = conc.shape
conc[0][0:int(x_points/2)] += 1
conc[1][int(x_points/2):] += 1
c_init = conc.copy()

# Main integration loop
ct = [[],[]]
alpha = []

for i,ti in enumerate(trange):

    # Periodic conditions
    conc[0][0] = conc[0][-1] = conc[0][-2]
    conc[1][0] = conc[1][-1] = conc[1][-2]

    # Save frames
    if i%int(tp/animation_frames) == 0 : 
        ct[0].append(conc[0].copy())
        ct[1].append(conc[1].copy())

        ai = -np.log(conc[0].mean()) / np.log(ti)
        alpha.append(ai)

    # Update concentrations
    conc[0] = np.linalg.inv(matrix.toarray())@conc[0] - dt*k*conc[0]*conc[1]
    conc[1] = np.linalg.inv(matrix.toarray())@conc[1] - dt*k*conc[1]*conc[1]

ct = np.array(ct,dtype=object)
alpha = np.array(alpha)

# Start animation
def GIF(frame):
    """Function that creates a frame for the GIF."""
    ax.clear()
    cA0, = ax.plot(xrange,c_init[0],c="red",alpha=0.2,label="$c_A^0(t)$")
    cB0, = ax.plot(xrange,c_init[1],c="blue",alpha=0.2,label="$c_A^0(t)$")
    cAt, = ax.plot(xrange,ct[0][frame],c="red",label="$c_A(t)$")
    cBt, = ax.plot(xrange,ct[1][frame],c="blue",label="$c_B(t)$")
    wall1 = ax.axvline(L,ymin=0,c="k",alpha=0.9)
    wall2 = ax.axvline(0,ymin=0,c="k",alpha=0.9)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$C(x,t)$")
    ax.set_title("Diffusion-Reaction system")
    ax.legend(loc="upper right")
    return cA0,cB0,cAt,cBt,wall1,wall2

fig,ax = plt.subplots(figsize=(6,5))
animation = FuncAnimation(fig,GIF,frames=animation_frames,interval=20,blit=True,repeat=True)
animation.save(animation_name,dpi=120,writer=PillowWriter(fps=25))
fig.tight_layout()
plt.show()

# Time evolution of alpha
alpha_t = trange[::int(tp/animation_frames)]
plt.figure(figsize=(6,5))
plt.plot(alpha_t[1:],alpha[1:],c="black")
plt.title('Time evolution of $alpha$')
plt.xlabel('$t$')
plt.ylabel('$alpha$')
plt.show()