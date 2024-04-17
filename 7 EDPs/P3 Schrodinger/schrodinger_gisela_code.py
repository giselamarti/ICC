import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
from scipy.sparse import diags
from matplotlib.animation import FuncAnimation, PillowWriter

# Constants
hbar = 1.0545718e-34           # Reduced Planck's constant in J*s
me = 9.1093837015e-31          # mass of an electron 
a0 = 5.29177210903e-11         # bohr radius
amu_to_me = 1836

# System parameters
mH = 1.0                       # amu
mH2 = 2.0                      # amu
Eh = hbar**2/(me*a0**2)        # Hartree energy
mu = mH*mH2 / (mH+mH2) * amu_to_me # Reduced mass (me)

# Initial conditions for Gaussian wave packet
x0 = 5.5                # au
delta_x_squared = 0.04  # au^2
px0_over_hbar = 5.5     # au^-1

# Grid parameters
Lx = 10.0               # length
Nx = 100                # Grid points
dx = Lx/Nx              # Subintervals
xrange = np.arange(0,Lx+dx,dx)

dt = 1.979e-16 *Eh/hbar # Time step (s)
Nt = 1000               # Time points
t = Nt*dt                # Total time
trange = np.arange(0,t+dt,dt)

# Animation parameters
animation_frames = 100
animation_name = "Schroedinger.gif"

# Create tridiagonal matrix with scipy
dim = len(xrange)
a = dt/(4*mu*dx**2)     # Parameter a 
off_diag = a
matrix_A = diags([-a,2j+2*a,-a], [-1, 0, 1], shape=(dim, dim), format='csr')
matrix_B = diags([a,2j-2*a,a], [-1, 0, 1], shape=(dim, dim), format='csr')

WF = np.array(1/(2*np.pi*delta_x_squared)**0.25 * np.exp(1*px0_over_hbar*(xrange-x0))* \
               np.exp(-(xrange-x0)**2/(4*delta_x_squared)))
WF_init = WF.copy()

WFt = []
for i,ti in enumerate(trange):
    WF[-1] = WF[0]
    if i%int(Nt/animation_frames) == 0 : 
        WFt.append(WF.copy())
    WF = inv(matrix_A.toarray())@matrix_B.toarray()@WF
WFt = np.array(WFt)

# Start animation
def GIF(frame):
    """Function that creates a frame for the GIF."""
    ax.clear()
    WF0_frame, = ax.plot(xrange,WF_init,c="blue",alpha=0.2,label="WF inicial")
    WFt_frame, = ax.plot(xrange,WFt[frame],c="red",label="WF final")
    wall1 = ax.axvline(Lx,ymin=0,c="k",alpha=0.9)
    wall2 = ax.axvline(0,ymin=0,c="k",alpha=0.9)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$\Psi$")
    ax.set_ylim([-1.5,5])
    ax.set_title("Wavefunction evolution")
    ax.legend(loc="upper right")
    return WF0_frame,WFt_frame,wall1,wall2

fig,ax = plt.subplots(figsize=(6,5))
animation = FuncAnimation(fig,GIF,frames=animation_frames,interval=20,blit=True,repeat=True)
animation.save(animation_name,dpi=120,writer=PillowWriter(fps=25))
fig.tight_layout()
plt.show()