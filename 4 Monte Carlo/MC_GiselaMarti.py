"""
Monte Carlo program for two water molecules.

Gisela Martí, ICC 2023
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#location of pdb files
path = "C:/Users/gisela/OneDrive/Documentos/Master/ICC2/MC/wat_gas_mol_origin.pdb" 

"""
PART 1: DEFINING FUNCTIONS AND PARAMETERS
"""

def read_pdb (pdb):
    """
    Reads atomic coordinates from a PDB file.
    
    Parameters:
    - pdb: Path to the PDB file.
    
    Returns:
    - coord: Numpy array containing atomic coordinates. Shape: (num_atoms, 3)
    """
    # Determine number of atoms
    file = open(pdb)
    num_atoms = 0
    for line in file:
        if 'ATOM' in line:
            num_atoms += 1
    # Read coordinates           
    coord = np.zeros((num_atoms,3))
    file = open(pdb)
    index = 0
    for line in file:
        if 'ATOM' in line:
            line = line.split() 
            coord[index,0] = float(line[5])
            coord[index,1] = float(line[6])
            coord[index,2] = float(line[7])
            index +=1
    return coord

class water():
    
    def __init__(self):
        """
        Initializes the water class. Reads atomic coordinates from a PDB file.
        """
        self.coord = read_pdb(path)
        
    def translation(self,Tx,Ty,Tz):
        """
        Translates the molecular coordinates by the specified translation vector (Tx, Ty, Tz).
        
        Parameters:
        - Tx: Translation along the x-axis.
        - Ty: Translation along the y-axis.
        - Tz: Translation along the z-axis.
        """
        self.coord = self.coord.T
        for i,T in enumerate([Tx,Ty,Tz]):
            self.coord[i] += T
        self.coord = self.coord.T
        
    def origin(self):
        """
        Moves the molecule to the origin of coordinates (centered on the Oxygen atom).
        """
        c = self.coord[0] 
        self.translation(-c[0],-c[1],-c[2])

    def rotation(self,Ax,Ay,Az,angle_units):
        """
        Rotates the molecule around the x, y, and z axes by the specified angles.
        
        Parameters:
        - Ax: Angle of rotation around the x-axis.
        - Ay: Angle of rotation around the y-axis.
        - Az: Angle of rotation around the z-axis.
        - angle_units: "deg" for degrees or "rad" for radians.
        """
        angle_units = angle_units.lower()
        if angle_units == "deg":
            Ax = np.radians(Ax)
            Ay = np.radians(Ay)
            Az = np.radians(Az)
        elif angle_units == "rad":
            pass
        
        # rotation
        Rx = np.array([[1,0,0],[0,np.cos(Ax),-np.sin(Ax)],[0,np.sin(Ax),np.cos(Ax)]])
        Ry = np.array([[np.cos(Ay),0,np.sin(Ay)],[0,1,0],[-np.sin(Ay),0,np.cos(Ay)]])
        Rz = np.array([[np.cos(Az),-np.sin(Az),0],[np.sin(Az),np.cos(Az),0],[0,0,1]])

        # moves to origin
        temp = [self.coord[0][i] for i in range(3)]
        self.origin()

        # rotates the molecule
        for i in range(3):
            for R in [Rx,Ry,Rz]:
                ctemp = self.coord[i].copy()
                self.coord[i] = np.dot(R,ctemp) 

        # Returns the rotated molecule to the same place
        self.translation(*temp) 

    def calculate_energies(self,other):
        """
        Calculates interaction energies (Van der Waals and Coulomb) between two molecules.
        
        Parameters:
        - other: Another instance of the water class (second water molecule)
        
        Returns:
        - Total energy, Coulombic energy, and Van der Waals energy.
        """ 
        coord1,coord2 = self.coord,other.coord
        d_ijs = coord1-coord2[:,None]
        d_ij = np.sqrt((d_ijs**2).sum(axis=-1)).T
        
        vdw_energy = (A/d_ij**12 - B/d_ij**6).sum()     
        coulomb_energy = (332.0*Q/d_ij).sum()
        E_total = vdw_energy + coulomb_energy
        self.E_total = E_total
        return E_total,coulomb_energy,vdw_energy

    def reset(self):
        """
        Resets the total energy (E_total) to zero.
        """
        self.E_total = 0

def dE_MC_step(dE):
    """
    Accepts or rejects a Monte Carlo step based on the change in energy (dE).
    
    Parameters:
    - dE: Change in energy for the proposed step.
    
    Returns:
    - boolean value: True if the step is accepted, False otherwise.
    """
    if dE < 0: 
        return True
    else:
        rand = np.random.random()
        if np.exp(-beta*dE) > rand:
            return True
        else:
            return False

# Simulation parameters
T = 300.                # Temperature
kb = 0.00119872041      # Boltzman constant (kcal/(mol⋅K))
beta = 1./(kb*T)        # Beta factor
coord = read_pdb(path)  # coordinates of the first molecule

# Energy parameters
A = np.array([[5.81935564e+05, 3.28388768e+02, 3.28388768e+02],
              [3.28388768e+02, 3.12143700e-06, 3.28388768e+02],
              [3.28388768e+02, 3.28388768e+02, 3.12143700e-06]])
B = np.array([[5.94825035e+02, 1.04803190e+01, 1.04803190e+01],
              [1.04803190e+01, 7.58000000e-04, 1.04803190e+01],
              [1.04803190e+01, 1.04803190e+01, 7.58000000e-04]])

# Coulomb charges
q = np.array([-0.834,+0.417,+0.417])
Q = np.dot(q[:,None],q[:,None].T)

"""
PART 2: MONTE CARLO SAMPLING
"""

# Parameters
n = 2       # Number of waters (1 fixed and 1 moving)
N = 100000  # Number of iterations
L = 5       # Half the box length (angstroms)

# Initialization of fixed water molecule at center
fixed_water = water()

# Lists to save the data
waters = [water() for _ in range(n-1)]
energies = [[] for _ in range(n-1)]
coords = [[] for _ in range(n-1)]

# Sampling loop
for j in range(N):
    for i,w in enumerate(waters):
        
        # start from random values for the positions and angles
        random_T = np.random.uniform(-L,L, size=3)
        random_A = np.random.uniform(-180,180,size=3) # degrees
        
        w.translation(*random_T)
        w.rotation(*random_A,"deg")
        E_total,coulomb_energy,vdw_energy = w.calculate_energies(fixed_water)
        energies[i].append(E_total)
        coords[i].append(w.coord.copy())

        w.reset()
        w.origin()

energies = np.array(energies)
minimum_energy_sampling = np.min(energies,axis=-1)        # minimum E found in the sampling
minimum_energy_i_sampling = np.argmin(energies,axis=-1)   # position in the array
minimum_energy = np.min(minimum_energy_sampling)         
minimum_energy_i = np.argmin(minimum_energy_sampling)
coords = np.array(coords)
minimum_coords = []
for i in range(n-1):
    w_min = water() # store information about the water molecule configuration that
                    # corresponds to the minimum energy found during the simulation.
    c_min = coords[i][minimum_energy_i_sampling[i]]   # coordinates of the water molecule with minimum 
                                                      # energy in the current iteration of the loop
    w_min.coord = c_min
    minimum_coords.append(c_min)
    
w_min = water()
w_min.coord = minimum_coords[minimum_energy_i]
    
print("MONTE CARLO SAMPLING")    
print(" ")
print("Minimum energy found in the sampling: ", *minimum_energy_sampling, "(kcal/mol)")
print("--------------------------------------------------------")

"""
PART 3: MONTE CARLO METROPOLIS ALGORITHM
"""

# Parameters
N_MC = 1000000     # number of iterations
step = 0.1         # translation step size
angle = 10         # angle step size
acceptance = [0,0] # acceptance ([accepted, rejected])
acc_steps = [0]    # total accepted steps
w = w_min

# Lists to save the data
energy_MC = [minimum_energy]
coulomb_energy_MC = []
vdw_energy_MC = []
coords_MC = [w_min.coord.copy()]

# Metropolis algorithm loop
for i in range(N_MC):

    random_T = np.random.uniform(-step,step,size=3)
    random_A = np.random.uniform(-angle,angle,size=3)
    w.translation(*random_T)
    w.rotation(*random_A,"deg")
    E_total,coulomb_energy,vdw_energy = w.calculate_energies(fixed_water)   
    
    dE = E_total - energy_MC[-1] 
    if not dE_MC_step(dE):
        acceptance[1] += 1              # Rejected step
        w.coord = coords_MC[-1].copy()  # Return to the last configuration
        continue                        # continue to next step
    else: 
        acceptance[0] += 1              # Accepted step

	# save the energies and coordinates of only the accepted steps
        energy_MC.append(E_total)           
        coulomb_energy_MC.append(coulomb_energy) 
        vdw_energy_MC.append(vdw_energy)
        coords_MC.append(w.coord.copy())
        acc_steps.append(i+1)

    w.reset()

print("MONTE CARLO METROPOLIS ALGORITHM RESULTS")  
print(" ")

E_min = np.min(energy_MC)       # minimum energy found
E_min_i = np.argmin(energy_MC)  # index of the minimum energy
coord_min = coords_MC[E_min_i]

print("Minimum energy found:",E_min,"(kcal/mol)")
print(" ")
print("Contributions of the minimum energy:")
print("Coulomb Energy:",coulomb_energy_MC[E_min_i-1],"(kcal/mol)")
print("Van der Waals energy:",vdw_energy_MC[E_min_i-1],"(kcal/mol)")
print(" ")
print("Acceptance percentage:",100*acceptance[0]/N)
print("Number of accepted steps:",acceptance[0])
print("Number of rejected steps:",acceptance[1])

# Downsampling for better visibility
downsample_factor = 100
downsampled_steps = acc_steps[::downsample_factor]
downsampled_energy = energy_MC[::downsample_factor]

# Plot of the energy
plt.figure(figsize=(12, 6))
plt.plot(downsampled_steps, downsampled_energy,c="black")
plt.title('Total energy variation during the simulation (downsampled)')
plt.xlabel('MC Steps')
plt.ylabel('Energy (kcal/mol)')
plt.xlim([0,N_MC])
plt.grid(True)
plt.show()

# Plotting the configurations of the water molecules in a 3D box

def plot_configuration(ax, coords, color, label):
    oxygen = coords[::3]
    hydrogen1 = coords[1::3]
    hydrogen2 = coords[2::3]
    ax.scatter(oxygen[:, 0], oxygen[:, 1], oxygen[:, 2], color=color[0], s=100, label=label[0])
    ax.scatter(hydrogen1[:, 0], hydrogen1[:, 1], hydrogen1[:, 2], color=color[1], s=50, label=label[1])
    ax.scatter(hydrogen2[:, 0], hydrogen2[:, 1], hydrogen2[:, 2], color=color[1], s=50)

# Plot initial configuration
fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(121, projection='3d')
ax1.set_title('Initial Configuration')
ax1.set_xlim([-L, L])
ax1.set_ylim([-L, L])
ax1.set_zlim([-L, L])
plot_configuration(ax1, fixed_water.coord, ['grey', 'black'], ['Fixed Oxygen', 'Fixed Hydrogen'])
for w in waters:
    plot_configuration(ax1, w.coord, ['red', 'blue'], ['Oxygen', 'Hydrogen'])
ax1.legend(loc='upper right', fontsize='small')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')

# Plot final configuration
ax2 = fig.add_subplot(122, projection='3d')
ax2.set_title('Final Configuration')
ax2.set_xlim([-L, L])
ax2.set_ylim([-L, L])
ax2.set_zlim([-L, L])
plot_configuration(ax2, fixed_water.coord, ['grey', 'black'], ['Fixed Oxygen', 'Fixed Hydrogen'])
plot_configuration(ax2, w_min.coord, ['red', 'blue'], ['Oxygen', 'Hydrogen'])
ax2.legend(loc='upper right', fontsize='small')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
plt.show()

# Save the final coordinates to a PDB file
output_pdb_path = "output_final_coordinates.pdb"

def write_pdb(file_path, coords):
    """
    Write atomic coordinates to a PDB file.

    Parameters:
    - file_path: Path to the PDB file.
    - coords: Numpy array containing atomic coordinates. Shape: (num_atoms, 3)
    """
    with open(file_path, 'w') as f:
        f.write("REMARK\n")
        for i, coord in enumerate(coords):
            atom_type = 'OW' if i < 3 else 'HW'  # Assuming the first three lines correspond to the fixed water molecule
            f.write(f"ATOM      {i + 1}  {atom_type}  WAT     0      {coord[0]:8.3f}  {coord[1]:8.3f}  {coord[2]:8.3f}  0.00  0.00           {atom_type[0]}\n")
        f.write("TER\n")

write_pdb(output_pdb_path, np.vstack((fixed_water.coord, w_min.coord)))

print("PDB file written:", output_pdb_path)