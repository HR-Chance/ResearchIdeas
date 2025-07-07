import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as con

# Define the molecular weights (g/mol) as global variables
lu_mass = 174.9668      # Lutetium
al_mass = 26.981538      # Aluminum
o_mass = 15.999          # Oxygen

def calculate_activity():
    pass

def plot_results():
    pass

def main():
    global lu_mass, al_mass, o_mass
    xLu = 0.002  # Desired mole fraction of Lutetium in the final sample
    xAl = 1.998  # Desired mole fraction of Aluminum in the final sample
    xO = 3       # Desired mole fraction of Oxygen in the final sample
    luAlO_mass = lu_mass * xLu + al_mass * xAl + o_mass * xO  # Molar mass of Lu-Alumina in g/mol

    lu175_cross = 23.3*1e-28  # Neutron capture cross-section for Lu-175 in m^2
    lu176_cross = 2057*1e-28  # Neutron capture cross-section for Lu-176 in m^2
    thermal_neutron_flux = 1e17  # Thermal neutron flux in neutrons/m^2/s
    epithermal_neutron_flux = 1e16  # Epithermal neutron flux in neutrons/m^2/s
    
    ## Initial Quantities
    # Mass of sample and number of Lu Atoms
    massLuAlO = 25  # Mass of Lu-Alumina sample in grams
    nLu = lu_mass * xLu * (massLuAlO/luAlO_mass) * con.N_A # Number of Lutetium atoms
    
    # Isotopic Abundances
    xLu175 = 0.97401  # Abundance of Lu-175
    xLu176 = 0.02599  # Abundance of Lu-176
    
    # Initialize quantities
    lu175_atoms = [xLu175 * nLu]   # Lu-175
    lu176_atoms = [xLu176 * nLu]   # Lu-176
    lu176m_atoms = [0]             # Lu-176m
    lu177_atoms = [0]              # Lu-177
    lu177m_atoms = [0]             # Lu-177m
    lu178_atoms = [0]              # Lu-178

    # Duration of irradiation and step
    days = 60  # Duration of irradiation in days
    time = np.linspace(0, days*86400, num=1000)  # Time in seconds
    dt = time[1] - time[0]  # Time step in seconds
    print(f"Time step: {dt} seconds")
    # Calculate the number of atoms over time
    # for t in time:
    #     lu175_atoms.append(lu175_atoms[-1] - lu175_cross * thermal_neutron_flux * lu175_atoms[-1] * t)


if __name__ == "__main__":
    main()