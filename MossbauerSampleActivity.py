import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as con

# Code written by Henry R. Chance
# This code is intended to calculate the activity of a Mossbauer sample containing Lutetium isot

# Define the molecular weights (g/mol) as global variables
lu_mass = 174.9668      # Lutetium
al_mass = 26.981538      # Aluminum
o_mass = 15.999          # Oxygen

# Neutron capture cross-sections in m^2
s175_176 = 6.6*1e-28  # Neutron capture cross-section for Lu-175 in m^2
i175_176 = 620*1e-28  # Fast neutron capture cross-section for Lu-175 in m^2
s175_176m = 16.7*1e-28  # Neutron capture cross-section to Lu-176m in m^2
i175_176m = 550*1e-28  # Fast neutron capture cross-section to Lu-176m in m^2
s176_177m = 2.8*1e-28  # Neutron capture cross-section for Lu-176m in m^2
i176_177m = 4.7*1e-28  # Fast neutron capture cross-section for Lu-176m in m^2
s176_177 = 2020*1e-28  # Neutron capture cross-section for Lu-176 in m^2
i176_177 = 1087*1e-28  # Fast neutron capture cross-section for Lu-176 in m^2
s177_178 = 1000*1e-28  # Neutron capture cross-section for Lu-177 in m^2
s177m_178 = 417*1e-28  # Neutron capture cross-section for Lu-177m in m^2
s177m_177 = 209*1e-28  # Neutron capture cross-section for Lu-177m in m^2

# Decay constants
d176m = np.log(2)/(3.664*3600)  # decay constant of Lu-176m in s^-1
d177 = np.log(2)/(6.6443*24*3600)  # decay constant of Lu-177 in s^-1
d177m = np.log(2)/(160.4*24*3600)  # decay constant of Lu-177m in s^-1
d178 = np.log(2)/(28.4*60)  # decay constant of Lu-178 in s^-1


# Reactor Flux in neutrons per m^2 per second
therm_flux = 1e20  # Thermal neutron flux in neutrons/m^2/s
epi_flux = 1e16  # Epithermal neutron flux in neutrons/m^2/s

def dN175(n175):
    global s175_176, i175_176, s175_176m, i175_176m, therm_flux, epi_flux
    change = -n175*((s175_176+s175_176m)*therm_flux + (i175_176+i175_176m)*epi_flux)
    # print(f"dN175: {change}")
    return change

def dN176(n175, n176):
    global s175_176, s176_177, s176_177m, i175_176, i176_177, i176_177m, therm_flux, epi_flux
    change = n175*(s175_176*therm_flux + i175_176*epi_flux) - n176*((s176_177+s176_177m)*therm_flux + (i176_177+i176_177m)*epi_flux)
    return change

def dN176m(n175, n176m):
    global s175_176m, i175_176m, d176m, therm_flux, epi_flux
    change = n175*(s175_176m*therm_flux + i175_176m*epi_flux) - n176m*d176m
    return change

def dN177(n176, n177, n177m):
    global s176_177, s177_178, s177m_178, i176_177, s177m_177, d177, d177m, therm_flux, epi_flux
    change = n176*(s176_177*therm_flux + i176_177*epi_flux) - n177*(s177_178*therm_flux + d177) + n177m*(s177m_178*therm_flux + 0.217*d177m)
    # Note: The 0.217 factor is derived from the ratio of decays toward Lu-177 from Lu-177m
    return change

def dN177m(n176, n177m):
    global s176_177m, i176_177m, s177m_177, s177m_178, d177m, therm_flux, epi_flux
    change = n176*(s176_177m*therm_flux + i176_177m*epi_flux) - n177m*((s177m_178+s177m_177)*therm_flux + d177m)
    return change

def dN178(n177, n177m, n178):
    global s177_178, s177m_178, d178, therm_flux, epi_flux
    change = therm_flux*(n177*s177_178 + n177m*s177m_178) - n178*d178
    return change

def main():
    global lu_mass, al_mass, o_mass
    xLu = 0.002  # Desired mole fraction of Lutetium in the final sample
    xAl = 1.998  # Desired mole fraction of Aluminum in the final sample
    xO = 3       # Desired mole fraction of Oxygen in the final sample
    luAlO_mass = lu_mass * xLu + al_mass * xAl + o_mass * xO  # Molar mass of Lu-Alumina in g/mol
    
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
    time = np.linspace(0, days*86400, num=1000000)  # Time in seconds
    dt = time[1] - time[0]  # Time step in seconds
    print(f"Time step: {dt} seconds")
    # Calculate the number of atoms over time
    for t in time:
        lu175_atoms.append(lu175_atoms[-1]+dN175(lu175_atoms[-1]))
        lu176_atoms.append(lu176_atoms[-1]+dN176(lu175_atoms[-1], lu176_atoms[-1]))
        lu176m_atoms.append(lu176m_atoms[-1]+dN176m(lu175_atoms[-1], lu176m_atoms[-1]))
        lu177_atoms.append(lu177_atoms[-1]+dN177(lu176_atoms[-1], lu177_atoms[-1], lu177m_atoms[-1]))
        lu177m_atoms.append(lu177m_atoms[-1]+dN177m(lu176_atoms[-1], lu177m_atoms[-1]))
        lu178_atoms.append(lu178_atoms[-1]+dN178(lu177_atoms[-1], lu177m_atoms[-1], lu178_atoms[-1]))

    
    plt.figure(figsize=(10, 6))
    plt.plot(time, lu175_atoms[:-1], label='Lu-175')
    plt.plot(time, lu176_atoms[:-1], label='Lu-176')
    plt.plot(time, lu176m_atoms[:-1], label='Lu-176m')
    plt.plot(time, lu177_atoms[:-1], label='Lu-177')
    plt.plot(time, lu177m_atoms[:-1], label='Lu-177m')
    plt.plot(time, lu178_atoms[:-1], label='Lu-178')
    plt.xlabel('Time (days)')
    plt.ylabel('Number of Atoms')
    plt.title('Number of Atoms of Each Isotope Over Time')
    plt.legend()
    plt.tight_layout()
    plt.yscale('log')
    plt.show()
    

if __name__ == "__main__":
    main()