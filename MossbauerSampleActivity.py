import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as con
from scipy.integrate import solve_ivp

# Code written by Ahmed Helal, modified by Henry R. Chance
# Date: 05/10/2025  v3
# This code is intended to calculate the activity of a Mossbauer sample containing Lutetium isotopes

# Define the molecular weights (g/mol) as global variables
lu_mass = 174.9668      # Lutetium
al_mass = 26.981538      # Aluminum
o_mass = 15.999          # Oxygen

therm_cross = {
    "175_to_176": 6.6e-28,    # Neutron capture cross-section for Lu-175
    "175_to_176m": 16.7e-28,  # Neutron capture cross-section to Lu-176m
    "176_to_177m": 2.8e-28,  # Neutron capture cross-section for Lu-176m
    "176_to_177": 2020e-28,  # Neutron capture cross-section for Lu-176
    "177_to_178": 1000e-28,  # Neutron capture cross-section for Lu-177
    "177m_to_178": 417e-28,  # Neutron capture cross-section for Lu-177m
    "177m_to_177": 209e-28,  # Inelastic scattering cross-section for Lu-177m
}

epi_cross = {
    "175_to_176": 620e-28,  # Fast neutron capture cross-section for Lu-175
    "175_to_176m": 550e-28,  # Fast neutron capture cross-section to Lu-176m
    "176_to_177m": 4.7e-28,  # Fast neutron capture cross-section for Lu-176m
    "176_to_177": 1087e-28,  # Fast neutron capture cross-section for Lu-176
}

decay_constants = {
    "177mLu": np.log(2)/(160.4*24*3600),  # decay constant of Lu-177m in s^-1
    "177Lu": np.log(2)/(6.6443*24*3600),  # decay constant of Lu-177 in s^-1
    "176mLu": np.log(2)/(3.664*3600),  # decay constant of Lu-176m in s^-1
    "178Lu": np.log(2)/(28.4*60)  # decay constant of Lu-178 in s^-1
}

# Reactor Flux in neutrons per m^2 per second
therm_flux = 1e20  # Thermal neutron flux in neutrons/m^2/s
epi_flux = 1e16  # Epithermal neutron flux in neutrons/m^2/s

def dN_dt(t, y):
    N175, N176, N177m, N177, N176m, N178 = y

    # 175Lu
    dN175 = -N175 * (therm_flux * therm_cross["175_to_176"] + epi_flux * epi_cross["175_to_176"]) \
            -N175 * (therm_flux * therm_cross["175_to_176m"] + epi_flux * epi_cross["175_to_176m"])

    # 176Lu
    dN176 = +N175 * (therm_flux * therm_cross["175_to_176"] + epi_flux * epi_cross["175_to_176"]) \
            -N176 * (therm_flux * therm_cross["176_to_177"] + epi_flux * epi_cross["176_to_177"]) \
            -N176 * (therm_flux * therm_cross["176_to_177m"] + epi_flux * epi_cross["176_to_177m"])

    # 177mLu
    dN177m = +N176 * (therm_flux * therm_cross["176_to_177m"] + epi_flux * epi_cross["176_to_177m"]) \
             -decay_constants["177mLu"] * N177m \
             -therm_flux * therm_cross["177m_to_178"] * N177m

    # 177Lu
    dN177 = +N176 * (therm_flux * therm_cross["176_to_177"] + epi_flux * epi_cross["176_to_177"]) \
            +decay_constants["177mLu"] * N177m * 0.214 \
            +therm_flux * therm_cross["177m_to_177"] * N177m \
            -therm_flux * therm_cross["177_to_178"] * N177 \
            -decay_constants["177Lu"] * N177

    # 176mLu
    dN176m = +N175 * (therm_flux * therm_cross["175_to_176m"] + epi_flux * epi_cross["175_to_176m"]) \
             -decay_constants["176mLu"] * N176m

    # 178Lu
    dN178 = +therm_flux * therm_cross["177_to_178"] * N177 \
            +therm_flux * therm_cross["177m_to_178"] * N177m \
            -decay_constants["178Lu"] * N178

    return [dN175, dN176, dN177m, dN177, dN176m, dN178]

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
    days = 30  # Duration of irradiation in days
    t_span = (0, days * 86400)  # Time span for the ODE solver
    t_eval = np.linspace(t_span[0], t_span[1], 1000)  # 1000 steps

    # Initial conditions
    N0 = [lu175_atoms[0], lu176_atoms[0], lu177m_atoms[0], lu177_atoms[0], lu176m_atoms[0], lu178_atoms[0]]
    sol = solve_ivp(dN_dt, t_span, N0, t_eval=t_eval, method="Radau")
    time_days = sol.t / (24 * 3600)

    ####################################
    # Plot number of Atoms of Each Isotope Over Time
    #####################################
    plt.figure(figsize=(10, 6))
    # remove '175Lu' from labels and sol.y for plotting
    labels = ['176Lu', '177mLu', '177Lu', '176mLu', '178Lu']
    soln = sol.y[1:]  # Skip the first isotope (175Lu)
    indmax = np.argmax(soln[1])  # Find the index of the maximum value in 177mLu
    for i in range(5):
        plt.plot(time_days, soln[i], label=labels[i])

    plt.axvline(x=time_days[indmax], color='r', linestyle='--', label='Max 177mLu at ' + str(np.round(time_days[indmax]*24,decimals=0)) + ' hours')
    plt.yscale('log')
    plt.xlabel('Irradiation Time (days)')
    plt.ylabel('Number of Atoms (log scale)')
    plt.title('Neutron Irradiation of ' + str(massLuAlO) + 'g Lu:Al2O3\n(Thermal + Epithermal Flux)')
    plt.legend()
    plt.grid(True, which="both", ls="-")
    plt.tight_layout()
    plt.show()



    ######################################################################################
    #### Plotting Activity of Each Isotope Over Time After Irradiation to Max Lu-177m ####
    ######################################################################################
    # Calculate activities in Ci
    pCi = 3.7e10  # 1 Ci = 3.7e10 decays per second
    activity = {}
    activityMax = {}
    solnMax = sol.y[2:, indmax]  # Skip the first two isotopes (175Lu)
    decayTime = np.linspace(0, 60*24*3600, 1000)  # 5 Days of Decay time in seconds
    decayTimeDays = decayTime / (24 * 3600)  # Convert to days
    for i, iso in enumerate(['177mLu', '177Lu', '176mLu', '178Lu']):
        # Activity in Ci at max time
        activityMax[iso] = (solnMax[i] * decay_constants[iso]) / pCi  
        
    for i, iso in enumerate(['177mLu', '177Lu', '176mLu', '178Lu']):
        # Calculate activity over decay time
        activity[iso] = activityMax[iso] * np.exp(-decay_constants[iso] * decayTime)

    # Plot mass of each isotope over time
    plt.figure(figsize=(10, 6))

    for i, iso in enumerate(['177mLu', '177Lu', '176mLu', '178Lu']):
        plt.plot(decayTimeDays, activity[iso], label=f'{iso}')
    plt.yscale('log')
    plt.xlabel('Decay Time (days)')
    plt.ylabel('Activity (Ci)')
    plt.title('Activity of Each Lu Isotope vs Time')
    plt.grid(True)
    plt.ylim(0.1, 1e6)  # Set y-axis limits for better visibility
    plt.legend()
    plt.tight_layout()
    plt.show()   
    

if __name__ == "__main__":
    main()