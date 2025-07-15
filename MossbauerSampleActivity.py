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
therm_flux = 1e17  # Thermal neutron flux in neutrons/m^2/s
epi_flux = 1e16  # Epithermal neutron flux in neutrons/m^2/s

def dN_dt(t, y):
    N175, N176, N176m, N177, N177m, N178 = y

    # 175Lu
    dN175 = -N175 * (therm_flux * therm_cross["175_to_176"] + epi_flux * epi_cross["175_to_176"]) \
            -N175 * (therm_flux * therm_cross["175_to_176m"] + epi_flux * epi_cross["175_to_176m"])

    # 176Lu
    dN176 = +N175 * (therm_flux * therm_cross["175_to_176"] + epi_flux * epi_cross["175_to_176"]) \
            -N176 * (therm_flux * therm_cross["176_to_177"] + epi_flux * epi_cross["176_to_177"]) \
            -N176 * (therm_flux * therm_cross["176_to_177m"] + epi_flux * epi_cross["176_to_177m"])
    
    # 176mLu
    dN176m = +N175 * (therm_flux * therm_cross["175_to_176m"] + epi_flux * epi_cross["175_to_176m"]) \
             -decay_constants["176mLu"] * N176m

    # 177Lu
    dN177 = +N176 * (therm_flux * therm_cross["176_to_177"] + epi_flux * epi_cross["176_to_177"]) \
            +decay_constants["177mLu"] * N177m * 0.214 \
            +therm_flux * therm_cross["177m_to_177"] * N177m \
            -therm_flux * therm_cross["177_to_178"] * N177 \
            -decay_constants["177Lu"] * N177

    # 177mLu
    dN177m = +N176 * (therm_flux * therm_cross["176_to_177m"] + epi_flux * epi_cross["176_to_177m"]) \
             -decay_constants["177mLu"] * N177m \
             -therm_flux * therm_cross["177m_to_178"] * N177m

    # 178Lu
    dN178 = +therm_flux * therm_cross["177_to_178"] * N177 \
            +therm_flux * therm_cross["177m_to_178"] * N177m \
            -decay_constants["178Lu"] * N178

    return [dN175, dN176, dN176m, dN177, dN177m, dN178]

def main():
    global lu_mass, al_mass, o_mass
    # Define stoichiometric parameters for Lu:Al2O3
    xLu = 0.002  # Desired mole fraction of Lutetium in the final sample
    xAl = 1.998  # Desired mole fraction of Aluminum in the final sample
    xO = 3       # Desired mole fraction of Oxygen in the final sample
    sampleSize = con.pi*(4**2)*0.1 # Volume of Sample in cm^3 (80mm diameter 1mm thick disk)

    # Natural abundances of Lutetium isotopes
    xLu175 = 0.97401  # Abundance of Lu-175
    xLu176 = 0.02599  # Abundance of Lu-176
    xLu176m = 0.0    # Abundance of Lu-176m (metastable state)
    xLu177 = 0.0     # Abundance of Lu-177
    xLu177m = 0.0    # Abundance of Lu-177m (metastable state)
    xLu178 = 0.0     # Abundance of Lu-178

    # For Plotting Later
    colors = {'LuAl2O3': 'blue', 'LuAG': 'green', 'Lu2O3': 'orange'}
    linestyles = {'177mLu': '-', '177Lu': '--'}

    molarMass = { # Molar masses of the materials in g/mol
        "LuAl2O3": lu_mass * xLu + al_mass * xAl + o_mass * xO,
        "LuAG": 3 * lu_mass + 5 * al_mass + 12 * o_mass,
        "Lu2O3": 2 * lu_mass + 3 * o_mass
    }

    density = {
        "LuAl2O3": 3.987,  # Density of Lu-Alumina in g/cm^3
        "LuAG": 6.71,  # Density of LuAG in g/cm^3
        "Lu2O3": 9.42  # Density of Lu2O3 in g/cm^3
    }

    nLu = {
        "LuAl2O3": xLu * (sampleSize * density["LuAl2O3"] / molarMass["LuAl2O3"]) * con.N_A,  # Number of Lutetium atoms in Lu-Alumina
        "LuAG": 3 * (sampleSize * density["LuAG"] / molarMass["LuAG"]) * con.N_A,  # Number of Lutetium atoms in LuAG
        "Lu2O3": 2 * (sampleSize * density["Lu2O3"] / molarMass["Lu2O3"]) * con.N_A  # Number of Lutetium atoms in Lu2O3
    }

    n_isotopes = {
        lattice: {
            "Lu175": xLu175 * nLu[lattice],
            "Lu176": xLu176 * nLu[lattice],
            "Lu176m": xLu176m * nLu[lattice],
            "Lu177": xLu177 * nLu[lattice],
            "Lu177m": xLu177m * nLu[lattice],
            "Lu178": xLu178 * nLu[lattice]
        } for lattice in nLu
    }
    
    # Duration of irradiation and step
    days = 5  # Duration of irradiation in days
    hoursPerDay = 16  # Hours of irradiation per day
    t_span = (0, days * hoursPerDay * 3600)  # Time span for the ODE solver
    t_eval = np.linspace(t_span[0], t_span[1], 1000)  # 1000 steps
    
    #########################################
    ## Transumutation of Lutetium isotopes ##
    #########################################
    results = {}
    for lattice in nLu:
        # Initial conditions for the ODE solver
        N0 = [
            n_isotopes[lattice]["Lu175"],
            n_isotopes[lattice]["Lu176"],
            n_isotopes[lattice]["Lu176m"],
            n_isotopes[lattice]["Lu177"],
            n_isotopes[lattice]["Lu177m"],
            n_isotopes[lattice]["Lu178"]
        ]
        
        # Solve the ODEs
        sol = solve_ivp(dN_dt, t_span, N0, t_eval=t_eval, method="Radau")
        time_days = sol.t / (24 * 3600)  # Convert time to days
        # Store results
        results[lattice] = {
            "solution": sol,
            "time_days": time_days
        }
    
    ######################################################################################
    #### Plotting Activity of Each Isotope Over Time After Irradiation to Max Lu-177m ####
    ######################################################################################
    # Calculate activities in Ci
    pCi = 3.7e10  # 1 Ci = 3.7e10 decays per second
    # Initialize activityMax and activity as nested dictionaries
    activityMax = {lattice: {iso: 0.0 for iso in ['176mLu', '177Lu', '177mLu', '178Lu']} for lattice in nLu}
    activity = {lattice: {iso: None for iso in ['176mLu', '177Lu', '177mLu', '178Lu']} for lattice in nLu}
    decayDays = 60  # Decay time in days after irradiation
    decayTime = np.linspace(0, decayDays*24*3600, 1000)  # 60 Days in seconds
    decayTimeDays = decayTime / (24 * 3600)  # Convert to days
    for lattice in nLu:
        tmpSol = results[lattice]["solution"]
        indmax = np.argmax(tmpSol.y[4])  # Index of max 177mLu
        solnMax = tmpSol.y[2:, indmax]  # Skip the first two isotopes (175Lu, 176Lu)
        for i, iso in enumerate(['176mLu', '177Lu', '177mLu', '178Lu']):
            # Activity in Ci at max time
            activityMax[lattice][iso] = (solnMax[i] * decay_constants[iso]) / pCi
            # print(f'Maximum Activity of ', iso, ' in host lattice ', lattice, ':', activityMax[lattice][iso])
            activity[lattice][iso] = activityMax[lattice][iso] * np.exp(-decay_constants[iso] * decayTime)

    # Plot activity of each isotope over time after irradiation
    plt.figure(2, figsize=(10, 6))
    for lattice in nLu:
        plt.plot(decayTimeDays, activity[lattice]['177mLu'], label=f'177mLu ({lattice})', color=colors[lattice], linestyle=linestyles['177mLu'])
        plt.plot(decayTimeDays, activity[lattice]['177Lu'], label=f'177Lu ({lattice})', color=colors[lattice], linestyle=linestyles['177Lu'])
    plt.yscale('log')
    plt.xlabel('Decay Time (days)')
    plt.ylabel('Activity (Ci)')
    plt.title('Activity of Each Lu Isotope After 5 Days of Neutron Irradiation\n'
    '[Thermal Flux = 1e13 neutrons/cm²/s, Epithermal Flux = 1e12 neutrons/cm²/s]\n'
    '(Lu:Al2O3 - ' + str(np.round(sampleSize * density["LuAl2O3"],1)) + 'g, LuAG - ' + str(np.round(sampleSize * density["LuAG"],1)) + 'g, Lu2O3 - ' + str(np.round(sampleSize * density["Lu2O3"],1)) + 'g)')
    plt.grid(True)
    # plt.ylim(1e-4, 1e3)  # Set y-axis limits for better visibility
    plt.legend()
    plt.tight_layout()
    plt.show()  

    # Print the maximum activity of each isotope
    for lattice in nLu:
        tmpSol = results[lattice]["solution"]
        indmax = np.argmax(tmpSol.y[4])
        tmp177 = tmpSol.y[3]
        tmp177m = tmpSol.y[4]
        print(f"Maximum Activity of 177mLu in {lattice}: {activityMax[lattice]['177mLu']*pCi:.2e} Bq\n Maximum Number of 177mLu in {lattice}: {tmp177m[indmax]:.2e} Atoms")
        print(f"Maximum Activity of 177Lu in {lattice}: {activityMax[lattice]['177Lu']*pCi:.2e} Bq\n Maximum Number of 177Lu in {lattice}: {tmp177[indmax]:.2e} Atoms")

if __name__ == "__main__":
    main()

# Excluded functionality:
#     # Original ODE Solver
#     # Initial conditions
    # N0 = [lu175_atoms[0], lu176_atoms[0], lu177m_atoms[0], lu177_atoms[0], lu176m_atoms[0], lu178_atoms[0]]
    # sol = solve_ivp(dN_dt, t_span, N0, t_eval=t_eval, method="Radau")
    # time_days = sol.t / (24 * 3600)
    ####################################################
    ## Plot number of Atoms of Each Isotope Over Time ##
    ####################################################
    # plt.figure(1, figsize=(10, 6))
    # labels = ['176Lu', '177mLu', '177Lu', '176mLu', '178Lu'] # remove '175Lu' from labels and sol.y for plotting
    # for lattice in nLu:
    #     plotSol = results[lattice]["solution"]
    #     plotSol.y = plotSol.y[1:]  # Skip the first isotope (175Lu)
    #     time_days = results[lattice]["time_days"]
    #     for i, iso in enumerate(labels):
    #         plt.plot(time_days, plotSol.y[i], label=f'{lattice} {iso}', color=colors[lattice], linestyle=linestyles.get(iso, '-'))
    # plt.yscale('log')
    # plt.xlabel('Irradiation Time (days)')
    # plt.ylabel('Number of Atoms (log scale)')
    # plt.title('Neutron Irradiation Various 20mm diameter by 1mm thick Lu Disks\n(Thermal + Epithermal Flux)')
    # plt.legend()
    # plt.grid(True, which="both", ls="-")
    # plt.tight_layout()
    # plt.show()
#     ###########################################################
#     ## Calculate Transmutation for Various Lu-176 Abundances ##
#     ###########################################################
#     # Varying Lu-176 abundance from 2.599% to 50%
#     lu176Abundances = np.linspace(0.02599, 0.5, num=100)  # Varying Lu-176 abundance from 2.599% to 50%
#     lis177mLu = []  # List to store the maximum number of Lu-177m atoms for each abundance
#     lis177Lu = []  # List to store the maximum number of Lu-177 atoms for each abundance
#     for abund in lu176Abundances:
#         x176Lu = abund
#         x175Lu = 1 - x176Lu

#         # Update initial conditions
#         init_175Lu = [x175Lu * nLu]
#         init_176Lu = [x176Lu * nLu]
#         init_176mLu = [0]
#         init_177Lu = [0]
#         init_177mLu = [0]
#         init_178Lu = [0]

#         # Solve the ODE
#         tmpN0 = [init_175Lu[0], init_176Lu[0], init_177mLu[0], init_177Lu[0], init_176mLu[0], init_178Lu[0]]
#         tmpSol = solve_ivp(dN_dt, t_span, tmpN0, t_eval=t_eval, method="Radau")
        
#         lis177mLu.append(np.max(tmpSol.y[2]))  # Maximum number of Lu-177m atoms
#         lis177Lu.append(np.max(tmpSol.y[3]))  # Maximum number of Lu-177 atoms

    # #######################################################################
    # ## Plot Activity of Lu-177m and Lu-177 for Varying Lu-176 Abundances ##
    # #######################################################################
    # act177mLu = np.array(lis177mLu) * decay_constants['177mLu'] / pCi  # Activity of Lu-177m in Ci
    # act177Lu = np.array(lis177Lu) * decay_constants['177Lu'] / pCi  # Activity of Lu-177 in Ci

    # plt.figure(3, figsize=(10, 6))
    # plt.yscale('log')  # Set y-axis to logarithmic scale
    # plt.plot(lu176Abundances * 100, act177mLu, label='177mLu')
    # plt.plot(lu176Abundances * 100, act177Lu, label='177Lu')
    # plt.xlabel('Lu-176 Abundance (%)')
    # plt.ylabel('Activity (Ci)')
    # plt.title('Activity of Lu-177m and Lu-177 for Varying Lu-176 Abundances\n (5 days of irradiation)')
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.show()