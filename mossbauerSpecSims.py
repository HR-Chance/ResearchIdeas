import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as con

# Code written by Henry R. Chance
# This code is intended to simulate Mossbauer Spectroscopy Emission/Absorption Signals

def recoilEnergy(mass, energy):
    # Source: Brent  Fultz, "Mössbauer Spectrometry", in  Characterization  of  Materials . Elton  Kaufmann, Editor (John  Wiley, New  York, 2011) .
    """
    Calculate the recoil energy of a nucleus after emitting or absorbing a gamma photon.
    Parameters:
    mass (float): Mass of the nucleus in kg.
    energy (float): Energy of the gamma photon in eV.
    Returns:
    float: Recoil energy in eV.
    """
    return ((energy*con.e)**2 / (2 * mass * con.c**2))/con.e  # Convert to eV

def recoillessFraction(mass, energy, debeye, temperature):
    # Source: Eq. 19 from Brent Fultz, "Mössbauer Spectrometry", in Characterization of Materials. Elton Kaufmann, Editor (John Wiley, New York, 2011).
    # Eq. 2.8 from Kiejna, A., and K. F. Wojciechowski. "The surface of real metals." Met. Surf. Electron Phys (1996): 19-32.
    """
    Calculate the recoilless fraction for a given nucleus.
    Parameters:
    mass (float): Mass of the nucleus in kg.
    energy (float): Energy of the gamma photon in eV.
    debeye (float): Debye temperature in K.
    temperature (float): Temperature in K.
    Returns:
    float: Recoilless fraction.
    """
    return np.exp(-(((energy*con.e)/(con.hbar*con.c))**2) * ((3*temperature*con.hbar**2)/(mass*con.k*debeye**2)))  # Recoilless fraction formula

def maximumAbsorptionCross(Ia, Ib, alpha, waveLength, F):
    # Source: Eq. 10, 11, 15 from Hans Frauenfelder, "The Mossbauer Effect"
    """
    Calculate the maximum absorption cross-section for a given nucleus.
    Parameters:
    Ia (float): Spin of the ground state.
    Ib (float): Spin of the excited state.
    alpha (float): Internal conversion fraction.
    wavelength (float): Wavelength of the gamma photon in m.
    F (float): Recoilless fraction.
    Returns:
    float: Maximum absorption cross-section in m^2.
    """
    return 2*con.pi * (waveLength**2) * ((2*Ib+1)/(2*Ia+1)) * (1/(2*(1+alpha))) * F  # Maximum absorption cross-section formula

def absorptionCross(detune, linewidth, sigMax):
    # Source: Eq. 10 from Hans Frauenfelder, "The Mossbauer Effect"
    """
    Calculate the absorption cross-section for a given detuning and linewidth.
    Parameters:
    detune (float): Detuning energy in eV.
    linewidth (float): Linewidth in eV.
    sigMax (float): Maximum absorption cross-section in m^2.
    Returns:
    float: Absorption cross-section in m^2.
    """
    return sigMax * (linewidth/2)**2 / ((linewidth/2)**2 + detune**2)  # Absorption cross-section formula

def velocityDetuning(vel,freq):
    """
    Calculate the detuning energy based on the velocity of the nucleus.
    Parameters:
    vel (float): Velocity in m/s.
    freq (float): Frequency of transition in Hz.
    Returns:
    float: Detuning energy in eV.
    """
    return 2*con.pi*(con.hbar/con.e)*(vel/con.c)*freq  # Detuning energy formula

def main():
    # Define parameters
    lu_mass = 174.9668      # Lutetium
    al_mass = 26.981538      # Aluminum
    o_mass = 15.999          # Oxygen
    xLu = 0.002  # Desired mole fraction of Lutetium in the final sample
    xAl = 1.998  # Desired mole fraction of Aluminum in the final sample
    xO = 3       # Desired mole fraction of Oxygen in the final sample
    mLu = 176.9437636 * con.atomic_mass # Lutetium mass in kg
    Egamma = 121.62e3  # Gamma energy in eV
    freqGamma = Egamma * con.e / con.h  # Frequency of gamma transition in Hz
    finDecayHL = 117e-12  # Final decay half-life in seconds
    finDecayRate = finDecayHL / np.log(2)  # Final decay parameter in s 
    IC_Coefficient = 0.773 # Internal conversion coefficient for 177mLu
    I_a = 7/2  # Spin of the ground state for 177mLu
    I_b = 9/2  # Spin of the excited state for 177mLu
    eta_Lu = 0.0629 # Decay fraction at desired energy (121.621 keV)
    # Energy Uncertainty of Transition
    dE_Lu = con.hbar / (finDecayRate*con.e) # Energy uncertainty in eV
    detectorThickness = 0.1 # 1mm listed in cm
    detector = con.pi*(4**2)*detectorThickness # Detector Size (0.5mm thick and 8cm in diameter) cm^3

    debye = {
        'LuAl2O3': 950,
        'LuAG': 750, 
        'Lu2O3': 400
    }  # Debye temperatures in K for different materials

    activity = {
        'LuAl2O3': 2.87e6,
        'LuAG':8.7e8, 
        'Lu2O3': 1.74e9
    }

    n177Lu = { # Calculated using transmutation code following detector size listed here
        'LuAl2O3': 3.17e16,
        'LuAG': 9.62e18, 
        'Lu2O3': 1.93e19
    }
    ###################################################################
    #### Plot the Recoilless Fraction as a Function of Temperature ####
    ###################################################################
    temps = np.linspace(4, 300, num=100)  # Temperature range
    frac = {
        'LuAG': [],
        'LuAl2O3': [],
        'Lu2O3': []
    }
    for material, debeye_temp in debye.items():
        frac[material] = [recoillessFraction(mLu, Egamma, debeye_temp, T) for T in temps] # Recoilless fraction for each temperature
    plt.figure(figsize=(10, 6))
    for material in debye:
        plt.plot(temps,frac[material], label = material)
        print("Recoilless Fraction for ",material," at 300K:", round(recoillessFraction(mLu, Egamma, debye[material], 300),4))
    plt.xlabel('Temperature (K)')
    plt.ylabel('Recoilless Fraction')
    plt.title('Recoilless Fraction of Various Lutetium Bearing Lattices')
    plt.legend()
    plt.grid()
    plt.show()
    
    #####################################################################
    #### Plot the Absorption Cross-Section as a Function of Velocity ####
    #####################################################################
    vel = np.linspace(-0.05,0.05,num=1000)  # Velocity range in m/s
    detunes = [velocityDetuning(v, freqGamma) for v in vel]  # Detuning energies in eV
    maxCross = {}
    absCross = {}
    plt.figure(figsize=(10, 6))
    for material in debye:
        maxCross[material] = maximumAbsorptionCross(I_a, I_b, IC_Coefficient, con.c/freqGamma, recoillessFraction(mLu, Egamma, debye[material], 300))  # Maximum absorption cross-section
        absCross[material] = [absorptionCross(detune, dE_Lu, maxCross[material]) for detune in detunes]  # Absorption cross-section for each detuning (m^2)
        absCross[material] = np.array(absCross[material]) * 1e4  # Convert to numpy array for plotting and to cm^2
        print(f"{material} Max Absorption Cross-Section: {maxCross[material]:.2e} cm^2")
        plt.plot(vel,absCross[material], label = material)
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Absorption Cross-Section (cm^2)')
    plt.title('Absorption Cross-Section of Various Lutetium Lattices')
    plt.legend()
    plt.grid()
    plt.show()

    ###################################################
    #### Plot the Signal as a Function of Detuning ####
    ###################################################
    fluxFrac = 0.1  # Estimated fraction of flux at absorber
    emissions = {}
    signal = {}
    plt.figure(figsize=(10,6))
    for material in ['LuAl2O3', 'LuAG', 'Lu2O3']: #  'LuAl2O3', 'LuAG', 'Lu2O3'
        emissions[material] = activity[material]*eta_Lu*(1-IC_Coefficient)*fluxFrac
        print(f'{material} Count Rate: {emissions[material]:.2e}')
        signal[material] = emissions[material] * np.exp(-absCross[material]*n177Lu[material]*detectorThickness)
        plt.plot(vel,signal[material], label = material)
    plt.yscale('log')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Detected Count Rate (counts/s)')
    plt.title('Mossbauer Absorption Signal in Various Lutetium Lattices')
    # plt.ylim(1e5, 1e7)  # Set y-axis limits for better visibility
    plt.legend()
    plt.grid()
    plt.show() 

if __name__ == "__main__":
    main()