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
    mLu = 176.9437636 * con.atomic_mass # Lutetium mass in kg
    Egamma = 121.62e3  # Gamma energy in eV
    freqGamma = Egamma * con.e / con.h  # Frequency of gamma transition in Hz
    finDecayHL = 117e-12  # Final decay half-life in seconds
    finDecayRate = finDecayHL / np.log(2)  # Final decay parameter in s 
    decayCon177mLu = np.log(2)/(160.35*24*3600)  # Decay constant for 177mLu in s^-1
    ci177mLu = 0.002  # Initial activity in Ci
    activity177mLu = ci177mLu * 3.7e10  # Convert to Bq
    IC_Coefficient = 0.773 # Internal conversion coefficient for 177mLu
    I_a = 7/2  # Spin of the ground state for 177mLu
    I_b = 9/2  # Spin of the excited state for 177mLu
    eta_Lu = 0.0629 # Decay fraction at desired energy (121.621 keV)
    # Debye Temperatures
    debye = {'LuAG': 750, 'LuAl2O3': 950, 'Lu2O3': 400}  # Debye temperatures in K for different materials
    # Activity of 177mLu in Bq
    activity = {'LuAg':}

    debeye_LuAG = 750  # Debye temperature for LuAG in K
    # Energy Uncertainty of Transition
    dE_Lu = con.hbar / (finDecayRate*con.e) # Energy uncertainty in eV

    ###################################################################
    #### Plot the Recoilless Fraction as a Function of Temperature ####
    ###################################################################
    temps = np.linspace(4, 300, num=100)  # Temperature range
    frac = [recoillessFraction(mLu, Egamma, debeye_LuAG, T) for T in temps]  # Recoilless fraction for each temperature
    plt.figure(figsize=(10, 6))
    plt.plot(temps, frac, label='Recoilless Fraction')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Recoilless Fraction')
    plt.title('Recoilless Fraction of LuAG Lattice')
    # plt.legend()
    plt.grid()
    plt.show()
    print("Recoilless Fraction at 300K:", round(recoillessFraction(mLu, Egamma, debeye_LuAG, 300),4))


    #####################################################################
    #### Plot the Absorption Cross-Section as a Function of Velocity ####
    #####################################################################
    vel = np.linspace(-0.02,0.02,num=1000)  # Velocity range in m/s
    detunes = [velocityDetuning(v, freqGamma) for v in vel]  # Detuning energies in eV
    maxCross = maximumAbsorptionCross(I_a, I_b, IC_Coefficient, con.c/freqGamma, recoillessFraction(mLu, Egamma, debeye_LuAG, 300))  # Maximum absorption cross-section
    absCross = [absorptionCross(detune, dE_Lu, maxCross) for detune in detunes]  # Absorption cross-section for each detuning (m^2)
    absCross = np.array(absCross) * 1e4  # Convert to numpy array for plotting and to cm^2
    plt.figure(figsize=(10, 6))
    plt.plot(vel, absCross, label='Absorption Cross-Section')
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Absorption Cross-Section (cm^2)')
    plt.title('Absorption Cross-Section of LuAG Lattice')
    plt.legend()
    plt.grid()
    plt.show()
    print("Maximum Absorption Cross-Section:", round(maxCross*1e4, 4), "cm^2")

    ###################################################
    #### Plot the Signal as a Function of Detuning ####
    ###################################################
    solidAngle = 0.1  # Estimated solid angle in steradians
    emissions = activity177mLu * eta_Lu * IC_Coefficient * recoillessFraction(mLu, Egamma, debeye_LuAG, 300)

if __name__ == "__main__":
    main()