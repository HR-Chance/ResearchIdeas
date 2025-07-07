# LuAluminaSamplePrep.py
# This file is intended for calculating the precursors for producing Lu-Alumina samples
# Reference paper for the procedure: "Structural Effects of Lanthanide Dopants on Alumina"
# DOI: https://doi.org/10.1038/srep39946

# Define the molecular weights (g/mol) as global variables
lu_mass = 174.9668      # Lutetium
al_mass = 26.981538      # Aluminum
no3_mass = 62.0049       # Nitrate
h2o_mass = 18.01528      # Water
h_mass = 1.00784         # Hydrogen
o_mass = 15.999          # Oxygen
mg_mass = 24.305         # Magnesium
# Calculate the molar masses of the precursors
mgo_mass = mg_mass + o_mass  # Molar mass of Magnesium Oxide in g/mol
aluminumNitrateNonahydrate = (al_mass + no3_mass * 3 + h2o_mass * 9)  # g/mol
magnesiumNitrateHexahydrate = (mg_mass + no3_mass * 2 + h2o_mass * 6)  # g/mol
lutetiumNitrate = (lu_mass + no3_mass * 3)  # g/mol


def solutionA(X_Lu, X_Al, X_MgO, X_O, Fin_Mass):
    global lu_mass, al_mass, o_mass, mgo_mass, aluminumNitrateNonahydrate, magnesiumNitrateHexahydrate, lutetiumNitrate
    luAlO = lu_mass*X_Lu + al_mass*X_Al + o_mass*X_O  # Molar mass of Lu-Alumina in g/mol
    # Calculate the moles of each precursor needed
    molesLu = (X_Lu * Fin_Mass) / luAlO  # Moles of Lutetium
    molesAl = (X_Al * Fin_Mass) / luAlO  # Moles of Aluminum
    molesMgO = (X_MgO * Fin_Mass) / mgo_mass  # Moles of Magnesium Oxide
    # Calculate the mass of each precursor needed
    massLuNitrate = molesLu * lutetiumNitrate  # Mass of Lutetium Nitrate in grams
    massAlNitrate = molesAl * aluminumNitrateNonahydrate  # Mass of Aluminum Nitrate in grams
    massMg = molesMgO * magnesiumNitrateHexahydrate  # Mass of Magnesium Nitrate in grams
    return massLuNitrate, massAlNitrate, massMg

def SolutionB(X_Carb, X_Hydroxide, X_H2O, mlH2O):
    mTotal = mlH2O/X_H2O  # Total mass of the solution in grams
    mCarbonate = mTotal * X_Carb  # Mass of Ammonium Bicarbonate in grams
    mHydroxide = 4*(mTotal * X_Hydroxide)  # Mass of Ammonium Hydroxide in grams (25% solution)
    mDI = 500 - 0.75*mHydroxide  # Mass of Water in grams
    return mCarbonate, mHydroxide, mDI

def main():
    global lu_mass, al_mass, no3_mass, h2o_mass, h_mass, o_mass, mg_mass
    # Solution A Goal Parameters
    xLu = 0.002 # Desired mole fraction of Lutetium in the final sample
    xAl = 1.998 # Desired mole fraction of Aluminum in the final sample
    xO = 3 # Desired mole fraction of Oxygen in the final sample
    xMgO = 0.00025 # Desired mole fraction of Magnesium-Oxide in the final sample (250 ppm)
    finMass = 25  # Final mass of Lu-Alumina sample in grams
    massesSolnA = solutionA(xLu, xAl, xMgO, xO, finMass)  # Calculate the masses of precursors for solution A

    # Solution B Goal Parameters
    wtCarbonate = 0.11 # Desired weight fraction of Ammonium Bicarbonate in the solution
    wtHydroxide = 0.03 # Desired weight fraction of Ammonium Hydroxide in the solution
    wtH2O = 0.86 # Desired weight fraction of Water in the solution
    mlH2O = 500 # Total volume of water in mL
    massesSolnB = SolutionB(wtCarbonate, wtHydroxide, wtH2O, mlH2O)  # Calculate the masses of precursors for solution B

    # Print the results
    print("Masses of precursors for " + str(finMass) + "g of Lu-Alumina:")
    print("Solution A Parameters:")
    print(f"Lutetium Nitrate [Lu(NO3)3·xH2O]: {massesSolnA[0]:.4f} g")
    print(f"Aluminum Nitrate [Al(NO3)3·9H2O]: {massesSolnA[1]:.4f} g")
    print(f"Magnesium Nitrate [Mg(NO3)2·6H2O]: {massesSolnA[2]:.4f} g")

    print("\nSolution B Parameters:")
    print(f"Powdered Ammonium Bicarbonate [NH4HCO3]: {massesSolnB[0]:.4f} g")
    print(f"25% Ammonium Hydroxide [NH4OH]: {massesSolnB[1]:.4f} g")
    print(f"Deionized Water [H2O]: {massesSolnB[2]:.4f} g")

if __name__ == "__main__":
    main()