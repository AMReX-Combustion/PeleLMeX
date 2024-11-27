import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Line width for plots
line_width = 2.5  

# Get the data
file_dir = os.path.dirname(os.path.abspath(__file__))
file_name = os.path.join(file_dir, "temporals/tempExtremas")
col_names = ["time", "max_density", "max_rho.Y(AR)", "max_rho.Y(N2)", "max_rho.Y(CO2)"]
data = pd.read_csv(file_name, usecols=col_names)

# Nicer column labels for each species
var_names = ["AR", "N2", "CO2"]

# Molar masses (kg/mol)
molar_masses = {
    "AR":  0.040, 
    "N2":  0.028,  
    "CO2": 0.044  
}

# Get density and mass data for each species
data_rho = data["max_density"]
data_rhoY = data[["max_rho.Y(AR)", "max_rho.Y(N2)", "max_rho.Y(CO2)"]]
data_rhoY.columns = var_names

# Calculate mass fraction data and rename columns
data_Y = data_rhoY.div(data_rho, axis=0)  # mass fraction = rhoY / density
data_Y["Sum(Y_m)"] = data_Y.sum(axis=1)

# Calculate moles = mass / molar mass
data_moles = data_rhoY.div([molar_masses[comp] for comp in var_names], axis=1)

# Exact solutions
time = data['time']
exact_moles = pd.DataFrame({
    "AR": data_moles["AR"].iloc[0],
    "N2": data_moles["N2"].iloc[0],
    "CO2": np.maximum(0, (time - 0.5) / molar_masses["CO2"])
}, index=data.index)

# Calculate errors between data_moles and exact_moles
error_moles = np.abs(data_moles - exact_moles)

# Calculate maximum error for each species
max_error_moles = error_moles.max()

# Display the maximum error for each species
print("\n===========================================")
print("Maximum error in moles for each species:")
for species, error in max_error_moles.items():
    print(f"{species}: {error:.16e} moles")
print("===========================================\n")

# Create subplots (2 rows, 2 columns)
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
fig.subplots_adjust(hspace=0.4, wspace=0.3)

# Plot 1: Mass Fractions vs. Time
for var in var_names:
    axs[0, 0].plot(time, data_Y[var], label=var, linewidth=line_width)
axs[0, 0].plot(time, data_Y["Sum(Y_m)"], label="Sum(Y$_m$)", linestyle="--", color="black", linewidth=line_width)
axs[0, 0].set_xlabel('Time (s)')
axs[0, 0].set_ylabel('Mass Fraction')
axs[0, 0].legend()
axs[0, 0].grid(True)

# Plot 2: Density vs. Time
axs[0, 1].plot(time, data_rho, label="Density", color="black", linewidth=line_width)
axs[0, 1].set_xlabel('Time (s)')
axs[0, 1].set_ylabel('Density (kg/m$^3$)')
axs[0, 1].legend()
axs[0, 1].grid(True)

# Plot 3: Number of Moles vs. Time
for var in var_names:
    axs[1, 0].plot(time, data_moles[var], label=f'Moles of {var}', linewidth=line_width)
axs[1, 0].set_xlabel('Time (s)')
axs[1, 0].set_ylabel('Number of Moles')
axs[1, 0].legend()
axs[1, 0].grid(True)

# Plot 4: Error in Moles vs. Time
for var in var_names:
    axs[1, 1].plot(time, error_moles[var], label=f'Error in {var}', linewidth=line_width)
axs[1, 1].set_xlabel('Time (s)')
axs[1, 1].set_ylabel('Error in Number of Moles')
axs[1, 1].legend()
axs[1, 1].grid(True)

# Show the plots
plt.tight_layout()
plt.show()
