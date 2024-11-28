import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# File paths
file_dir = os.path.dirname(os.path.abspath(__file__))
data_file = os.path.join(file_dir, "temporals/tempExtremas")
input_file = os.path.join(file_dir, "composition-test-2d.inp")

# Load data
col_names = ["time", "max_density", "max_rho.Y(AR)", "max_rho.Y(N2)", "max_rho.Y(CO2)"]
data = pd.read_csv(data_file, usecols=col_names, delimiter=',')
var_names = ["AR", "N2", "CO2"]
time = data['time']

# Constants
molar_masses = {"AR": 0.040, "N2": 0.028, "CO2": 0.044}
line_width = 2.5

# Parse input file for necessary solution parameters
with open(input_file, 'r') as file:
    content = file.read()
    soln_params = {
        "ts": float(re.search(r"prob\.extRhoYCO2_ts\s*=\s*([\d\.]+)", content).group(1)),
        "src_strength": float(re.search(r"prob\.extRhoYCO2\s*=\s*([\d\.]+)", content).group(1))
    }
ts, k = soln_params["ts"], soln_params["src_strength"]

# Density and species mass data
data_rho = data["max_density"]
data_rhoY = data.iloc[:, 2:].copy()
data_rhoY.columns = var_names

# Initial values
AR_0, N2_0, CO2_0 = data_rhoY.iloc[0].to_dict().values()
rho_0 = data_rho.iloc[0]

# Mass fractions and number of moles
data_Y = data_rhoY.div(data_rho, axis=0)
data_Y["Sum(Y_m)"] = data_Y.sum(axis=1)
data_moles = data_rhoY.div([molar_masses[comp] for comp in var_names], axis=1)

# Exact solutions
exact_rhoY = {
    "AR": AR_0,
    "N2": N2_0,
    "CO2": CO2_0 * (1 - np.heaviside(time - ts, 0)) 
           + CO2_0 * np.exp(k * (time - ts)) * np.heaviside(time - ts, 0)
}
exact_density = rho_0 + CO2_0 * (np.exp(k * (time - ts)) - 1) * np.heaviside(time - ts, 0)
exact_Y = pd.DataFrame({
    key: exact_rhoY[key] / exact_density for key in var_names}, 
    index=data.index)
exact_moles = pd.DataFrame({
    key: exact_rhoY[key] / molar_masses[key] for key in var_names}
    , index=data.index)

# Error calculations
max_error_density = np.abs(data_rho - exact_density).max()
max_error_Y = np.abs(data_Y.iloc[:, :-1] - exact_Y).max()
max_error_moles = np.abs(data_moles - exact_moles).max()

# Print errors
print("\n=====================================================")
print("Max absolute errors:")
print(f"rho: {max_error_density:.16e} kg/m^3")
for species in max_error_Y.index:
    print(f"Y({species}): {max_error_Y[species]:.16e}")
for species in max_error_moles.index:
    print(f"{species}: {max_error_moles[species]:.16e} moles")
print("=====================================================\n")

# Plot results
fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))
xticks = np.linspace(time.min(), time.max(), 5)

# Plot 1: Mass fractions
for var in var_names:
    axs[0].plot(time, data_Y[var], label=f"Y({var})", linewidth=line_width)
axs[0].plot(time, data_Y["Sum(Y_m)"], label="Sum(Y$_m$)", linestyle="--", color="black", linewidth=line_width)
axs[0].set_prop_cycle(None)
for i, var in enumerate(var_names):
    axs[0].plot(time, exact_Y[var], 'o', markerfacecolor='none', label='Exact solutions' if i == 0 else None)

axs[0].set_xticks(xticks)
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('Mass Fraction')
axs[0].legend(loc='upper left')
axs[0].grid(True)

# Plot 2: Density
axs[1].plot(time, data_rho, label="Density", color="black", linewidth=line_width)
axs[1].plot(time, exact_density, 'ko', markerfacecolor='none', label='Exact solution')
axs[1].set_xticks(xticks)
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Density (kg/m$^3$)')
axs[1].legend(loc='upper left')
axs[1].grid(True)

# Plot 3: Moles
for var in var_names:
    axs[2].plot(time, data_moles[var], label=f"Moles of {var}", linewidth=line_width)
axs[2].set_prop_cycle(None)
for i, var in enumerate(var_names):
    axs[2].plot(time, exact_moles[var], 'o', markerfacecolor='none', label='Exact solutions' if i == 0 else None)
axs[2].set_xticks(xticks)
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Number of Moles')
axs[2].legend(loc='upper left')
axs[2].grid(True)

plt.tight_layout()
plt.show()
