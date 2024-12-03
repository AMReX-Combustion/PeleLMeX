import os
import re
import numpy as np
import pandas as pd
import unittest

class CompTestCase(unittest.TestCase):
    """Test composition of species with external sources"""

    def test_composition(self):
        """Are the number of moles, mass fractions, and density correct?"""

        # Molar masses (kg/mol)
        molar_masses = {"AR": 0.040, "N2": 0.028, "CO2": 0.044}

        # Load the data
        file_dir = os.path.abspath(".")
        file_name = os.path.join(file_dir, "temporals/tempExtremas")
        col_names = ["time", "max_density", "max_rho.Y(AR)", "max_rho.Y(N2)", "max_rho.Y(CO2)"]
        var_names = ["AR", "N2", "CO2"]
        data = pd.read_csv(file_name, usecols=col_names, delimiter=',')
        print("Raw data:\n", data)

        # Ensure numeric types for all data columns
        data = data.apply(pd.to_numeric, errors="coerce")

        # Check for NaN values and raise an error if found
        if data.isnull().any().any():
            raise ValueError("Data contains NaN values after conversion.")
        
        time = data["time"]

        # Parse input file for necessary solution parameters
        input_file = os.path.join(file_dir, "composition-test-2d.inp")
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

        # Ensure numeric types for all data columns
        data_rho = pd.to_numeric(data_rho, errors="coerce")
        data_rhoY = data_rhoY.apply(pd.to_numeric, errors="coerce")

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
            key: exact_rhoY[key] / molar_masses[key] for key in var_names}, 
            index=data.index)

        # Error calculations
        max_error_density = np.abs(data_rho - exact_density).max()
        max_error_Y = np.abs(data_Y.iloc[:, :-1] - exact_Y).max()
        max_error_moles = np.abs(data_moles - exact_moles).max()

        # Expected max errors ["AR", "N2", "CO2"]
        expected_errors_moles = np.array([1e-15, 1e-15, 1.0])
        expected_errors_massfrac = np.array([1.032e-02, 4.3e-03, 1.47e-02])
        expected_error_rho = 4e-02

        # Perform tests using np.testing to make sure errors haven't increased
        np.testing.assert_array_less(
            np.array([max_error_moles[var] for var in var_names]),
            expected_errors_moles,
            err_msg="Maximum mole errors exceed specified tolerance."
        )

        np.testing.assert_array_less(
            np.array([max_error_Y[var] for var in var_names]),
            expected_errors_massfrac,
            err_msg="Maximum mass fraction errors exceed specified tolerance."
        )

        np.testing.assert_array_less(
            np.array([max_error_density]),
            np.array([expected_error_rho]),
            err_msg="Maximum density error exceeds specified tolerance."
        )

if __name__ == "__main__":
    unittest.main()
