import os
import numpy as np
import pandas as pd
import unittest

class CompTestCase(unittest.TestCase):
    """Tests composition of species with external sources"""

    def test_composition(self):
        """Are the number of moles for each species correct?"""

        # Molar masses (kg/mol)
        molar_masses = {"AR":  0.040, "N2":  0.028, "CO2": 0.044}

        # Load the data
        file_dir = os.path.dirname(os.path.abspath(__file__))
        file_name = os.path.join(file_dir, "temporals/tempExtremas")
        col_names = ["time", "max_rho.Y(AR)", "max_rho.Y(N2)", "max_rho.Y(CO2)"]
        var_names = ["AR", "N2", "CO2"]
        data = pd.read_csv(file_name, usecols=col_names)
        
        # Calculate moles = mass / molar mass
        data_moles = data[["max_rho.Y(AR)", "max_rho.Y(N2)", "max_rho.Y(CO2)"]]
        data_moles.columns = var_names
        data_moles = data_moles[var_names].div([molar_masses[comp] for comp in var_names], axis=1)

        # Exact solutions
        exact_moles = pd.DataFrame({
            "AR": data_moles["AR"].iloc[0],
            "N2": data_moles["N2"].iloc[0],
            "CO2": np.maximum(0, (data['time'] - 0.5) / molar_masses["CO2"])
        }, index=data.index)

        # Calculate errors and maximum error
        max_error_moles = np.abs(data_moles - exact_moles).max()

        # Expected max errors using forward Euler
        expected_errors = {
            "AR": 2.5e-10, 
            "N2": 3.6e-11, 
            "CO2": 2.3e-02
        }

        # Convert arrays for np.testing comparison
        expected_errors = np.array([expected_errors[species] for species in var_names])
        calculated_errors = np.array([max_error_moles[species] for species in var_names])

        print("Testing")
        np.testing.assert_array_less(
            calculated_errors,
            expected_errors,
            err_msg="Maximum errors are not less than specified tolerance."
        )

if __name__ == "__main__":
    unittest.main()