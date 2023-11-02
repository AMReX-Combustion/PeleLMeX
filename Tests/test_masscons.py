# ========================================================================
#
# Imports
#
# ========================================================================
import os
import numpy.testing as npt
import pandas as pd
import unittest


# ========================================================================
#
# Test definitions
#
# ========================================================================
class ConsTestCase(unittest.TestCase):
    """Tests for conservation in Pele."""

    def test_conservation(self):
        """Is mass conserved?"""

        # Load the data
        fdir = os.path.abspath(".")
        fname = os.path.join(fdir, "temporals/tempMass")
        cnames=['step','time','mass','dmdt','massFluxBal','balance']
        df = pd.read_csv(fname, names=cnames, delim_whitespace=True)
        init_mass = df.mass[0]
        npt.assert_allclose(df.mass, init_mass, rtol=1e-13)


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":
    unittest.main()
