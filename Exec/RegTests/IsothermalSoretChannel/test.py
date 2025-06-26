import os
import numpy.testing as npt
import pandas as pd
import unittest

class SpeciesBalTestCase(unittest.TestCase):
    """Test species balance with isothermal walls and soret"""

    def test_composition(self):
        """Verify species conservation"""

        # Load the data
        fdir = os.path.abspath(".")
        fname = os.path.join(fdir, "temporals/tempSpecies")
        df = pd.read_csv(fname)
        print(df)
        for col in df.columns:
            if col.startswith('rhoYnew'):
                init_value = df[col][0]
                if init_value != 0.0:
                    print('testing (relative)', col)
                    npt.assert_allclose(df[col], init_value, rtol=1e-13)
                else:
                    print('testing (absolute)', col)
                    npt.assert_allclose(df[col], init_value, atol=1e-13)
            if col.startswith('netFlux'):
                print('testing (absolute)', col)
                npt.assert_allclose(df[col], 0.0, atol=1e-13)

if __name__ == "__main__":
    unittest.main()
