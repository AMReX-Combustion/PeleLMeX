import os
import numpy.testing as npt
import pandas as pd
import unittest

def read_last_dataset(path, **read_csv_kwargs):
    """Read only the most-recently-appended block from a header-repeating CSV.

    The file is assumed to have the same column-header line at the start
    of every appended block. This finds the last occurrence of that line
    and returns a DataFrame of just the rows after it, with the header
    taken from that line. Extra kwargs are forwarded to pandas.read_csv.
    """
    with open(path) as f:
        header = f.readline().rstrip('\n')
        last_header_lineno = 0  # 0 = the very first header at the top
        for lineno, line in enumerate(f, start=1):
            if line.rstrip('\n') == header:
                last_header_lineno = lineno
    return pd.read_csv(path, skiprows=last_header_lineno, **read_csv_kwargs)

class SpeciesBalTestCase(unittest.TestCase):
    """Test species balance with isothermal walls and soret"""

    def test_composition(self):
        """Verify species conservation"""

        # Load the data
        fdir = os.path.abspath(".")
        fname = os.path.join(fdir, "temporals/tempSpecies")
        df = read_last_dataset(fname)
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
